"""
Stage 3b - Per-pixel soil prediction maps over Beni Suef AOI

Strategy (locked after Stage 3a v2):
  * Mean map         = IDW(q50) on the 60 profiles directly   (best CV RMSE)
  * Uncertainty band = QRF q05/q95 with 50 trimmed covariates (calibrated PICP)
  * Cross-check      = |IDW_q50 - QRF_q50| / max(QRF q95-q05, eps)

For each (property, depth) we write 4 single-band GeoTIFFs:
    {prop}_{depth}_q50_idw.tif
    {prop}_{depth}_q05_qrf.tif
    {prop}_{depth}_q50_qrf.tif
    {prop}_{depth}_q95_qrf.tif

Inputs  : outputs/stage2b_local_stack/covariate_stack_30m.tif
          outputs/stage1/soil_profiles_0_{30,100}cm.gpkg
          outputs/stage3a/predictors_used.txt
Outputs : outputs/stage3b/
"""
from __future__ import annotations
from pathlib import Path
import warnings
import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio
from rasterio.windows import Window
from rasterio.features import geometry_mask
from scipy.spatial import cKDTree
from quantile_forest import RandomForestQuantileRegressor

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

ROOT      = Path(__file__).resolve().parents[1]
STAGE1    = ROOT / "outputs" / "stage1"
STAGE2B   = ROOT / "outputs" / "stage2b_local_stack"
STAGE3A   = ROOT / "outputs" / "stage3a"
OUT       = ROOT / "outputs" / "stage3b"
OUT.mkdir(parents=True, exist_ok=True)

STACK_TIF  = STAGE2B / "covariate_stack_30m.tif"
AOI_PATH   = ROOT / "data" / "AOI_BeniSuef.gpkg"   # if present
RANDOM_STATE = 42
N_TREES = 500
CHUNK = 512   # window size

TARGETS = ["EC", "ESP", "CaCO3", "pH", "CEC", "Clay", "SOM",
           "Ks", "AW", "N_avail", "P_avail", "K_avail", "Sand_total"]
DEPTHS = {
    "0_100cm": STAGE1 / "soil_profiles_0_100cm.gpkg",
    "0_30cm":  STAGE1 / "soil_profiles_0_30cm.gpkg",
}

# ---------------------------------------------------------------------------
def load_profiles(depth_key, predictor_cols):
    soil = gpd.read_file(DEPTHS[depth_key], engine="pyogrio").to_crs(32636)
    soil["Sand_total"] = soil["Coarse_sand"] + soil["Fine_sand"]
    soil["X_utm"] = soil.geometry.x.values
    soil["Y_utm"] = soil.geometry.y.values
    cov = pd.read_csv(ROOT / "outputs" / "stage2" / "covariates_at_profiles.csv")
    df = soil.drop(columns="geometry").merge(cov, on="Profile", how="inner")
    return df

def idw_window(xy_train, y_train, xy_test, power=2, k=8):
    tree = cKDTree(xy_train)
    dist, idx = tree.query(xy_test, k=min(k, len(xy_train)))
    if dist.ndim == 1:
        dist = dist[:, None]; idx = idx[:, None]
    dist = np.where(dist == 0, 1e-9, dist)
    w = 1.0 / dist ** power
    return (w * y_train[idx]).sum(axis=1) / w.sum(axis=1)

# ---------------------------------------------------------------------------
def main():
    if not STACK_TIF.exists():
        raise SystemExit(f"Missing covariate stack: {STACK_TIF}\n"
                         "Run stage2b_download_stack.py first.")

    predictor_cols = (STAGE3A / "predictors_used.txt").read_text().splitlines()
    qrf_bands = [c for c in predictor_cols if c not in ("X_utm", "Y_utm")]

    with rasterio.open(STACK_TIF) as src:
        stack_bands = list(src.descriptions)
        if not all(b for b in stack_bands):
            # fall back to band_names.txt if descriptions not embedded
            stack_bands = (STAGE2B / "band_names.txt").read_text().splitlines()
        band_idx = {n: i + 1 for i, n in enumerate(stack_bands)}
        missing = [b for b in qrf_bands if b not in band_idx]
        if missing:
            raise SystemExit(f"Stack missing bands: {missing[:5]}...")

        height, width = src.height, src.width
        transform = src.transform
        crs = src.crs
        print(f"Stack: {width} x {height}, CRS={crs}, bands={src.count}")

        # AOI mask (optional)
        aoi_mask = None
        if AOI_PATH.exists():
            aoi = gpd.read_file(AOI_PATH, engine="pyogrio").to_crs(crs)
            aoi_mask = geometry_mask([g for g in aoi.geometry],
                                     out_shape=(height, width),
                                     transform=transform, invert=True)

        # Pre-compute pixel-centre coordinates per row (X same, Y per row)
        cols = np.arange(width)
        x_centres = transform.c + (cols + 0.5) * transform.a  # transform.a = pixel width
        # row centre: y = transform.f + (row + 0.5) * transform.e (e is negative)

        out_meta = dict(driver="GTiff", height=height, width=width, count=1,
                        dtype="float32", crs=crs, transform=transform,
                        nodata=np.float32(np.nan), compress="deflate",
                        predictor=3, tiled=True, blockxsize=512, blockysize=512,
                        BIGTIFF="IF_NEEDED")

        for depth_key in DEPTHS:
            df = load_profiles(depth_key, predictor_cols)
            print(f"\n=== Depth {depth_key} (n={len(df)}) ===")

            for target in TARGETS:
                sub = df.dropna(subset=[target] + predictor_cols)
                if len(sub) < 30:
                    print(f"  skip {target}: n<30"); continue
                print(f"  -> {target} (n={len(sub)})", flush=True)

                X_train = sub[predictor_cols].to_numpy()
                y_train = sub[target].to_numpy()
                xy_train = sub[["X_utm", "Y_utm"]].to_numpy()

                qrf = RandomForestQuantileRegressor(
                    n_estimators=N_TREES, min_samples_leaf=3,
                    max_features="sqrt", random_state=RANDOM_STATE, n_jobs=-1,
                ).fit(X_train, y_train)

                paths = {q: OUT / f"{target}_{depth_key}_{q}.tif"
                         for q in ("q50_idw", "q05_qrf", "q50_qrf", "q95_qrf")}
                writers = {q: rasterio.open(p, "w", **out_meta)
                           for q, p in paths.items()}

                try:
                    for row0 in range(0, height, CHUNK):
                        rh = min(CHUNK, height - row0)
                        for col0 in range(0, width, CHUNK):
                            cw = min(CHUNK, width - col0)
                            win = Window(col0, row0, cw, rh)

                            # read only QRF bands
                            arr = np.stack([src.read(band_idx[b], window=win)
                                            for b in qrf_bands], axis=0)  # (B,h,w)
                            valid = np.all(np.isfinite(arr), axis=0)

                            if aoi_mask is not None:
                                valid &= aoi_mask[row0:row0+rh, col0:col0+cw]

                            out = {q: np.full((rh, cw), np.nan, np.float32)
                                   for q in paths}
                            if valid.any():
                                # build feature matrix in qrf_bands order, then add X,Y
                                X = arr.reshape(arr.shape[0], -1).T  # (N,B)
                                xs = np.broadcast_to(x_centres[col0:col0+cw], (rh, cw))
                                ys = transform.f + (row0 + np.arange(rh) + 0.5) * transform.e
                                ys = np.broadcast_to(ys[:, None], (rh, cw))
                                X = np.column_stack([X, xs.ravel(), ys.ravel()])
                                # reorder to predictor_cols
                                # currently columns = qrf_bands + [X_utm, Y_utm]
                                # predictor_cols == keep_bands + [X_utm, Y_utm]
                                # qrf_bands == keep_bands  -> already aligned
                                v = valid.ravel()
                                if v.any():
                                    q = qrf.predict(X[v], quantiles=[0.05, 0.50, 0.95])
                                    idw = idw_window(xy_train, y_train,
                                                     X[v][:, -2:])
                                    flat = {k: np.full(v.size, np.nan, np.float32)
                                            for k in paths}
                                    flat["q05_qrf"][v] = q[:, 0]
                                    flat["q50_qrf"][v] = q[:, 1]
                                    flat["q95_qrf"][v] = q[:, 2]
                                    flat["q50_idw"][v] = idw
                                    for k in paths:
                                        out[k] = flat[k].reshape(rh, cw)

                            for k, w in writers.items():
                                w.write(out[k], 1, window=win)
                finally:
                    for w in writers.values():
                        w.close()

    print(f"\nWrote: {OUT}")

if __name__ == "__main__":
    main()
