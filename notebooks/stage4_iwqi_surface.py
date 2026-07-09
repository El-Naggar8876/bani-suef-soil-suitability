"""
Stage 4 - Irrigation Water Quality Index (IWQI) surface for Beni Suef AOI

Inputs:
    outputs/stage1/water_samples_iwqi.gpkg  (20 samples, IWQI_meireles)
    outputs/stage2b_local_stack/covariate_stack_30m.tif (target grid)

Method:
    Ordinary Kriging (spherical variogram) on IWQI, fitted on a coarse
    prediction grid (250 m for speed) then NN-upsampled to the 30 m
    target grid. Uncertainty intervals from the OK kriging variance:
        q05 = z50 - 1.645 * sqrt(var)
        q95 = z50 + 1.645 * sqrt(var)
    (Bootstrap was tried but n=20 with near-co-located samples produces
    ill-conditioned variogram matrices and inflated quantiles.)

Outputs in outputs/stage4/:
    iwqi_q50.tif, iwqi_q05.tif, iwqi_q95.tif
    iwqi_kriging_var.tif    (point OK variance, single fit)
    iwqi_class.tif          (Meireles 5-class on q50)
    stage4_summary.txt
"""
from __future__ import annotations
from pathlib import Path
import warnings
import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio
from rasterio.transform import Affine
from pykrige.ok import OrdinaryKriging

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

ROOT     = Path(__file__).resolve().parents[1]
STAGE1   = ROOT / "outputs" / "stage1"
STAGE2B  = ROOT / "outputs" / "stage2b_local_stack"
OUT      = ROOT / "outputs" / "stage4"
OUT.mkdir(parents=True, exist_ok=True)

WATER    = STAGE1 / "water_samples_iwqi.gpkg"
STACK    = STAGE2B / "covariate_stack_30m.tif"

RANDOM_STATE = 42
N_BOOT = 200
COARSE_RES = 250.0   # metres for prediction grid (NN-upsampled to 30 m)

# Meireles (2010) IWQI classes
IWQI_CLASSES = [
    (85, 100, 1, "No restriction"),
    (70, 85,  2, "Low restriction"),
    (55, 70,  3, "Moderate restriction"),
    (40, 55,  4, "High restriction"),
    (0,  40,  5, "Severe restriction"),
]

def classify(iwqi):
    out = np.full(iwqi.shape, np.uint8(0))
    for lo, hi, cls, _ in IWQI_CLASSES:
        m = (iwqi >= lo) & (iwqi < hi if cls != 1 else iwqi <= hi)
        out[m] = cls
    return out

def fit_predict_ok(xs, ys, vals, gx, gy, model="spherical"):
    ok = OrdinaryKriging(xs, ys, vals, variogram_model=model,
                         verbose=False, enable_plotting=False, nlags=8)
    z, ss = ok.execute("grid", gx, gy)
    return np.asarray(z, float), np.asarray(ss, float)

# ---------------------------------------------------------------------------
def main():
    samples = gpd.read_file(WATER, engine="pyogrio").to_crs(32636)
    iwqi = samples["IWQI_meireles"].astype(float).to_numpy()
    sx = samples.geometry.x.to_numpy()
    sy = samples.geometry.y.to_numpy()
    print(f"Water samples: n={len(samples)}, IWQI range "
          f"{iwqi.min():.1f}-{iwqi.max():.1f}, mean {iwqi.mean():.1f}")

    # --- target 30 m grid from stack ---
    with rasterio.open(STACK) as src:
        H, W = src.height, src.width
        T = src.transform
        crs = src.crs
        finite_mask = np.isfinite(src.read(1))   # AOI footprint
        # 30 m pixel-centre coords
        x30 = T.c + (np.arange(W) + 0.5) * T.a
        y30 = T.f + (np.arange(H) + 0.5) * T.e   # T.e is negative

    # --- coarse 250 m grid ---
    px = abs(T.a)
    step = max(1, int(round(COARSE_RES / px)))
    cx = T.c + (np.arange(0, W, step) + step/2) * T.a
    cy = T.f + (np.arange(0, H, step) + step/2) * T.e
    cy_sorted = np.sort(cy)   # pykrige expects ascending
    print(f"Coarse grid: {len(cx)} x {len(cy)} (step {step} px = {step*px:.0f} m)")

    # --- single OK fit (mean + variance) ---
    z_coarse, var_coarse = fit_predict_ok(sx, sy, iwqi, cx, cy_sorted)
    # pykrige z_coarse: rows correspond to cy_sorted (ascending). Our raster
    # rows go top-to-bottom (descending y). Flip.
    z_coarse = np.flipud(z_coarse)
    var_coarse = np.flipud(var_coarse)

    # NN upsample to 30 m
    row_idx = np.clip(((y30 - cy[0]) / (T.e * step)).astype(int), 0, len(cy)-1)
    col_idx = np.clip(((x30 - cx[0]) / (T.a * step)).astype(int), 0, len(cx)-1)
    z50 = z_coarse[row_idx[:, None], col_idx[None, :]].astype(np.float32)
    var = var_coarse[row_idx[:, None], col_idx[None, :]].astype(np.float32)
    z50[~finite_mask] = np.nan
    var[~finite_mask] = np.nan

    # --- bootstrap for q05/q95 ---
    # Use OK kriging variance directly (standard geostatistical PI):
    # q05/q95 = z +/- 1.645 * sigma. Clip to plausible IWQI range [0, 100].
    sigma_coarse = np.sqrt(np.clip(var_coarse, 0, None))
    z05_c = z_coarse - 1.645 * sigma_coarse
    z95_c = z_coarse + 1.645 * sigma_coarse
    q05 = z05_c[row_idx[:, None], col_idx[None, :]].astype(np.float32)
    q95 = z95_c[row_idx[:, None], col_idx[None, :]].astype(np.float32)
    q05 = np.clip(q05, 0, 100)
    q95 = np.clip(q95, 0, 100)
    q05[~finite_mask] = np.nan
    q95[~finite_mask] = np.nan

    cls = classify(z50).astype(np.uint8)
    cls[~finite_mask] = 0   # 0 = nodata

    # --- write rasters ---
    meta = dict(driver="GTiff", height=H, width=W, count=1, crs=crs,
                transform=T, compress="deflate", predictor=3, tiled=True,
                blockxsize=512, blockysize=512, BIGTIFF="IF_NEEDED")
    fmeta = dict(meta, dtype="float32", nodata=np.float32(np.nan))
    cmeta = dict(meta, dtype="uint8", nodata=0, predictor=2)

    rasters = {
        "iwqi_q50.tif": (z50, fmeta),
        "iwqi_q05.tif": (q05, fmeta),
        "iwqi_q95.tif": (q95, fmeta),
        "iwqi_kriging_var.tif": (var, fmeta),
        "iwqi_class.tif": (cls, cmeta),
    }
    for name, (arr, m) in rasters.items():
        with rasterio.open(OUT / name, "w", **m) as dst:
            dst.write(arr, 1)
        print(f"  wrote {name}")

    # --- summary ---
    valid = finite_mask
    cls_counts = {c: int(((cls == c) & valid).sum()) for c, *_ in
                  [(c, *rest) for _, _, c, *rest in IWQI_CLASSES]}
    total = int(valid.sum())
    iwqi_pct = {c: 100 * n / total for c, n in cls_counts.items() if total}

    lines = [
        "Stage 4 - IWQI surface (Ordinary Kriging + bootstrap)",
        "=" * 60,
        f"Samples         : {len(samples)}",
        f"IWQI obs range  : {iwqi.min():.2f} - {iwqi.max():.2f} (mean {iwqi.mean():.2f})",
        f"Variogram       : spherical, nlags=8",
        f"PI method       : z +/- 1.645 * sqrt(OK_variance), clipped [0,100]",
        f"Coarse grid res : {step * px:.0f} m",
        f"Target grid     : {W} x {H} @ {px:.0f} m, EPSG:32636",
        f"Valid AOI px    : {total}",
        "",
        "Class shares on q50 surface:",
    ]
    for lo, hi, c, name in IWQI_CLASSES:
        lines.append(f"  {c} {name:22s} ({lo:>3}-{hi:>3}): "
                     f"{iwqi_pct.get(c, 0):5.1f} %")
    lines.append(f"\nMean uncertainty width (q95-q05): "
                 f"{np.nanmean(q95 - q05):.2f} IWQI units")
    text = "\n".join(lines)
    print("\n" + text)
    (OUT / "stage4_summary.txt").write_text(text)

if __name__ == "__main__":
    main()
