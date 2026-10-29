"""
Stage 12 - Spatial structure of classification uncertainty and its drivers.

Purpose (Reviewer 2, comment 4.2): the submitted manuscript reports aggregate
uncertainty but does not interpret where the classification is least confident,
or which properties drive that pattern.

Metrics per pixel:
  PIW      prediction-interval width of the index (q95 - q50 - q05 envelope)
  NCLS     number of suitability classes spanned by [q05, q95]
  AMBIG    binary-entropy of P(class >= S2); maximal where P = 0.5
  DBREAK   distance of the index q50 to the nearest class break
Candidate drivers:
  DPROF    Euclidean distance to the nearest of the 60 soil profiles
  KSSD     local standard deviation of predicted Ks (5 x 5 window) - soil-unit
           transition proxy
  KS       predicted Ks itself
  CLAYSD   local standard deviation of predicted clay (5 x 5 window)
  DEDGE    distance to the edge of the valid mask

Outputs -> 02_New_Analysis/outputs_rev/
"""
from __future__ import annotations
import sys, pathlib
import numpy as np, pandas as pd, rasterio
from scipy import ndimage
from scipy.stats import spearmanr

import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parent))
from revpaths import REPO as _REPO, OUT as _OUT, WAPOR as _WAPOR

HERE = pathlib.Path(__file__).resolve().parent
REPO = _REPO
sys.path.insert(0, str(REPO / "notebooks"))
OUT = _OUT
from stage5_ales_montecarlo import classify, STAGE3B, CLASS_BREAKS

S5 = REPO / "outputs" / "stage5"
BREAKS = np.array([12.5, 25.0, 50.0, 75.0])


def rd(p):
    with rasterio.open(p) as s:
        return s.read(1).astype("float32"), s.profile


def local_sd(a, valid, size=5):
    x = np.where(valid, a, np.nan).astype("float64")
    m = valid.astype("float64")
    xf = np.nan_to_num(x)
    s1 = ndimage.uniform_filter(xf, size) * size * size
    s2 = ndimage.uniform_filter(xf ** 2, size) * size * size
    cnt = ndimage.uniform_filter(m, size) * size * size
    cnt = np.where(cnt > 1, cnt, np.nan)
    var = s2 / cnt - (s1 / cnt) ** 2
    return np.sqrt(np.clip(var, 0, None)).astype("float32")


def main():
    # profile locations -> distance surface
    import geopandas as gpd
    gdf = gpd.read_file(REPO / "outputs" / "stage1" / "soil_profiles_0_100cm.gpkg")
    rows_corr, rows_zone = [], []

    for depth in ["0_30cm", "0_100cm"]:
        q05, prof = rd(S5 / f"ales_index_q05_{depth}.tif")
        q50, _ = rd(S5 / f"ales_index_q50_{depth}.tif")
        q95, _ = rd(S5 / f"ales_index_q95_{depth}.tif")
        p2, _ = rd(S5 / f"p_at_least_S2_{depth}.tif")
        valid = np.isfinite(q50)

        piw = q95 - q05
        ncls = np.abs(classify(np.where(valid, q95, 0)).astype(np.int16)
                      - classify(np.where(valid, q05, 0)).astype(np.int16)) + 1
        pc = np.clip(p2, 1e-6, 1 - 1e-6)
        ambig = (-(pc * np.log2(pc) + (1 - pc) * np.log2(1 - pc))).astype("float32")
        dbreak = np.min(np.abs(q50[..., None] - BREAKS[None, None, :]), axis=-1).astype("float32")

        ks, _ = rd(STAGE3B / f"Ks_{depth}_q50_qrf.tif")
        clay, _ = rd(STAGE3B / f"Clay_{depth}_q50_qrf.tif")
        kssd = local_sd(ks, valid)
        claysd = local_sd(clay, valid)

        # distance to nearest profile (in metres, on the 30 m grid)
        tr = prof["transform"]
        seed = np.ones(q50.shape, bool)
        for geom in gdf.geometry:
            c, r = ~tr * (geom.x, geom.y)
            r, c = int(round(r)), int(round(c))
            if 0 <= r < q50.shape[0] and 0 <= c < q50.shape[1]:
                seed[r, c] = False
        dprof = (ndimage.distance_transform_edt(seed) * 30.0).astype("float32")
        dedge = (ndimage.distance_transform_edt(valid) * 30.0).astype("float32")

        m = dict(prof, dtype="float32", count=1, nodata=np.float32(np.nan), compress="deflate")
        for nm, arr in [("piw", piw), ("ambiguity", ambig), ("n_classes_spanned", ncls.astype("float32")),
                        ("dist_to_profile_m", dprof), ("ks_local_sd", kssd)]:
            with rasterio.open(OUT / f"unc_{nm}_{depth}.tif", "w", **m) as d:
                d.write(np.where(valid, arr, np.nan).astype("float32"), 1)

        # driver correlations
        V = valid & np.isfinite(piw) & np.isfinite(kssd) & np.isfinite(claysd)
        drivers = {"dist_to_profile_m": dprof, "ks_local_sd_5x5": kssd,
                   "clay_local_sd_5x5": claysd, "ks_q50": ks,
                   "dist_to_mask_edge_m": dedge, "index_q50": q50,
                   "dist_to_class_break": dbreak}
        for target_name, target in [("PI_width", piw), ("ambiguity_entropy", ambig)]:
            for dn, dv in drivers.items():
                r, p = spearmanr(target[V].ravel(), dv[V].ravel())
                rows_corr.append(dict(depth=depth, target=target_name, driver=dn,
                                      spearman_rho=round(float(r), 3),
                                      p_value=f"{p:.2e}", n=int(V.sum())))
        # uncertainty by suitability class
        cls = classify(np.where(valid, q50, 0))
        lab = {1: "S1", 2: "S2", 3: "S3", 4: "N1", 5: "N2"}
        for c, nm in lab.items():
            sel = V & (cls == c)
            if sel.sum() < 50:
                continue
            rows_zone.append(dict(depth=depth, zone=f"class {nm}", n_pixels=int(sel.sum()),
                                  mean_PI_width=round(float(piw[sel].mean()), 2),
                                  mean_ambiguity=round(float(ambig[sel].mean()), 3),
                                  mean_dist_to_profile_m=round(float(dprof[sel].mean()), 0),
                                  mean_ks_local_sd=round(float(kssd[sel].mean()), 4)))
        # transition vs interior: top decile of local Ks sd = transition zone
        thr = np.nanpercentile(kssd[V], 90)
        for nm, sel in [("transition (top decile local Ks SD)", V & (kssd >= thr)),
                        ("interior (lower nine deciles)", V & (kssd < thr))]:
            rows_zone.append(dict(depth=depth, zone=nm, n_pixels=int(sel.sum()),
                                  mean_PI_width=round(float(piw[sel].mean()), 2),
                                  mean_ambiguity=round(float(ambig[sel].mean()), 3),
                                  mean_dist_to_profile_m=round(float(dprof[sel].mean()), 0),
                                  mean_ks_local_sd=round(float(kssd[sel].mean()), 4)))
        # share of pixels whose interval spans >1 class
        multi = float(((ncls > 1) & V).sum()) / float(V.sum()) * 100
        print(f"{depth}: pixels whose 90% interval spans >1 suitability class = {multi:.1f}%", flush=True)
        rows_zone.append(dict(depth=depth, zone="ALL (pct spanning >1 class)",
                              n_pixels=int(V.sum()), mean_PI_width=round(multi, 1),
                              mean_ambiguity=np.nan, mean_dist_to_profile_m=np.nan,
                              mean_ks_local_sd=np.nan))

    pd.DataFrame(rows_corr).to_csv(OUT / "unc_driver_correlations.csv", index=False)
    pd.DataFrame(rows_zone).to_csv(OUT / "unc_by_zone.csv", index=False)
    print("\n===== uncertainty drivers =====")
    print(pd.DataFrame(rows_corr).to_string(index=False))
    print("\n===== uncertainty by zone =====")
    print(pd.DataFrame(rows_zone).to_string(index=False))

if __name__ == "__main__":
    main()
