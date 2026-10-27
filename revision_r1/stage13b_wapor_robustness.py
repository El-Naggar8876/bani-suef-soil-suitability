"""
Stage 13b - Robustness checks on the WaPOR contrast.

The Stage 13 result (transpiration and biomass DECREASE as predicted Ks
increases) is the opposite of what a drainage-limitation reading predicts. Three
competing explanations are tested before the result is interpreted:

  C1  Mixed-pixel artefact. WaPOR is 100 m; cropped pixels near the desert
      margin contain bare soil, depressing T and TBP. Tested by eroding the
      cropped mask and by restricting to double-cropped pixels only.
  C2  Soil-type confounding. In this floodplain low Ks may coincide with heavy,
      fertile central alluvium and high Ks with coarse, calcareous desert-fringe
      soils. Tested by correlating Ks with texture and chemistry and with
      distance from the cultivated core, and by partial correlation.
  C3  Crop-composition confounding. The priority zones differ in crop mix.
      Tested at the 993 survey points, within single crops.

Outputs -> 02_New_Analysis/outputs_rev/
"""
from __future__ import annotations
import pathlib
import numpy as np, pandas as pd, rasterio
from scipy import ndimage
from scipy.stats import spearmanr

import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parent))
from revpaths import REPO as _REPO, OUT as _OUT, WAPOR as _WAPOR

HERE = pathlib.Path(__file__).resolve().parent
REPO = _REPO
OUT = _OUT
Wd = _WAPOR
S5, S6, S7, S3B = (REPO / "outputs" / s for s in ("stage5", "stage6", "stage7", "stage3b"))
DEPTH = "0_100cm"


def rd(p):
    with rasterio.open(p) as s:
        a = s.read(1).astype("float32")
        if s.nodata is not None and np.isfinite(s.nodata):
            a[a == s.nodata] = np.nan
        return a


def partial_spearman(x, y, Z):
    """Spearman partial correlation of x,y given covariates Z (rank-residual)."""
    from scipy.stats import rankdata
    R = np.column_stack([rankdata(z) for z in Z] + [np.ones(len(x))])
    def resid(v):
        v = rankdata(v).astype(float)
        beta, *_ = np.linalg.lstsq(R, v, rcond=None)
        return v - R @ beta
    rx, ry = resid(x), resid(y)
    return float(np.corrcoef(rx, ry)[0, 1])


def main():
    T = rd(Wd / "wapor_T_mean_2021_2025_30m.tif")
    TBP = rd(Wd / "wapor_TBP_mean_2021_2025_30m.tif")
    lu = rd(S7 / "landuse_class.tif")
    pri = rd(S7 / f"policy_priority_{DEPTH}.tif")
    ks = rd(S3B / f"Ks_{DEPTH}_q50_qrf.tif")
    clay = rd(S3B / f"Clay_{DEPTH}_q50_qrf.tif")
    sand = rd(S3B / f"Sand_total_{DEPTH}_q50_qrf.tif")
    ec = rd(S3B / f"EC_{DEPTH}_q50_qrf.tif")
    caco3 = rd(S3B / f"CaCO3_{DEPTH}_q50_qrf.tif")
    som = rd(S3B / f"SOM_{DEPTH}_q50_qrf.tif")
    idx = rd(S5 / f"ales_index_q50_{DEPTH}.tif")

    cropped = np.isin(lu, [2, 3])
    dbl = (lu == 2)
    fin = np.isfinite(T) & np.isfinite(TBP) & np.isfinite(ks) & np.isfinite(idx)

    # ---------------- C1 mixed-pixel artefact -------------------------------
    rows = []
    for name, m0 in [("all cropped (single+double)", cropped),
                     ("double-cropped only", dbl)]:
        for erode_px in (0, 2, 4, 7):
            m = ndimage.binary_erosion(m0, np.ones((2 * erode_px + 1,) * 2)) if erode_px else m0
            m = m & fin
            if m.sum() < 2000:
                continue
            kv, tv, bv = ks[m], T[m], TBP[m]
            rows.append(dict(subset=name, erosion_px=erode_px, erosion_m=erode_px * 30,
                             n_pixels=int(m.sum()),
                             rho_T_vs_Ks=round(spearmanr(tv, kv).statistic, 3),
                             rho_TBP_vs_Ks=round(spearmanr(bv, kv).statistic, 3),
                             T_lowKs_q1=round(float(tv[kv <= np.percentile(kv, 20)].mean()), 1),
                             T_highKs_q5=round(float(tv[kv >= np.percentile(kv, 80)].mean()), 1),
                             TBP_lowKs_q1=round(float(bv[kv <= np.percentile(kv, 20)].mean()), 0),
                             TBP_highKs_q5=round(float(bv[kv >= np.percentile(kv, 80)].mean()), 0)))
    c1 = pd.DataFrame(rows)
    c1.to_csv(OUT / "wapor_C1_mixed_pixel.csv", index=False)
    print("===== C1  mixed-pixel / edge-mixing test =====")
    print(c1.to_string(index=False))

    # ---------------- C2 soil-type confounding ------------------------------
    thin = np.zeros_like(cropped); thin[::4, ::4] = True
    m = cropped & fin & thin & np.isfinite(clay) & np.isfinite(ec) & np.isfinite(caco3) & np.isfinite(som) & np.isfinite(sand)
    # distance from the cultivated core: distance to nearest non-cropped pixel
    dcore = (ndimage.distance_transform_edt(cropped) * 30.0).astype("float32")
    rows = []
    for vn, v in [("Clay_pct", clay), ("Sand_pct", sand), ("CaCO3_pct", caco3),
                  ("EC_dS_m", ec), ("SOM_pct", som),
                  ("dist_to_cropland_edge_m", dcore), ("ALES_index_q50", idx)]:
        rows.append(dict(variable=vn,
                         rho_with_Ks=round(spearmanr(ks[m], v[m]).statistic, 3),
                         rho_with_T=round(spearmanr(T[m], v[m]).statistic, 3),
                         rho_with_TBP=round(spearmanr(TBP[m], v[m]).statistic, 3)))
    c2 = pd.DataFrame(rows)
    c2.to_csv(OUT / "wapor_C2_confounding.csv", index=False)
    print("\n===== C2  what co-varies with predicted Ks (cropped, thinned 120 m) =====")
    print(c2.to_string(index=False))

    Z = [clay[m], sand[m], caco3[m], ec[m], som[m], dcore[m]]
    pr_T = partial_spearman(ks[m], T[m], Z)
    pr_B = partial_spearman(ks[m], TBP[m], Z)
    print(f"\npartial rho(T, Ks | clay, sand, CaCO3, EC, SOM, edge-distance) = {pr_T:.3f}")
    print(f"partial rho(TBP, Ks | same covariates)                        = {pr_B:.3f}")
    pd.DataFrame([dict(pair="T vs Ks", raw_rho=round(spearmanr(T[m], ks[m]).statistic, 3),
                       partial_rho=round(pr_T, 3)),
                  dict(pair="TBP vs Ks", raw_rho=round(spearmanr(TBP[m], ks[m]).statistic, 3),
                       partial_rho=round(pr_B, 3))]).to_csv(
        OUT / "wapor_C2_partial.csv", index=False)

    # ---------------- C3 crop-composition confounding -----------------------
    pts = pd.read_csv(REPO / "outputs" / "stage7b" / "crop_points_extracted.csv")
    with rasterio.open(Wd / "wapor_T_mean_2021_2025_30m.tif") as s:
        tr = s.transform
        pts["T_mean"] = [T[int((~tr * (x, y))[1]), int((~tr * (x, y))[0])]
                         if 0 <= int((~tr * (x, y))[1]) < T.shape[0] and 0 <= int((~tr * (x, y))[0]) < T.shape[1]
                         else np.nan for x, y in zip(pts.x, pts.y)]
        pts["TBP_mean"] = [TBP[int((~tr * (x, y))[1]), int((~tr * (x, y))[0])]
                           if 0 <= int((~tr * (x, y))[1]) < TBP.shape[0] and 0 <= int((~tr * (x, y))[0]) < TBP.shape[1]
                           else np.nan for x, y in zip(pts.x, pts.y)]
        pts["Ks"] = [ks[int((~tr * (x, y))[1]), int((~tr * (x, y))[0])]
                     if 0 <= int((~tr * (x, y))[1]) < ks.shape[0] and 0 <= int((~tr * (x, y))[0]) < ks.shape[1]
                     else np.nan for x, y in zip(pts.x, pts.y)]
    pts = pts.dropna(subset=["T_mean", "TBP_mean", "Ks"])
    rows = []
    for crop, g in pts.groupby("CROP"):
        if len(g) < 15:
            continue
        hi = g[g.priority == 1]; rest = g[g.priority != 1]
        rows.append(dict(crop=crop, n=len(g), n_high=len(hi),
                         rho_T_vs_Ks=round(spearmanr(g.T_mean, g.Ks).statistic, 3),
                         rho_TBP_vs_Ks=round(spearmanr(g.TBP_mean, g.Ks).statistic, 3),
                         T_high=round(float(hi.T_mean.mean()), 1) if len(hi) >= 5 else np.nan,
                         T_rest=round(float(rest.T_mean.mean()), 1) if len(rest) >= 5 else np.nan,
                         TBP_high=round(float(hi.TBP_mean.mean()), 0) if len(hi) >= 5 else np.nan,
                         TBP_rest=round(float(rest.TBP_mean.mean()), 0) if len(rest) >= 5 else np.nan))
    c3 = pd.DataFrame(rows)
    c3.to_csv(OUT / "wapor_C3_crop_controlled.csv", index=False)
    print("\n===== C3  within-crop test at the 993 survey points =====")
    print(c3.to_string(index=False))

if __name__ == "__main__":
    main()
