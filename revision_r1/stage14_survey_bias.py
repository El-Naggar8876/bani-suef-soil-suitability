"""
Stage 14 - Quantifying the corridor-sampling bias in the 993-point crop survey.

Purpose (Reviewer 1, comment 6): the survey followed roads, canal banks and
field tracks rather than a probability design, so per-crop statistics may be
biased. This script measures the bias directly by comparing the surveyed points
with the full population of cropped pixels, and reports a reweighted estimate.

Outputs -> 02_New_Analysis/outputs_rev/
"""
from __future__ import annotations
import pathlib
import numpy as np, pandas as pd, rasterio
from scipy import ndimage
from scipy.stats import ks_2samp

import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parent))
from revpaths import REPO as _REPO, OUT as _OUT, WAPOR as _WAPOR

HERE = pathlib.Path(__file__).resolve().parent
REPO = _REPO
OUT = _OUT
S5, S6, S7, S3B = (REPO / "outputs" / s for s in ("stage5", "stage6", "stage7", "stage3b"))
DEPTH = "0_100cm"


def rd(p):
    with rasterio.open(p) as s:
        a = s.read(1).astype("float32")
        if s.nodata is not None and np.isfinite(s.nodata):
            a[a == s.nodata] = np.nan
        return a, s.transform


def main():
    idx, tr = rd(S5 / f"ales_index_q50_{DEPTH}.tif")
    cls, _ = rd(S5 / f"ales_class_q50_{DEPTH}.tif")
    cfd, _ = rd(S6 / f"cf_CF_Ks_{DEPTH}_delta_q50.tif")
    ks, _ = rd(S3B / f"Ks_{DEPTH}_q50_qrf.tif")
    pri, _ = rd(S7 / f"policy_priority_{DEPTH}.tif")
    lu, _ = rd(S7 / "landuse_class.tif")
    cropped = np.isin(lu, [2, 3]) & np.isfinite(idx)

    pts = pd.read_csv(REPO / "outputs" / "stage7b" / "crop_points_extracted.csv")
    rc = [(~tr * (x, y)) for x, y in zip(pts.x, pts.y)]
    pts["row"] = [int(r) for c, r in rc]
    pts["col"] = [int(c) for c, r in rc]
    inb = (pts.row.between(0, idx.shape[0] - 1)) & (pts.col.between(0, idx.shape[1] - 1))
    pts = pts[inb].copy()

    # distance from each cropped pixel to the nearest survey point -> corridor proxy
    seed = np.ones(idx.shape, bool)
    seed[pts.row.values, pts.col.values] = False
    dsurv = (ndimage.distance_transform_edt(seed) * 30.0).astype("float32")

    rows = []
    for vn, arr in [("ALES_index_q50", idx), ("CF_Ks_delta_q50", cfd), ("Ks_q50", ks)]:
        sv = arr[pts.row.values, pts.col.values]
        sv = sv[np.isfinite(sv)]
        pop = arr[cropped]
        pop = pop[np.isfinite(pop)]
        st = ks_2samp(sv, pop)
        rows.append(dict(variable=vn, n_survey=sv.size, n_cropped_pixels=pop.size,
                         survey_mean=round(float(sv.mean()), 3),
                         cropped_population_mean=round(float(pop.mean()), 3),
                         difference=round(float(sv.mean() - pop.mean()), 3),
                         pct_difference=round(100 * (sv.mean() - pop.mean()) / pop.mean(), 2),
                         ks_statistic=round(float(st.statistic), 4),
                         ks_p_value=f"{st.pvalue:.2e}"))
    b1 = pd.DataFrame(rows)
    b1.to_csv(OUT / "survey_bias_distributions.csv", index=False)
    print("===== survey points vs all cropped pixels =====")
    print(b1.to_string(index=False))

    # class and priority composition
    comp = []
    for nm, lab, codes in [("ALES class", cls, {1: "S1", 2: "S2", 3: "S3", 4: "N1", 5: "N2"}),
                           ("priority", pri, {1: "HIGH", 2: "MEDIUM", 3: "WATCH", 0: "none"})]:
        sv = lab[pts.row.values, pts.col.values]
        pop = lab[cropped]
        for code, name in codes.items():
            comp.append(dict(stratum=f"{nm} {name}",
                             survey_pct=round(100 * float(np.mean(sv == code)), 2),
                             cropped_population_pct=round(100 * float(np.mean(pop == code)), 2),
                             over_representation_pct_pts=round(
                                 100 * float(np.mean(sv == code)) - 100 * float(np.mean(pop == code)), 2)))
    b2 = pd.DataFrame(comp)
    b2.to_csv(OUT / "survey_bias_composition.csv", index=False)
    print("\n===== composition: surveyed vs cropped population =====")
    print(b2.to_string(index=False))

    # post-stratified reweighting: weight each survey point by the ratio of the
    # population share to the sample share of its ALES class
    sv_cls = cls[pts.row.values, pts.col.values]
    pop_cls = cls[cropped]
    w = np.ones(len(pts))
    for code in np.unique(sv_cls[np.isfinite(sv_cls)]):
        s_share = float(np.mean(sv_cls == code))
        p_share = float(np.mean(pop_cls == code))
        if s_share > 0:
            w[sv_cls == code] = p_share / s_share
    pts["w"] = w
    rows = []
    for crop, g in pts.groupby("CROP"):
        if len(g) < 5:
            continue
        n1 = (g.ales_cls >= 4).astype(float)
        hi = (g.priority == 1).astype(float)
        rows.append(dict(crop=crop, n=len(g),
                         pct_N1_unweighted=round(100 * float(n1.mean()), 1),
                         pct_N1_reweighted=round(100 * float(np.average(n1, weights=g.w)), 1),
                         pct_HIGH_unweighted=round(100 * float(hi.mean()), 1),
                         pct_HIGH_reweighted=round(100 * float(np.average(hi, weights=g.w)), 1),
                         mean_CFKs_unweighted=round(float(g.cf_Ks_d.mean()), 2),
                         mean_CFKs_reweighted=round(float(np.average(g.cf_Ks_d, weights=g.w)), 2)))
    b3 = pd.DataFrame(rows).sort_values("n", ascending=False)
    b3.to_csv(OUT / "survey_bias_reweighted_percrop.csv", index=False)
    print("\n===== per-crop statistics, unweighted vs post-stratified =====")
    print(b3.to_string(index=False))

    # how much of the cropped area lies close to a surveyed corridor
    rows = []
    for d in (250, 500, 1000, 2000, 5000):
        rows.append(dict(buffer_m=d,
                         pct_cropped_area_within=round(100 * float(np.mean(dsurv[cropped] <= d)), 1),
                         mean_index_within=round(float(idx[cropped & (dsurv <= d)].mean()), 2),
                         mean_index_beyond=round(float(idx[cropped & (dsurv > d)].mean()), 2) if (cropped & (dsurv > d)).sum() else np.nan))
    b4 = pd.DataFrame(rows)
    b4.to_csv(OUT / "survey_corridor_coverage.csv", index=False)
    print("\n===== corridor coverage of the cropped area =====")
    print(b4.to_string(index=False))

if __name__ == "__main__":
    main()
