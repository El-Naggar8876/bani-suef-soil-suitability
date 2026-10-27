"""
Stage 13 - Independent remote-sensing check of the modelled drainage constraint.

Purpose: the suitability surface and the counterfactual decomposition are both
model products built from the same 60 profiles. FAO WaPOR v3 (100 m, L2) supplies
water-consumption and biomass fields that enter neither the soil model nor the
ALES rating, and therefore provide an external observational test of the central
claim. If low saturated hydraulic conductivity restricts root-zone aeration and
drainage, cropped land in the Ks-limited priority zone should show lower
transpiration and lower biomass production than comparable cropped land outside
it, after crop type is held constant.

Because neighbouring 30 m pixels are strongly autocorrelated, nominal p-values
from pixel counts would be meaningless. Confidence intervals are therefore
obtained by a moving-block bootstrap over 1 km x 1 km spatial blocks.

Outputs -> 02_New_Analysis/outputs_rev/
"""
from __future__ import annotations
import sys, pathlib
import numpy as np, pandas as pd, rasterio
from scipy.stats import spearmanr

import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parent))
from revpaths import REPO as _REPO, OUT as _OUT, WAPOR as _WAPOR

HERE = pathlib.Path(__file__).resolve().parent
REPO = _REPO
OUT = _OUT
W = _WAPOR
S5, S6, S7, S3B = (REPO / "outputs" / s for s in ("stage5", "stage6", "stage7", "stage3b"))
DEPTH = "0_100cm"
BLOCK = 33           # 33 x 30 m ~ 1 km
NBOOT = 400
RNG = np.random.default_rng(11)


def rd(p):
    with rasterio.open(p) as s:
        a = s.read(1).astype("float32")
        if s.nodata is not None and np.isfinite(s.nodata):
            a[a == s.nodata] = np.nan
        return a


def block_ci(values, blocks, stat=np.mean, nboot=NBOOT, q=(2.5, 97.5)):
    """Moving-block bootstrap CI for a statistic of `values`, resampling block ids."""
    ub = np.unique(blocks)
    if ub.size < 5:
        return (np.nan, np.nan)
    idx_by_block = {b: np.where(blocks == b)[0] for b in ub}
    out = np.empty(nboot)
    for i in range(nboot):
        pick = RNG.choice(ub, ub.size, replace=True)
        sel = np.concatenate([idx_by_block[b] for b in pick])
        out[i] = stat(values[sel])
    return tuple(np.percentile(out, q))


def main():
    lay = {
        "T_2025_mm": rd(W / "wapor_T_2025_30m.tif"),
        "T_mean_mm": rd(W / "wapor_T_mean_2021_2025_30m.tif"),
        "AETI_mean_mm": rd(W / "wapor_AETI_mean_2021_2025_30m.tif"),
        "TBP_mean_kg_ha": rd(W / "wapor_TBP_mean_2021_2025_30m.tif"),
    }
    # biomass water productivity, kg dry matter per m3 of water consumed
    lay["WP_kg_m3"] = lay["TBP_mean_kg_ha"] / (lay["AETI_mean_mm"] * 10.0)

    pri = rd(S7 / f"policy_priority_{DEPTH}.tif")
    lu = rd(S7 / "landuse_class.tif")
    cls = rd(S5 / f"ales_class_q50_{DEPTH}.tif")
    cfd = rd(S6 / f"cf_CF_Ks_{DEPTH}_delta_q50.tif")
    ks = rd(S3B / f"Ks_{DEPTH}_q50_qrf.tif")
    idx = rd(S5 / f"ales_index_q50_{DEPTH}.tif")

    H, Wd = pri.shape
    rr, cc = np.mgrid[0:H, 0:Wd]
    blk = (rr // BLOCK) * 10000 + (cc // BLOCK)

    cropped = np.isin(lu, [2, 3])
    base = np.isfinite(idx) & cropped
    for v in lay.values():
        base &= np.isfinite(v)
    print(f"cropped pixels with complete WaPOR coverage: {int(base.sum())}")

    rows = []
    zones = {
        "HIGH priority (cropped, N1, top-quartile CF_Ks)": base & (pri == 1),
        "MEDIUM priority (cropped, N1, lower quartiles)": base & (pri == 2),
        "WATCH (cropped, not N1)": base & (pri == 3),
        "all cropped S3": base & (cls == 3),
        "all cropped N1": base & (cls == 4),
    }
    for zn, m in zones.items():
        if m.sum() < 500:
            continue
        r = dict(zone=zn, n_pixels=int(m.sum()), area_km2=round(m.sum() * 0.0009, 1))
        for vn, arr in lay.items():
            vals = arr[m]
            lo, hi = block_ci(vals, blk[m])
            r[f"{vn}_mean"] = round(float(vals.mean()), 2)
            r[f"{vn}_ci95"] = f"[{lo:.1f}, {hi:.1f}]"
        rows.append(r)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "wapor_by_priority_zone.csv", index=False)
    print("\n===== WaPOR by priority zone (cropped land only) =====")
    print(df[["zone", "n_pixels", "T_mean_mm_mean", "T_mean_mm_ci95",
              "AETI_mean_mm_mean", "TBP_mean_kg_ha_mean", "TBP_mean_kg_ha_ci95",
              "WP_kg_m3_mean", "WP_kg_m3_ci95"]].to_string(index=False))

    # ---- paired HIGH vs rest-of-cropped contrast, with block bootstrap on the difference
    con = []
    hi_m = base & (pri == 1)
    for other_name, other_m in [("MEDIUM+WATCH cropped", base & np.isin(pri, [2, 3])),
                                ("cropped S3", base & (cls == 3))]:
        for vn, arr in lay.items():
            a, b = arr[hi_m], arr[other_m]
            d = float(a.mean() - b.mean())
            ub_a, ub_b = np.unique(blk[hi_m]), np.unique(blk[other_m])
            ia = {x: np.where(blk[hi_m] == x)[0] for x in ub_a}
            ib = {x: np.where(blk[other_m] == x)[0] for x in ub_b}
            boot = np.empty(NBOOT)
            for i in range(NBOOT):
                pa = RNG.choice(ub_a, ub_a.size, replace=True)
                pb = RNG.choice(ub_b, ub_b.size, replace=True)
                boot[i] = (a[np.concatenate([ia[x] for x in pa])].mean()
                           - b[np.concatenate([ib[x] for x in pb])].mean())
            lo, hi = np.percentile(boot, [2.5, 97.5])
            con.append(dict(variable=vn, contrast=f"HIGH minus {other_name}",
                            difference=round(d, 2), ci95_lo=round(float(lo), 2),
                            ci95_hi=round(float(hi), 2),
                            pct_of_comparison=round(100 * d / float(b.mean()), 1),
                            excludes_zero=bool(lo > 0 or hi < 0)))
    cdf = pd.DataFrame(con)
    cdf.to_csv(OUT / "wapor_contrasts.csv", index=False)
    print("\n===== HIGH-priority contrasts (block-bootstrap 95% CI) =====")
    print(cdf.to_string(index=False))

    # ---- WaPOR response along the predicted Ks gradient, cropped land only
    g = []
    kv = ks[base]
    qs = np.percentile(kv, [0, 20, 40, 60, 80, 100])
    for i in range(5):
        m = base & (ks >= qs[i]) & (ks <= qs[i + 1] if i == 4 else ks < qs[i + 1])
        if m.sum() < 500:
            continue
        row = dict(ks_quintile=f"Q{i+1}", ks_range=f"{qs[i]:.2f}-{qs[i+1]:.2f}",
                   n_pixels=int(m.sum()), mean_ks=round(float(ks[m].mean()), 3))
        for vn, arr in lay.items():
            row[vn] = round(float(arr[m].mean()), 2)
        g.append(row)
    gdf = pd.DataFrame(g)
    gdf.to_csv(OUT / "wapor_by_ks_quintile.csv", index=False)
    print("\n===== WaPOR along the predicted Ks gradient (cropped land) =====")
    print(gdf.to_string(index=False))

    # ---- rank correlations on a thinned sample (every 5th pixel both ways)
    thin = np.zeros_like(base)
    thin[::5, ::5] = True
    m = base & thin
    cr = []
    for vn, arr in lay.items():
        for dn, dv in [("Ks_q50", ks), ("ALES_index_q50", idx), ("CF_Ks_delta_q50", cfd)]:
            r, p = spearmanr(arr[m], dv[m])
            cr.append(dict(wapor_variable=vn, model_variable=dn,
                           spearman_rho=round(float(r), 3), p_value=f"{p:.2e}",
                           n_thinned=int(m.sum())))
    crdf = pd.DataFrame(cr)
    crdf.to_csv(OUT / "wapor_model_correlations.csv", index=False)
    print("\n===== WaPOR vs model variables (thinned, 150 m spacing) =====")
    print(crdf.to_string(index=False))

if __name__ == "__main__":
    main()
