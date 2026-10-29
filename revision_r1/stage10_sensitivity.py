"""
Stage 10 - ALES structural sensitivity analysis.

Purpose (Reviewer 2, comment 4.1): establish whether the dominance of saturated
hydraulic conductivity (Ks) in the counterfactual decomposition is an artefact of
the multiplicative aggregation rule, or a genuine property of the input data.

Four independent tests:
  T1  Limiting-factor census     - which sub-rating is minimum at each pixel.
  T2  Normalised elasticity      - d(ln index)/d(ln property) by finite difference.
  T3  Aggregation-rule swap      - repeat the counterfactual ranking under four
                                   aggregation rules (multiplicative, geometric
                                   mean, arithmetic mean, minimum).
  T4  Threshold perturbation     - shift the r_wetness lower optimum by +/-50 %
                                   and re-rank the counterfactuals.

Outputs -> 02_New_Analysis/outputs_rev/
"""
from __future__ import annotations
import sys, json, pathlib
import numpy as np, pandas as pd, rasterio

import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parent))
from revpaths import REPO as _REPO, OUT as _OUT, WAPOR as _WAPOR

HERE = pathlib.Path(__file__).resolve().parent
REPO = _REPO
sys.path.insert(0, str(REPO / "notebooks"))
OUT = _OUT

from stage5_ales_montecarlo import (          # identical implementation to submission
    r_topo, r_wetness, r_physical, r_chemical, r_fertility,
    r_ec, r_esp, r_caco3, r_cec, r_som, piecewise, classify, ales_index,
    STAGE3B, STACK, STAGE2B, PROPS, DEPTHS,
)

PROP_KEY = {"EC": "EC", "ESP": "ESP", "CaCO3": "CaCO3", "Clay": "Clay",
            "Sand": "Sand_total", "CEC": "CEC", "SOM": "SOM", "Ks": "Ks"}

def load(prop, depth, key="q50_qrf"):
    with rasterio.open(STAGE3B / f"{prop}_{depth}_{key}.tif") as s:
        a = s.read(1).astype("float32")
        return a

def get_slope():
    with rasterio.open(STACK) as s:
        names = list(s.descriptions) or (STAGE2B / "band_names.txt").read_text().splitlines()
        return s.read(names.index("slope") + 1).astype("float32")

# ---------------------------------------------------------------- aggregation rules
def aggregate(rt, rw, rp, rc, rf, rule):
    R = np.stack([rt, rw, rp, rc, rf])
    if rule == "multiplicative":                     # as submitted (Sys-Verheye)
        return (rt * rw * rp * rc * rf) / 100.0 ** 4
    if rule == "geometric":
        return np.exp(np.mean(np.log(np.clip(R, 1e-6, None)), axis=0))
    if rule == "arithmetic":
        return R.mean(axis=0)
    if rule == "minimum":                            # Liebig / Storie-type
        return R.min(axis=0)
    raise ValueError(rule)

RULES = ["multiplicative", "geometric", "arithmetic", "minimum"]

# counterfactual targets, incl. threshold-consistent CaCO3 variant
SCEN = {
    "CF_EC":        {"EC": (2.0, "cap")},
    "CF_ESP":       {"ESP": (10.0, "cap")},
    "CF_CaCO3":     {"CaCO3": (15.0, "cap")},
    "CF_CaCO3_10":  {"CaCO3": (10.0, "cap")},        # NEW: at implemented no-penalty point
    "CF_Ks":        {"Ks": (2.0, "lift")},
    "CF_SOM":       {"SOM": (1.0, "lift")},
    "CF_chem":      {"EC": (2.0, "cap"), "ESP": (10.0, "cap"), "CaCO3": (15.0, "cap")},
    "CF_chem_10":   {"EC": (2.0, "cap"), "ESP": (10.0, "cap"), "CaCO3": (10.0, "cap")},
    "CF_full":      {"EC": (2.0, "cap"), "ESP": (10.0, "cap"), "CaCO3": (15.0, "cap"),
                     "Ks": (2.0, "lift"), "SOM": (1.0, "lift")},
}

def apply_scen(vals, spec):
    out = dict(vals)
    for p, (tgt, mode) in spec.items():
        out[p] = np.minimum(vals[p], tgt) if mode == "cap" else np.maximum(vals[p], tgt)
    return out

def subratings(v, slope, ks_lo=1.5):
    rw = r_wetness_shift(v["Ks"], ks_lo)
    return (r_topo(slope), rw, r_physical(v["Clay"], v["Sand_total"]),
            r_chemical(v["EC"], v["ESP"], v["CaCO3"]), r_fertility(v["CEC"], v["SOM"]))

def r_wetness_shift(ks, lo):
    """r_wetness with the lower optimum threshold moved from 1.5 to `lo`."""
    return np.where(ks < lo,
                    piecewise(ks, [0, 0.1 * lo / 1.5, 0.5 * lo / 1.5, lo], [30, 30, 60, 85]),
                    np.where(ks <= 15, 100,
                             np.where(ks <= 30, 85, np.where(ks <= 60, 60, 30)))).astype("float32")

# ---------------------------------------------------------------------------- main
def main():
    slope = get_slope()
    rows_T1, rows_T2, rows_T3, rows_T4 = [], [], [], []

    for depth in DEPTHS:
        v = {p: load(p, depth) for p in PROPS}
        valid = np.isfinite(v["Ks"]) & np.isfinite(slope)
        n = int(valid.sum())

        # ---- T1 limiting-factor census -------------------------------------
        rt, rw, rp, rc, rf = subratings(v, slope)
        S = np.stack([rt, rw, rp, rc, rf])
        names = ["topographic", "wetness (Ks)", "physical (texture)",
                 "chemical (EC/ESP/CaCO3)", "fertility (CEC/SOM)"]
        amin = np.argmin(np.where(np.isfinite(S), S, np.inf), axis=0)
        for i, nm in enumerate(names):
            share = float(((amin == i) & valid).sum()) / n * 100
            sub = S[i][valid]
            rows_T1.append(dict(depth=depth, sub_rating=nm, pct_pixels_limiting=round(share, 2),
                                mean_sub_rating=round(float(np.nanmean(sub)), 2),
                                min_sub_rating=round(float(np.nanmin(sub)), 2),
                                max_sub_rating=round(float(np.nanmax(sub)), 2)))

        base_idx = aggregate(rt, rw, rp, rc, rf, "multiplicative")

        # ---- T2 normalised elasticity --------------------------------------
        # e = d(ln I)/d(ln x) by central difference at +/-5 % of each property
        for label, key in PROP_KEY.items():
            up, dn = dict(v), dict(v)
            up[key] = v[key] * 1.05
            dn[key] = v[key] * 0.95
            iu = aggregate(*subratings(up, slope), "multiplicative")
            idn = aggregate(*subratings(dn, slope), "multiplicative")
            with np.errstate(divide="ignore", invalid="ignore"):
                e = (np.log(np.clip(iu, 1e-9, None)) - np.log(np.clip(idn, 1e-9, None))) / (np.log(1.05) - np.log(0.95))
            e = e[valid & np.isfinite(e)]
            rows_T2.append(dict(depth=depth, property=label,
                                mean_elasticity=round(float(e.mean()), 4),
                                median_elasticity=round(float(np.median(e)), 4),
                                p95_elasticity=round(float(np.percentile(e, 95)), 4),
                                pct_pixels_responsive=round(float((np.abs(e) > 0.01).mean() * 100), 2)))

        # ---- T3 aggregation-rule swap --------------------------------------
        for rule in RULES:
            b = aggregate(rt, rw, rp, rc, rf, rule)
            for sc, spec in SCEN.items():
                vv = apply_scen(v, spec)
                idx = aggregate(*subratings(vv, slope), rule)
                d = (idx - b)[valid]
                rows_T3.append(dict(depth=depth, rule=rule, scenario=sc,
                                    mean_delta=round(float(np.nanmean(d)), 3),
                                    baseline_mean=round(float(np.nanmean(b[valid])), 3)))

        # ---- T4 r_wetness threshold perturbation ---------------------------
        for lo in [0.75, 1.125, 1.5, 1.875, 2.25]:
            rt2, rw2, rp2, rc2, rf2 = subratings(v, slope, ks_lo=lo)
            b = aggregate(rt2, rw2, rp2, rc2, rf2, "multiplicative")
            for sc in ["CF_Ks", "CF_EC", "CF_SOM", "CF_full"]:
                vv = apply_scen(v, SCEN[sc])
                idx = aggregate(*subratings(vv, slope, ks_lo=lo), "multiplicative")
                d = (idx - b)[valid]
                rows_T4.append(dict(depth=depth, ks_lower_optimum=lo, scenario=sc,
                                    mean_delta=round(float(np.nanmean(d)), 3),
                                    baseline_mean=round(float(np.nanmean(b[valid])), 3)))
        print(f"  depth {depth} done ({n} valid pixels)", flush=True)

    for nm, rows in [("T1_limiting_factor", rows_T1), ("T2_elasticity", rows_T2),
                     ("T3_aggregation_rule", rows_T3), ("T4_threshold_perturbation", rows_T4)]:
        df = pd.DataFrame(rows); df.to_csv(OUT / f"sens_{nm}.csv", index=False)
        print(f"\n===== {nm} =====\n{df.to_string(index=False)}")

if __name__ == "__main__":
    main()
