"""
Stage 16 - does the suitability result depend on which median surface feeds the index?

Stage 3b produced two point predictions per property: an inverse-distance
weighted median (the best performer under spatial cross-validation for most
properties) and the quantile regression forest median (the one that is
internally consistent with the q05/q95 envelope the Monte Carlo requires). The
index and everything downstream use the QRF median. This script recomputes the
baseline index and the counterfactual decomposition from the IDW median instead,
so that the consequence of that choice is reported rather than assumed.
"""
from __future__ import annotations
import sys, pathlib
import numpy as np, pandas as pd, rasterio

import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parent))
from revpaths import REPO as _REPO, OUT as _OUT

REPO = _REPO
OUT = _OUT
sys.path.insert(0, str(REPO / "notebooks"))
from stage5_ales_montecarlo import (ales_index, classify, STAGE3B, STACK, STAGE2B,
                                    PROPS, DEPTHS)

SCEN = {"CF_EC": {"EC": (2.0, "cap")}, "CF_ESP": {"ESP": (10.0, "cap")},
        "CF_CaCO3": {"CaCO3": (15.0, "cap")}, "CF_Ks": {"Ks": (2.0, "lift")},
        "CF_SOM": {"SOM": (1.0, "lift")},
        "CF_full": {"EC": (2.0, "cap"), "ESP": (10.0, "cap"), "CaCO3": (15.0, "cap"),
                    "Ks": (2.0, "lift"), "SOM": (1.0, "lift")}}


def rd(p):
    with rasterio.open(p) as s:
        a = s.read(1).astype("float32")
        if s.nodata is not None and np.isfinite(s.nodata):
            a[a == s.nodata] = np.nan
        return a


def main():
    with rasterio.open(STACK) as s:
        names = list(s.descriptions) or (STAGE2B / "band_names.txt").read_text().splitlines()
        slope = s.read(names.index("slope") + 1).astype("float32")
    rows, lim_rows = [], []
    for depth in DEPTHS:
        for src, key in [("QRF median (as used)", "q50_qrf"), ("IDW median", "q50_idw")]:
            v = {p: rd(STAGE3B / f"{p}_{depth}_{key}.tif") for p in PROPS}
            valid = np.isfinite(v["Ks"]) & np.isfinite(slope)
            rt, rw, rp, rc, rf, idx = ales_index(
                slope, v["Ks"], v["Clay"], v["Sand_total"], v["EC"], v["ESP"],
                v["CaCO3"], v["CEC"], v["SOM"])
            cls = classify(np.where(valid, idx, 0))
            shares = {nm: 100 * float(((cls == c) & valid).sum()) / valid.sum()
                      for c, nm in [(1, "S1"), (2, "S2"), (3, "S3"), (4, "N1"), (5, "N2")]}
            S = np.stack([rt, rw, rp, rc, rf])
            amin = np.argmin(np.where(np.isfinite(S), S, np.inf), axis=0)
            names5 = ["topo", "wetness", "physical", "chemical", "fertility"]
            for i, nm in enumerate(names5):
                lim_rows.append(dict(depth=depth, source=src, sub_rating=nm,
                                     pct_limiting=round(100 * float(((amin == i) & valid).sum())
                                                        / valid.sum(), 2)))
            row = dict(depth=depth, source=src,
                       mean_index=round(float(np.nanmean(idx[valid])), 2), **{k: round(x, 1) for k, x in shares.items()})
            for sc, spec in SCEN.items():
                vv = dict(v)
                for pr, (tgt, mode) in spec.items():
                    vv[pr] = np.minimum(v[pr], tgt) if mode == "cap" else np.maximum(v[pr], tgt)
                _, _, _, _, _, i2 = ales_index(slope, vv["Ks"], vv["Clay"], vv["Sand_total"],
                                               vv["EC"], vv["ESP"], vv["CaCO3"], vv["CEC"],
                                               vv["SOM"])
                row[sc] = round(float(np.nanmean((i2 - idx)[valid])), 2)
            row["CF_Ks : CF_EC"] = round(row["CF_Ks"] / row["CF_EC"], 1) if row["CF_EC"] else np.inf
            rows.append(row)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "idw_vs_qrf_index.csv", index=False)
    ld = pd.DataFrame(lim_rows)
    ld.to_csv(OUT / "idw_vs_qrf_limiting.csv", index=False)
    print("===== baseline and counterfactuals from each median surface =====")
    print(df.to_string(index=False))
    print("\n===== limiting sub-rating (% of pixels) =====")
    print(ld.pivot_table(index=["depth", "sub_rating"], columns="source",
                         values="pct_limiting").round(2).to_string())


if __name__ == "__main__":
    main()
