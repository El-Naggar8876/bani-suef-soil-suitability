"""
Stage 11 - Monte Carlo propagation with an explicitly correlated draw structure.

Purpose (Reviewer 1, comment 4): the submitted implementation samples each soil
property independently from its own triangular predictive distribution. Here the
empirical dependence between the eight properties is measured from the 60
profiles, the propagation is repeated under a Gaussian copula reproducing that
dependence, and the effect on the published quantities is reported.

It also re-evaluates every counterfactual scenario *inside* the Monte Carlo loop,
so that each scenario is compared with its own draw-matched baseline. This is a
different quantity from the deterministic Delta reported in the submitted
manuscript, which perturbs the median-property surface only.

Usage:  python stage11_copula_mc.py <depth> <independent|correlated> [n_draws] [n_sample]
Writes: outputs_rev/mc_run_<depth>_<structure>.json  and the correlation matrix.
"""
from __future__ import annotations
import sys, json, pathlib
import numpy as np, pandas as pd, rasterio
from scipy.stats import spearmanr, norm

import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parent))
from revpaths import REPO as _REPO, OUT as _OUT, WAPOR as _WAPOR

HERE = pathlib.Path(__file__).resolve().parent
REPO = _REPO
sys.path.insert(0, str(REPO / "notebooks"))
OUT = _OUT
from stage5_ales_montecarlo import ales_index, STAGE3B, STACK, STAGE2B, PROPS

SEED = 42
CSV = {d: REPO / f"outputs/stage1/soil_per_profile_{d}.csv" for d in ("0_30cm", "0_100cm")}


def latent_corr(depth):
    d = pd.read_csv(CSV[depth])
    d["Sand_total"] = d["Coarse_sand"] + d["Fine_sand"]
    X = d[list(PROPS)].astype(float).values
    rho, _ = spearmanr(X, nan_policy="omit")
    rho = np.asarray(rho, float)
    R = 2 * np.sin(np.pi * rho / 6)              # Spearman -> latent Pearson
    w, V = np.linalg.eigh((R + R.T) / 2)         # nearest positive definite
    w = np.clip(w, 1e-8, None)
    A = V @ np.diag(w) @ V.T
    s = np.sqrt(np.diag(A))
    return rho, A / np.outer(s, s)


def tri_icdf(u, a, m, b):
    span = np.where(b - a > 0, b - a, 1.0)
    fc = np.clip((m - a) / span, 0.0, 1.0)
    return np.where(u < fc,
                    a + np.sqrt(np.maximum(u * fc, 0)) * span,
                    b - np.sqrt(np.maximum((1 - u) * (1 - fc), 0)) * span).astype("float32")


SCEN = {
    "CF_EC":       {"EC": (2.0, "cap")},
    "CF_ESP":      {"ESP": (10.0, "cap")},
    "CF_CaCO3":    {"CaCO3": (15.0, "cap")},
    "CF_CaCO3_10": {"CaCO3": (10.0, "cap")},
    "CF_Ks":       {"Ks": (2.0, "lift")},
    "CF_SOM":      {"SOM": (1.0, "lift")},
    "CF_full":     {"EC": (2.0, "cap"), "ESP": (10.0, "cap"), "CaCO3": (15.0, "cap"),
                    "Ks": (2.0, "lift"), "SOM": (1.0, "lift")},
}


def main():
    depth, structure = sys.argv[1], sys.argv[2]
    n_draws = int(sys.argv[3]) if len(sys.argv) > 3 else 200
    n_sample = int(sys.argv[4]) if len(sys.argv) > 4 else 150000
    corr = structure == "correlated"

    with rasterio.open(STACK) as s:
        names = list(s.descriptions) or (STAGE2B / "band_names.txt").read_text().splitlines()
        slope_full = s.read(names.index("slope") + 1).astype("float32")
    q = {}
    for p in PROPS:
        q[p] = {}
        for k in ("q05_qrf", "q50_qrf", "q95_qrf"):
            with rasterio.open(STAGE3B / f"{p}_{depth}_{k}.tif") as s:
                q[p][k] = s.read(1).astype("float32")
    valid = np.isfinite(q[PROPS[0]]["q50_qrf"]) & np.isfinite(slope_full)
    N_all = int(valid.sum())
    sel_rng = np.random.default_rng(7)
    if n_sample and n_sample < N_all:
        pick = np.zeros(N_all, bool)
        pick[sel_rng.choice(N_all, n_sample, replace=False)] = True
    else:
        pick = np.ones(N_all, bool)
    N = int(pick.sum())
    slope = slope_full[valid][pick]
    QA = {p: {k: q[p][k][valid][pick] for k in q[p]} for p in PROPS}
    del q, slope_full

    rho, R = latent_corr(depth)
    L = np.linalg.cholesky(R).astype("float32")
    if not pathlib.Path(OUT / "mc_property_spearman.csv").exists():
        rows = [dict(prop_a=PROPS[i], prop_b=PROPS[j], spearman=round(float(rho[i, j]), 3))
                for i in range(len(PROPS)) for j in range(len(PROPS)) if j > i]
        pd.DataFrame(rows).to_csv(OUT / f"mc_property_spearman_{depth}.csv", index=False)

    def idx_of(s):
        return ales_index(slope, s["Ks"], s["Clay"], s["Sand_total"], s["EC"],
                          s["ESP"], s["CaCO3"], s["CEC"], s["SOM"])[5]

    rng = np.random.default_rng(SEED)
    k = len(PROPS)
    draws = np.empty((n_draws, N), "float16")
    nge = np.zeros(N, "uint16")
    nge_cf = {sc: np.zeros(N, "uint16") for sc in SCEN}
    acc = {sc: 0.0 for sc in SCEN}
    draw_means = {sc: [] for sc in SCEN}
    for d in range(n_draws):
        e = rng.standard_normal((k, N)).astype("float32")
        z = (L @ e) if corr else e
        u = norm.cdf(z).astype("float32")
        s = {p: tri_icdf(u[i], QA[p]["q05_qrf"], QA[p]["q50_qrf"], QA[p]["q95_qrf"])
             for i, p in enumerate(PROPS)}
        idx = idx_of(s)
        draws[d] = idx.astype("float16")
        nge += (idx >= 50).astype("uint16")
        for sc, spec in SCEN.items():
            s2 = dict(s)
            for pr, (tgt, mode) in spec.items():
                s2[pr] = np.minimum(s[pr], tgt) if mode == "cap" else np.maximum(s[pr], tgt)
            icf = idx_of(s2)
            dl = icf - idx
            acc[sc] += float(dl.mean())
            draw_means[sc].append(float(dl.mean()))
            nge_cf[sc] += (icf >= 50).astype("uint16")
    D = draws.astype("float32")
    q05, q50, q95 = (np.percentile(D, p, axis=0) for p in (5, 50, 95))
    scen_out = {}
    for sc in SCEN:
        dm = np.array(draw_means[sc])
        scen_out[sc] = dict(mean_delta=round(acc[sc] / n_draws, 3),
                            ci90_lo=round(float(np.percentile(dm, 5)), 3),
                            ci90_hi=round(float(np.percentile(dm, 95)), 3),
                            P_ge_S2=round(float((nge_cf[sc] / n_draws).mean()), 4))
    res = dict(depth=depth, structure=structure, n_draws=n_draws,
               n_pixels_sampled=N, n_pixels_total=N_all,
               mean_index_mc_median=round(float(q50.mean()), 3),
               mean_PI_width=round(float((q95 - q05).mean()), 3),
               baseline_P_ge_S2=round(float((nge / n_draws).mean()), 4),
               scenarios=scen_out)
    json.dump(res, open(OUT / f"mc_run_{depth}_{structure}.json", "w"), indent=1)
    print(json.dumps(res, indent=1))

if __name__ == "__main__":
    main()
