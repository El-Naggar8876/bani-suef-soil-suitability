"""
Stage 6 - Counterfactual gap decomposition (THE HEADLINE NOVELTY)

Question: which soil constraints, if removed, would most close the gap
between current ALES-Arid suitability and full S1 suitability?

Method:
  Baseline                  : current ALES q50 surface from Stage 5
  Counterfactual scenarios  : replace one limiting factor with a "best-case"
                              level, recompute ALES, compare to baseline.
  Attribution               : Delta_index per scenario, Delta_class share,
                              pixels flipped UP (e.g. N1 -> S3).

Scenarios (per depth):
  CF_EC      : EC capped at 2 dS/m            (full leaching with good water)
  CF_ESP     : ESP capped at 10 %             (gypsum + leaching)
  CF_CaCO3   : CaCO3 capped at 15 %           (deep ploughing / dilution)
  CF_Ks      : Ks raised to >= 2 cm/h         (subsoiling / sand amendment)
  CF_SOM     : SOM raised to >= 1 %           (compost / cover crops)
  CF_chem    : EC + ESP + CaCO3 simultaneously
  CF_full    : all five interventions

The cap = min(baseline, target) for "lower is better" properties (EC, ESP,
CaCO3) and = max(baseline, target) for "higher is better" (Ks, SOM).

We also do this with Monte-Carlo (N draws) to attach uncertainty to each
intervention's effect.

Outputs in outputs/stage6/:
    Per (depth, scenario):
      cf_{scenario}_{depth}_index_q50.tif
      cf_{scenario}_{depth}_class.tif
      cf_{scenario}_{depth}_delta_q50.tif        (idx_cf - idx_base)
      cf_{scenario}_{depth}_p_S2plus.tif         (Monte-Carlo)
    Per depth:
      gap_decomposition_{depth}.csv              (mean delta by scenario)
      class_transition_{depth}_{scenario}.csv    (transition matrices)
    stage6_summary.txt
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
import rasterio

# Re-use Stage 5 helpers
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from stage5_ales_montecarlo import (
    ales_index, classify, CLASS_BREAKS, PROPS, DEPTHS,
    STACK, STAGE3B, STAGE2B, load_qrf,
)

ROOT = Path(__file__).resolve().parents[1]
OUT  = ROOT / "outputs" / "stage6"
OUT.mkdir(parents=True, exist_ok=True)

RANDOM_STATE = 42
N_DRAWS = 100   # smaller than Stage 5 (each scenario re-runs MC)

# Scenario definitions: dict of property -> (target_value, mode)
# mode = "cap"  -> new = min(baseline, target)   (lower better)
# mode = "lift" -> new = max(baseline, target)   (higher better)
SCENARIOS = {
    "CF_EC":     {"EC":    (2.0,  "cap")},
    "CF_ESP":    {"ESP":   (10.0, "cap")},
    "CF_CaCO3":  {"CaCO3": (15.0, "cap")},
    "CF_Ks":     {"Ks":    (2.0,  "lift")},
    "CF_SOM":    {"SOM":   (1.0,  "lift")},
    "CF_chem":   {"EC":    (2.0,  "cap"),
                  "ESP":   (10.0, "cap"),
                  "CaCO3": (15.0, "cap")},
    "CF_full":   {"EC":    (2.0,  "cap"),
                  "ESP":   (10.0, "cap"),
                  "CaCO3": (15.0, "cap"),
                  "Ks":    (2.0,  "lift"),
                  "SOM":   (1.0,  "lift")},
}

# ---------------------------------------------------------------------------
def apply_scenario(values, scenario):
    """Return shallow-copied dict with scenario substitutions."""
    out = dict(values)
    for prop, (target, mode) in scenario.items():
        if mode == "cap":
            out[prop] = np.minimum(values[prop], target).astype(np.float32)
        elif mode == "lift":
            out[prop] = np.maximum(values[prop], target).astype(np.float32)
    return out

def compute_index(values, slope):
    return ales_index(
        slope, values["Ks"], values["Clay"], values["Sand_total"],
        values["EC"], values["ESP"], values["CaCO3"],
        values["CEC"], values["SOM"],
    )[-1]

# ---------------------------------------------------------------------------
def main():
    with rasterio.open(STACK) as src:
        bands = list(src.descriptions) or (STAGE2B / "band_names.txt").read_text().splitlines()
        slope = src.read(bands.index("slope") + 1).astype(np.float32)
        H, W  = src.height, src.width
        T, crs = src.transform, src.crs

    fmeta = dict(driver="GTiff", height=H, width=W, count=1, crs=crs,
                 transform=T, dtype="float32", nodata=np.float32(np.nan),
                 compress="deflate", predictor=3, tiled=True,
                 blockxsize=512, blockysize=512, BIGTIFF="IF_NEEDED")
    cmeta = dict(fmeta, dtype="uint8", nodata=0, predictor=2)

    rng = np.random.default_rng(RANDOM_STATE)
    summary_blocks = []

    for depth in DEPTHS:
        print(f"\n=== Depth {depth} ===", flush=True)
        # Load q05/q50/q95 for every property
        q = {}
        for p in PROPS:
            q[p] = {k: load_qrf(p, depth, k)[0].astype(np.float32)
                    for k in ("q05_qrf", "q50_qrf", "q95_qrf")}

        valid = np.isfinite(q[PROPS[0]]["q50_qrf"]) & np.isfinite(slope)

        # Baseline (deterministic)
        base_vals = {p: q[p]["q50_qrf"] for p in PROPS}
        base_idx  = compute_index(base_vals, slope)
        base_idx[~valid] = np.nan
        base_class = classify(np.where(valid, base_idx, 0))
        base_class[~valid] = 0

        deco_rows = []
        for sname, sdef in SCENARIOS.items():
            print(f"  -> {sname}", flush=True)
            # Deterministic q50 counterfactual
            cf_vals = apply_scenario(base_vals, sdef)
            cf_idx  = compute_index(cf_vals, slope)
            cf_idx[~valid] = np.nan
            delta = (cf_idx - base_idx).astype(np.float32)
            cf_class = classify(np.where(valid, cf_idx, 0))
            cf_class[~valid] = 0

            # Monte-Carlo P(>=S2) under counterfactual
            n_geS2 = np.zeros((H, W), dtype=np.uint16)
            for d in range(N_DRAWS):
                samp = {}
                for p in PROPS:
                    a = q[p]["q05_qrf"]; m = q[p]["q50_qrf"]; b = q[p]["q95_qrf"]
                    span = b - a
                    fc = np.where(span > 0, (m - a) / np.where(span > 0, span, 1), 0.5)
                    u = rng.random((H, W), dtype=np.float32)
                    left = u < fc
                    samp[p] = np.where(
                        left,
                        a + np.sqrt(np.maximum(u * fc, 0)) * span,
                        b - np.sqrt(np.maximum((1 - u) * (1 - fc), 0)) * span,
                    ).astype(np.float32)
                samp = apply_scenario(samp, sdef)
                idx_d = compute_index(samp, slope)
                idx_d[~valid] = np.nan
                n_geS2 += (idx_d >= 50).astype(np.uint16)
            p_geS2 = (n_geS2.astype(np.float32) / N_DRAWS)
            p_geS2[~valid] = np.nan

            # Write rasters
            for name, (arr, m) in {
                f"cf_{sname}_{depth}_index_q50.tif": (cf_idx.astype(np.float32), fmeta),
                f"cf_{sname}_{depth}_delta_q50.tif": (delta, fmeta),
                f"cf_{sname}_{depth}_class.tif":     (cf_class, cmeta),
                f"cf_{sname}_{depth}_p_S2plus.tif":  (p_geS2, fmeta),
            }.items():
                with rasterio.open(OUT / name, "w", **m) as dst:
                    dst.write(arr, 1)

            # Class transition matrix (5x5)
            tm = np.zeros((5, 5), dtype=np.int64)
            mask = valid & (base_class > 0) & (cf_class > 0)
            for i in range(1, 6):
                for j in range(1, 6):
                    tm[i-1, j-1] = int(((base_class == i) & (cf_class == j) & mask).sum())
            tm_df = pd.DataFrame(
                tm, index=[c for _, _, c, _ in CLASS_BREAKS],
                columns=[c for _, _, c, _ in CLASS_BREAKS],
            )
            tm_df.index.name = "from"; tm_df.columns.name = "to"
            tm_df.to_csv(OUT / f"class_transition_{depth}_{sname}.csv")

            mean_delta   = float(np.nanmean(delta))
            mean_p_geS2  = float(np.nanmean(p_geS2))
            cf_class_share = {c: float(((cf_class == c) & valid).sum() / valid.sum())
                              for _, _, c, _ in CLASS_BREAKS}
            n_flip_up = int(((cf_class > 0) & (base_class > 0)
                             & (cf_class < base_class)).sum())
            # class number lower = better, so cf_class < base_class means flip up

            deco_rows.append(dict(
                scenario=sname, depth=depth,
                mean_delta_index=round(mean_delta, 3),
                mean_p_at_least_S2=round(mean_p_geS2, 3),
                share_S1=round(cf_class_share[1], 4),
                share_S2=round(cf_class_share[2], 4),
                share_S3=round(cf_class_share[3], 4),
                share_N1=round(cf_class_share[4], 4),
                share_N2=round(cf_class_share[5], 4),
                pixels_improved_class=n_flip_up,
            ))

        deco_df = pd.DataFrame(deco_rows)
        deco_df.to_csv(OUT / f"gap_decomposition_{depth}.csv", index=False)

        # Compute baseline shares for context
        base_share = {c: float(((base_class == c) & valid).sum() / valid.sum())
                      for _, _, c, _ in CLASS_BREAKS}
        block = [f"\n=== Depth {depth} ==="]
        block.append(f"  Baseline mean ALES q50 : {np.nanmean(base_idx):.2f}")
        block.append(f"  Baseline class shares  : "
                     f"S1={base_share[1]*100:.1f}% S2={base_share[2]*100:.1f}% "
                     f"S3={base_share[3]*100:.1f}% N1={base_share[4]*100:.1f}% "
                     f"N2={base_share[5]*100:.1f}%")
        block.append("  Counterfactual mean DELTA index (rank order):")
        for r in sorted(deco_rows, key=lambda x: -x["mean_delta_index"]):
            block.append(f"    {r['scenario']:10s}  "
                         f"Delta={r['mean_delta_index']:+6.2f}  "
                         f"P(>=S2)={r['mean_p_at_least_S2']:.3f}  "
                         f"S2={r['share_S2']*100:5.1f}%  "
                         f"S3={r['share_S3']*100:5.1f}%  "
                         f"flipped_up={r['pixels_improved_class']}")
        summary_blocks.append("\n".join(block))

    text = ("Stage 6 - Counterfactual gap decomposition\n"
            + "=" * 60 + "\n"
            + f"Monte-Carlo draws per scenario : {N_DRAWS}\n"
            + f"Scenarios                      : {list(SCENARIOS)}\n"
            + "Targets (cap=upper bound for bad properties; "
            + "lift=lower bound for good properties):\n"
            + "  EC<=2 dS/m, ESP<=10%, CaCO3<=15%, Ks>=2 cm/h, SOM>=1%\n"
            + "\n".join(summary_blocks))
    print("\n" + text)
    (OUT / "stage6_summary.txt").write_text(text)

if __name__ == "__main__":
    main()
