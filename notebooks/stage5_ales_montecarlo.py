"""
Stage 5 - ALES-Arid parametric land suitability with Monte Carlo uncertainty

Framework: Sys & Verheye (1992) / Sys, Van Ranst & Debaveye (1991) parametric
land evaluation, adapted for arid lands.

Suitability rating per pixel = (R_topo * R_wetness * R_physical * R_chemical
                                * R_fertility) / 100^4

Sub-ratings (linear interpolation between thresholds):

  R_topo (slope, %): 100 if <=2, 95 if <=4, 85 if <=8, 60 if <=16,
                     30 if <=30, 10 otherwise

  R_wetness (Ks, cm/h):  100 if 1.5-15, 85 if 0.5-1.5 or 15-30,
                         60 if 0.1-0.5 or 30-60, 30 if <0.1 or >60

  R_physical (Clay %, Sand %): texture rating from clay class
                          (loam/clay-loam/sandy-loam = 100; clay = 75;
                           sand = 50; sandy/clay extremes = 30)

  R_chemical (EC dS/m, ESP %, CaCO3 %):
        EC: 100 if <=2, 90 if <=4, 75 if <=8, 50 if <=16, 25 if <=30, 10 otherwise
        ESP: 100 if <=10, 90 if <=15, 70 if <=25, 50 if <=40, 25 otherwise
        CaCO3: 100 if <=10, 95 if <=20, 85 if <=30, 70 if <=50, 50 otherwise
        R_chemical = min(EC, ESP, CaCO3) (limiting factor)

  R_fertility (CEC, SOM):
        CEC: 100 if >=20, 90 if >=15, 75 if >=10, 60 if >=5, 40 otherwise
        SOM: 100 if >=2, 90 if >=1, 70 if >=0.5, 50 otherwise
        R_fertility = sqrt(CEC * SOM)   (geometric mean)

Class breaks on final index (Sys et al.):
  S1 (highly suitable):     >=75
  S2 (moderately suitable): 50-75
  S3 (marginally suitable): 25-50
  N1 (currently not):       12.5-25
  N2 (permanently not):     <12.5

Monte-Carlo uncertainty propagation:
  For each soil property, treat the (q05, q50, q95) trio from Stage 3b as a
  triangular distribution. Draw N=200 realisations per pixel, compute the
  suitability index N times, report q05/q50/q95 of the index AND P(class=S1+S2)
  (probability of being at least moderately suitable).

Inputs:
    outputs/stage3b/{prop}_{depth}_{q05_qrf,q50_qrf,q95_qrf}.tif
    outputs/stage2b_local_stack/covariate_stack_30m.tif (slope band)

Outputs in outputs/stage5/:
    For each depth (0_30cm, 0_100cm):
        ales_index_q05.tif, ales_index_q50.tif, ales_index_q95.tif
        ales_class_q50.tif       uint8 (1=S1..5=N2)
        p_at_least_S2.tif        float32 in [0,1]
        ales_subratings_q50.tif  6-band stack (topo, wet, phys, chem, fert, total)
    stage5_summary.txt
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import rasterio

ROOT     = Path(__file__).resolve().parents[1]
STAGE2B  = ROOT / "outputs" / "stage2b_local_stack"
STAGE3B  = ROOT / "outputs" / "stage3b"
OUT      = ROOT / "outputs" / "stage5"
OUT.mkdir(parents=True, exist_ok=True)

STACK    = STAGE2B / "covariate_stack_30m.tif"
RANDOM_STATE = 42
N_DRAWS  = 200
DEPTHS   = ["0_30cm", "0_100cm"]
PROPS    = ["EC", "ESP", "CaCO3", "Clay", "Sand_total", "CEC", "SOM", "Ks"]

CLASS_BREAKS = [(75, 100, 1, "S1"), (50, 75, 2, "S2"),
                (25, 50, 3, "S3"), (12.5, 25, 4, "N1"),
                (0,  12.5, 5, "N2")]

# ---------------------------------------------------------------------------
def piecewise(x, knots, vals):
    """Linear interpolation, flat outside the knot range."""
    return np.interp(x, knots, vals, left=vals[0], right=vals[-1])

def r_topo(slope_pct):
    return piecewise(slope_pct, [0, 2, 4, 8, 16, 30, 100],
                     [100, 100, 95, 85, 60, 30, 10])

def r_wetness(ks):
    # double-sided: optimum ~5 cm/h, falls off both sides
    out = np.where(
        ks < 1.5,
        piecewise(ks, [0, 0.1, 0.5, 1.5], [30, 30, 60, 85]),
        np.where(
            ks <= 15, 100,
            np.where(ks <= 30, 85,
                     np.where(ks <= 60, 60, 30))))
    return out.astype(np.float32)

def r_physical(clay, sand):
    # Generic textural rating (target: medium texture)
    silt = np.clip(100 - clay - sand, 0, 100)
    # distance from "ideal" loam (clay 20, sand 40, silt 40)
    d = np.sqrt((clay - 20) ** 2 + (sand - 40) ** 2 + (silt - 40) ** 2)
    return np.clip(100 - 0.9 * d, 30, 100)

def r_ec(ec):       return piecewise(ec, [0, 2, 4, 8, 16, 30, 100], [100, 100, 90, 75, 50, 25, 10])
def r_esp(esp):     return piecewise(esp, [0, 10, 15, 25, 40, 100], [100, 100, 90, 70, 50, 25])
def r_caco3(caco3): return piecewise(caco3, [0, 10, 20, 30, 50, 100], [100, 100, 95, 85, 70, 50])

def r_chemical(ec, esp, caco3):
    return np.minimum(np.minimum(r_ec(ec), r_esp(esp)), r_caco3(caco3))

def r_cec(cec): return piecewise(cec, [0, 5, 10, 15, 20, 100], [40, 40, 60, 75, 90, 100])
def r_som(som): return piecewise(som, [0, 0.5, 1, 2, 10],     [50, 50, 70, 90, 100])

def r_fertility(cec, som):
    return np.sqrt(np.clip(r_cec(cec), 1, None) * np.clip(r_som(som), 1, None))

def ales_index(slope, ks, clay, sand, ec, esp, caco3, cec, som):
    rt = r_topo(slope)
    rw = r_wetness(ks)
    rp = r_physical(clay, sand)
    rc = r_chemical(ec, esp, caco3)
    rf = r_fertility(cec, som)
    idx = (rt * rw * rp * rc * rf) / (100.0 ** 4)
    return rt, rw, rp, rc, rf, idx

def classify(idx):
    out = np.full(idx.shape, np.uint8(0))
    for lo, hi, c, _ in CLASS_BREAKS:
        m = (idx >= lo) & (idx < hi if c != 1 else idx <= hi + 1e-6)
        out[m] = c
    return out

# ---------------------------------------------------------------------------
def triangular_draw(q05, q50, q95, n, rng):
    """Sample n values per pixel from triangular(q05, q50, q95)."""
    # broadcast: out shape (n, *q50.shape)
    shp = (n,) + q50.shape
    u = rng.random(shp, dtype=np.float32)
    span = (q95 - q05).astype(np.float32)
    span = np.where(span > 0, span, 1.0)
    fc = ((q50 - q05) / span).astype(np.float32)   # mode location in [0,1]
    left = u < fc
    a = q05.astype(np.float32); b = q95.astype(np.float32)
    out = np.where(left,
                   a + np.sqrt(u * fc) * span,
                   b - np.sqrt((1 - u) * (1 - fc)) * span)
    return out

# ---------------------------------------------------------------------------
def load_qrf(prop, depth, key):
    """Load a single GeoTIFF, return (array, profile)."""
    with rasterio.open(STAGE3B / f"{prop}_{depth}_{key}.tif") as src:
        return src.read(1), src.profile

def main():
    # Load slope from stack
    with rasterio.open(STACK) as src:
        bands = list(src.descriptions) or (STAGE2B / "band_names.txt").read_text().splitlines()
        slope = src.read(bands.index("slope") + 1).astype(np.float32)
        H, W = src.height, src.width
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
        # Load all 3 quantiles for each property
        q = {}
        for p in PROPS:
            q[p] = {k: load_qrf(p, depth, k)[0].astype(np.float32)
                    for k in ("q05_qrf", "q50_qrf", "q95_qrf")}
            print(f"  loaded {p}")

        valid = np.isfinite(q[PROPS[0]]["q50_qrf"]) & np.isfinite(slope)

        # --- median surface (single deterministic run) ---
        rt, rw, rp, rc, rf, idx_med = ales_index(
            slope,
            q["Ks"]["q50_qrf"], q["Clay"]["q50_qrf"], q["Sand_total"]["q50_qrf"],
            q["EC"]["q50_qrf"], q["ESP"]["q50_qrf"], q["CaCO3"]["q50_qrf"],
            q["CEC"]["q50_qrf"], q["SOM"]["q50_qrf"],
        )
        for arr in (rt, rw, rp, rc, rf, idx_med):
            arr[~valid] = np.nan

        cls_med = classify(np.where(valid, idx_med, 0)).astype(np.uint8)
        cls_med[~valid] = 0

        # --- Monte-Carlo for q05/q95 of index and P(class<=2) ---
        # Stream draws to stay memory-friendly: keep running quantile estimators
        # via reservoir of all draws (200 * H*W * 4B = ~1 GB, acceptable here).
        print(f"  Monte Carlo {N_DRAWS} draws (slope deterministic) ...", flush=True)
        idx_draws = np.empty((N_DRAWS, H, W), dtype=np.float32)
        # Pre-sample all properties (200 x H x W x 8 = ~6.7 GB float32).
        # That's too large; instead loop draws one at a time.
        n_at_least_S2 = np.zeros((H, W), dtype=np.uint16)
        for d in range(N_DRAWS):
            samp = {}
            for p in PROPS:
                u = rng.random((H, W), dtype=np.float32)
                a = q[p]["q05_qrf"]; m = q[p]["q50_qrf"]; b = q[p]["q95_qrf"]
                span = b - a
                fc = np.where(span > 0, (m - a) / np.where(span > 0, span, 1), 0.5)
                left = u < fc
                samp[p] = np.where(
                    left,
                    a + np.sqrt(np.maximum(u * fc, 0)) * span,
                    b - np.sqrt(np.maximum((1 - u) * (1 - fc), 0)) * span,
                ).astype(np.float32)
            _, _, _, _, _, idx_d = ales_index(
                slope,
                samp["Ks"], samp["Clay"], samp["Sand_total"],
                samp["EC"], samp["ESP"], samp["CaCO3"],
                samp["CEC"], samp["SOM"],
            )
            idx_d[~valid] = np.nan
            idx_draws[d] = idx_d
            n_at_least_S2 += (idx_d >= 50).astype(np.uint16)
            if (d + 1) % 25 == 0:
                print(f"    draw {d+1}/{N_DRAWS}", flush=True)

        idx_q05 = np.nanpercentile(idx_draws, 5, axis=0).astype(np.float32)
        idx_q95 = np.nanpercentile(idx_draws, 95, axis=0).astype(np.float32)
        idx_q05[~valid] = np.nan; idx_q95[~valid] = np.nan
        p_geS2 = (n_at_least_S2.astype(np.float32) / N_DRAWS)
        p_geS2[~valid] = np.nan
        del idx_draws

        # --- write outputs ---
        outs = {
            f"ales_index_q05_{depth}.tif": (idx_q05, fmeta),
            f"ales_index_q50_{depth}.tif": (idx_med.astype(np.float32), fmeta),
            f"ales_index_q95_{depth}.tif": (idx_q95, fmeta),
            f"ales_class_q50_{depth}.tif": (cls_med, cmeta),
            f"p_at_least_S2_{depth}.tif":  (p_geS2, fmeta),
        }
        for name, (arr, m) in outs.items():
            with rasterio.open(OUT / name, "w", **m) as dst:
                dst.write(arr, 1)
            print(f"  wrote {name}")

        # 6-band sub-rating stack
        sub_meta = dict(fmeta, count=6)
        with rasterio.open(OUT / f"ales_subratings_q50_{depth}.tif", "w", **sub_meta) as dst:
            for i, (arr, name) in enumerate(zip(
                    [rt, rw, rp, rc, rf, idx_med],
                    ["topo", "wetness", "physical", "chemical",
                     "fertility", "index"])):
                dst.write(arr.astype(np.float32), i + 1)
                dst.set_band_description(i + 1, name)

        # --- per-depth summary ---
        total = int(valid.sum())
        cls_counts = {c: int(((cls_med == c) & valid).sum())
                      for _, _, c, _ in CLASS_BREAKS}
        block = [f"\n=== Depth {depth} ==="]
        block.append(f"  Valid pixels         : {total}")
        block.append(f"  Mean ALES index q50  : {np.nanmean(idx_med):.2f}")
        block.append(f"  Mean PI width q95-q05: {np.nanmean(idx_q95 - idx_q05):.2f}")
        block.append(f"  Mean P(>=S2)         : {np.nanmean(p_geS2):.3f}")
        block.append("  Class shares (deterministic q50):")
        for lo, hi, c, name in CLASS_BREAKS:
            pct = 100 * cls_counts[c] / total if total else 0
            block.append(f"    {c} {name} ({lo:>4}-{hi:>4}): "
                         f"{pct:5.1f} %  ({cls_counts[c]} px)")
        summary_blocks.append("\n".join(block))

    text = ("Stage 5 - ALES-Arid suitability with Monte-Carlo uncertainty\n"
            + "=" * 60
            + f"\nMonte-Carlo draws : {N_DRAWS} (triangular q05/q50/q95)\n"
            + f"Properties varied : {PROPS}\n"
            + f"Slope source       : DEM band (deterministic)\n"
            + "\n".join(summary_blocks))
    print("\n" + text)
    (OUT / "stage5_summary.txt").write_text(text)

if __name__ == "__main__":
    main()
