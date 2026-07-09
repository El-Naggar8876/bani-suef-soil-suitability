"""
Stage 7 - Cropping mismatch & policy priority

Use Sentinel-2 winter+summer NDVI/MNDWI from the Stage 2b stack to label each
pixel as one of:
    1 = water        (MNDWI_sum > 0.30 OR MNDWI_win > 0.30)
    2 = double crop  (NDVI_win > 0.40 AND NDVI_sum > 0.40)
    3 = single crop  (only one season > 0.40)
    4 = bare/desert  (both seasons < 0.20)
    5 = sparse veg   (else)

Then cross-tabulate with Stage 5 ALES-class (q50) per depth and produce a
"policy priority" map = pixels that are CURRENTLY CROPPED (class 2 or 3) AND
sit on N1 substrate AND show large CF_Ks delta (top quartile of delta in the
Ks counterfactual).

Outputs in outputs/stage7/:
    landuse_class.tif                       uint8 (1..5)
    cropping_x_ales_{depth}.csv             contingency tables
    policy_priority_{depth}.tif             uint8 (1=high, 2=medium, 0=other)
    suitability_realised_{depth}.csv        per ALES class: cropped fraction
    stage7_summary.txt
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
import rasterio

ROOT     = Path(__file__).resolve().parents[1]
STAGE2B  = ROOT / "outputs" / "stage2b_local_stack"
STAGE5   = ROOT / "outputs" / "stage5"
STAGE6   = ROOT / "outputs" / "stage6"
OUT      = ROOT / "outputs" / "stage7"
OUT.mkdir(parents=True, exist_ok=True)

STACK = STAGE2B / "covariate_stack_30m.tif"
DEPTHS = ["0_30cm", "0_100cm"]

LU_LABELS = {1: "water", 2: "double_crop", 3: "single_crop",
             4: "bare_desert", 5: "sparse_veg"}
ALES_LABELS = {1: "S1", 2: "S2", 3: "S3", 4: "N1", 5: "N2"}

def classify_landuse(ndvi_win, ndvi_sum, mndwi_win, mndwi_sum):
    out = np.full(ndvi_win.shape, 5, dtype=np.uint8)   # default sparse
    water  = (mndwi_sum > 0.30) | (mndwi_win > 0.30)
    bare   = (ndvi_win < 0.20) & (ndvi_sum < 0.20)
    double = (ndvi_win > 0.40) & (ndvi_sum > 0.40)
    single = ((ndvi_win > 0.40) ^ (ndvi_sum > 0.40))   # XOR
    out[single] = 3
    out[double] = 2
    out[bare]   = 4
    out[water]  = 1
    return out

def main():
    with rasterio.open(STACK) as src:
        bands = list(src.descriptions) or (STAGE2B / "band_names.txt").read_text().splitlines()
        ndvi_win  = src.read(bands.index("NDVI_win") + 1).astype(np.float32)
        ndvi_sum  = src.read(bands.index("NDVI_sum") + 1).astype(np.float32)
        mndwi_win = src.read(bands.index("MNDWI_win") + 1).astype(np.float32)
        mndwi_sum = src.read(bands.index("MNDWI_sum") + 1).astype(np.float32)
        H, W = src.height, src.width
        T, crs = src.transform, src.crs

    valid_stack = np.isfinite(ndvi_win) & np.isfinite(ndvi_sum)
    lu = classify_landuse(np.nan_to_num(ndvi_win, nan=0),
                          np.nan_to_num(ndvi_sum, nan=0),
                          np.nan_to_num(mndwi_win, nan=0),
                          np.nan_to_num(mndwi_sum, nan=0))
    lu[~valid_stack] = 0

    fmeta = dict(driver="GTiff", height=H, width=W, count=1, crs=crs,
                 transform=T, compress="deflate", predictor=2, tiled=True,
                 blockxsize=512, blockysize=512, BIGTIFF="IF_NEEDED",
                 dtype="uint8", nodata=0)
    with rasterio.open(OUT / "landuse_class.tif", "w", **fmeta) as dst:
        dst.write(lu, 1)
        dst.update_tags(**{f"CLASS_{k}": v for k, v in LU_LABELS.items()})
    print(f"Wrote landuse_class.tif")

    # Overall LU shares
    total = int(valid_stack.sum())
    lu_shares = {LU_LABELS[k]: 100 * int((lu == k).sum()) / total
                 for k in LU_LABELS}
    print("Land use shares (%):", {k: round(v, 1) for k, v in lu_shares.items()})

    summary_blocks = []
    for depth in DEPTHS:
        with rasterio.open(STAGE5 / f"ales_class_q50_{depth}.tif") as r:
            ales = r.read(1)
        with rasterio.open(STAGE6 / f"cf_CF_Ks_{depth}_delta_q50.tif") as r:
            cf_delta = r.read(1).astype(np.float32)

        mask = valid_stack & (ales > 0) & np.isfinite(cf_delta)

        # Contingency: rows = LU, cols = ALES
        ct = np.zeros((5, 5), dtype=np.int64)
        for i in range(1, 6):
            for j in range(1, 6):
                ct[i-1, j-1] = int(((lu == i) & (ales == j) & mask).sum())
        ct_df = pd.DataFrame(ct,
                             index=[LU_LABELS[i] for i in range(1, 6)],
                             columns=[ALES_LABELS[j] for j in range(1, 6)])
        ct_df.index.name = "land_use"; ct_df.columns.name = "ales_class"
        ct_df.to_csv(OUT / f"cropping_x_ales_{depth}.csv")

        # Suitability realisation: per ALES class, what fraction is cropped?
        rows = []
        cropped = (lu == 2) | (lu == 3)
        for j in range(1, 6):
            cls_mask = (ales == j) & mask
            n_cls = int(cls_mask.sum())
            n_crop = int((cls_mask & cropped).sum())
            rows.append(dict(
                ales_class=ALES_LABELS[j],
                pixels=n_cls,
                cropped_pixels=n_crop,
                cropped_fraction=round(n_crop / n_cls if n_cls else 0, 3),
            ))
        sr_df = pd.DataFrame(rows)
        sr_df.to_csv(OUT / f"suitability_realised_{depth}.csv", index=False)

        # Policy priority: cropped AND N1 substrate AND large CF_Ks delta
        delta_top_q = np.nanpercentile(cf_delta[mask], 75)
        priority = np.zeros((H, W), dtype=np.uint8)
        # 1 = HIGH:   cropped, N1, delta >= top quartile
        # 2 = MEDIUM: cropped, N1, delta < top quartile
        # 3 = WATCH:  cropped on S3 with delta >= top quartile (improvable)
        is_n1 = (ales == 4)
        is_s3 = (ales == 3)
        priority[mask & cropped & is_s3 & (cf_delta >= delta_top_q)] = 3
        priority[mask & cropped & is_n1] = 2
        priority[mask & cropped & is_n1 & (cf_delta >= delta_top_q)] = 1

        with rasterio.open(OUT / f"policy_priority_{depth}.tif", "w", **fmeta) as dst:
            dst.write(priority, 1)
            dst.update_tags(CLASS_1="HIGH_cropped_N1_largeKsGain",
                            CLASS_2="MEDIUM_cropped_N1",
                            CLASS_3="WATCH_cropped_S3_largeKsGain")

        n_high = int((priority == 1).sum())
        n_med  = int((priority == 2).sum())
        n_watch = int((priority == 3).sum())
        n_cropped_total = int(((cropped) & mask).sum())

        block = [f"\n=== Depth {depth} ==="]
        block.append("  Land-use x ALES contingency (pixels):")
        block.append(ct_df.to_string())
        block.append("\n  Suitability realisation (cropped fraction per ALES):")
        block.append(sr_df.to_string(index=False))
        block.append(f"\n  CF_Ks delta top-quartile threshold: {delta_top_q:.2f}")
        block.append(f"  Policy priority pixels:")
        block.append(f"    HIGH   (cropped on N1, large Ks gain) : {n_high} "
                     f"= {100*n_high/n_cropped_total:.1f} % of cropped")
        block.append(f"    MEDIUM (cropped on N1, smaller gain)  : {n_med} "
                     f"= {100*n_med/n_cropped_total:.1f} % of cropped")
        block.append(f"    WATCH  (cropped on S3, large Ks gain) : {n_watch} "
                     f"= {100*n_watch/n_cropped_total:.1f} % of cropped")
        summary_blocks.append("\n".join(block))

    text = ("Stage 7 - Cropping mismatch & policy priority\n"
            + "=" * 60 + "\n"
            + f"Land-use thresholds  : NDVI 0.40 (crop), MNDWI 0.30 (water), "
              f"NDVI 0.20 (bare)\n"
            + f"Land-use shares (%): "
            + ", ".join(f"{k}={v:.1f}" for k, v in lu_shares.items())
            + "\n"
            + "\n".join(summary_blocks))
    print("\n" + text)
    (OUT / "stage7_summary.txt").write_text(text)

if __name__ == "__main__":
    main()
