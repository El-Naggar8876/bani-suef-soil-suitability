"""
Stage 7b - Independent crop-survey validation (n = 993, March 2025).

Cross-walks the 993 ground-observed crop points (Crop_March2025_11.shp) against:
  - Sentinel-2 land-use class    (outputs/stage7/landuse_class.tif)
  - ALES q50 index               (outputs/stage5/ales_index_q50_0_100cm.tif)
  - ALES modal class             (outputs/stage5/ales_class_q50_0_100cm.tif)
  - CF_Ks delta q50              (outputs/stage6/cf_CF_Ks_0_100cm_delta_q50.tif)
  - Policy priority              (outputs/stage7/policy_priority_0_100cm.tif)

Outputs (outputs/stage7b/):
  crop_points_extracted.csv
  crop_x_landuse_confusion.csv  + accuracy / kappa scalar
  crop_x_ales_class_0_100cm.csv
  crop_x_priority_0_100cm.csv
  crop_kpis.csv                 (per-crop n, mean q50, %N1+N2, %HIGH-priority, mean dKs)
  stage7b_summary.txt
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio

ROOT = Path(__file__).resolve().parents[1]
SHP  = ROOT / "Crop_March2025_11" / "Crop_March2025_11.shp"
S5   = ROOT / "outputs" / "stage5"
S6   = ROOT / "outputs" / "stage6"
S7   = ROOT / "outputs" / "stage7"
OUT  = ROOT / "outputs" / "stage7b"
OUT.mkdir(parents=True, exist_ok=True)

# Land-use legend (from stage7)
LU_LABELS = {1: "water", 2: "double crop", 3: "single crop",
             4: "bare/desert", 5: "sparse veg."}
ALES_LABELS = {1: "S1", 2: "S2", 3: "S3", 4: "N1", 5: "N2"}
PRI_LABELS  = {1: "HIGH", 2: "MEDIUM", 3: "WATCH"}

# Crop -> "vegetated cropland" mapping for the LU confusion matrix.
# All 15 crops are field/perennial vegetation in March 2025 (winter season:
# wheat, clover, sugar beet, garlic, onion, beans, peas, potatoes, cabbage,
# tomato + perennials citrus/banana/grape/sugar cane + medical plants).
# The S2 LU mask collapses these into "double crop" or "single crop"; both are
# correct cropland classes. "water", "bare/desert" and "sparse veg." are
# misclassifications relative to the ground truth.
CROP_TRUE_CROPLAND = True  # every survey point is cropland


def extract_at_points(gdf, raster_paths: dict) -> pd.DataFrame:
    """Sample each raster at every point. Returns a DataFrame with one column per raster."""
    df = pd.DataFrame({"CROP": gdf["CROP"].values})
    df["x"] = gdf.geometry.x.values
    df["y"] = gdf.geometry.y.values
    coords = list(zip(df["x"], df["y"]))
    for name, path in raster_paths.items():
        with rasterio.open(path) as r:
            vals = np.array([v[0] for v in r.sample(coords)], dtype=np.float64)
            nd = r.nodata
            if nd is not None:
                vals[vals == nd] = np.nan
        df[name] = vals
    return df


def cohen_kappa_binary(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    """Two-class kappa: cropland (1) vs non-cropland (0)."""
    n = len(y_true)
    po = (y_true == y_pred).mean()
    p1t = (y_true == 1).mean();  p0t = 1 - p1t
    p1p = (y_pred == 1).mean();  p0p = 1 - p1p
    pe = p1t * p1p + p0t * p0p
    if pe >= 1.0:
        return float("nan")
    return float((po - pe) / (1 - pe))


def main():
    print(f"Loading {SHP.name} ...")
    gdf = gpd.read_file(SHP, engine="pyogrio")
    print(f"  N = {len(gdf)}, CRS = {gdf.crs}")

    # AOI rasters are EPSG:32636 (matches shapefile) - confirm
    with rasterio.open(S5 / "ales_index_q50_0_100cm.tif") as r:
        if str(gdf.crs).lower() != str(r.crs).lower():
            print(f"  Re-projecting points {gdf.crs} -> {r.crs}")
            gdf = gdf.to_crs(r.crs)

    rasters = {
        "landuse":   S7 / "landuse_class.tif",
        "ales_idx":  S5 / "ales_index_q50_0_100cm.tif",
        "ales_cls":  S5 / "ales_class_q50_0_100cm.tif",
        "cf_Ks_d":   S6 / "cf_CF_Ks_0_100cm_delta_q50.tif",
        "priority":  S7 / "policy_priority_0_100cm.tif",
    }
    df = extract_at_points(gdf, rasters)
    df.to_csv(OUT / "crop_points_extracted.csv", index=False)
    print(f"  Wrote crop_points_extracted.csv ({len(df)} rows)")

    # ----- 1. Land-use confusion ---------------------------------------------
    df_lu = df.dropna(subset=["landuse"]).copy()
    df_lu["lu_label"] = df_lu["landuse"].astype(int).map(LU_LABELS)
    # Cropland = double crop OR single crop
    df_lu["pred_cropland"] = df_lu["landuse"].astype(int).isin([2, 3]).astype(int)
    df_lu["true_cropland"] = 1  # all survey points are cropland
    n = len(df_lu)
    n_correct = int(df_lu["pred_cropland"].sum())
    oa = n_correct / n
    kappa = cohen_kappa_binary(df_lu["true_cropland"].values,
                                df_lu["pred_cropland"].values)
    # Per-class breakdown of the predicted LU at survey points
    conf = (df_lu.groupby("lu_label").size()
                  .reindex(list(LU_LABELS.values()), fill_value=0)
                  .rename("n_points").to_frame())
    conf["share_%"] = (100 * conf["n_points"] / conf["n_points"].sum()).round(2)
    conf.to_csv(OUT / "crop_x_landuse_confusion.csv")
    print(f"  Sentinel-2 cropland mask: user's accuracy = {oa:.3f}  kappa = {kappa:.3f}")
    print(conf)

    # ----- 2. Crop x ALES class contingency ----------------------------------
    df_a = df.dropna(subset=["ales_cls"]).copy()
    df_a["ales_label"] = df_a["ales_cls"].astype(int).map(ALES_LABELS)
    ct = pd.crosstab(df_a["CROP"], df_a["ales_label"]) \
            .reindex(columns=list(ALES_LABELS.values()), fill_value=0)
    ct.to_csv(OUT / "crop_x_ales_class_0_100cm.csv")
    print("  crop_x_ales_class_0_100cm.csv:")
    print(ct)

    # ----- 3. Crop x policy priority -----------------------------------------
    df_p = df.dropna(subset=["priority"]).copy()
    df_p["pri_label"] = df_p["priority"].astype(int).map(PRI_LABELS)
    cp = pd.crosstab(df_p["CROP"], df_p["pri_label"]) \
            .reindex(columns=list(PRI_LABELS.values()), fill_value=0)
    cp.to_csv(OUT / "crop_x_priority_0_100cm.csv")

    # ----- 4. Per-crop KPIs --------------------------------------------------
    rows = []
    for crop, sub in df.groupby("CROP"):
        n_c   = len(sub)
        m_q50 = float(np.nanmean(sub["ales_idx"]))
        cls   = sub["ales_cls"].dropna().astype(int)
        pct_un = float(100 * cls.isin([4, 5]).mean()) if len(cls) else float("nan")
        pri   = sub["priority"].dropna().astype(int)
        pct_hi = float(100 * (pri == 1).mean()) if len(pri) else float("nan")
        m_dKs = float(np.nanmean(sub["cf_Ks_d"]))
        rows.append({
            "crop": crop, "n": n_c,
            "mean_ALES_q50":   round(m_q50, 2),
            "pct_N1_or_N2_%":  round(pct_un, 1),
            "pct_HIGH_pri_%":  round(pct_hi, 1),
            "mean_delta_Ks":   round(m_dKs, 2),
        })
    kpis = pd.DataFrame(rows).sort_values("n", ascending=False)
    kpis.to_csv(OUT / "crop_kpis.csv", index=False)
    print("  crop_kpis.csv:")
    print(kpis.to_string(index=False))

    # ----- 5. Summary --------------------------------------------------------
    with open(OUT / "stage7b_summary.txt", "w", encoding="utf-8") as f:
        f.write("Stage 7b - Independent crop-survey validation\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"Survey: Crop_March2025_11.shp  (n = {len(gdf)}, CRS = EPSG:32636)\n")
        f.write(f"Crop classes: {sorted(gdf['CROP'].unique().tolist())}\n\n")
        f.write(f"Sentinel-2 cropland mask vs survey:\n")
        f.write(f"  Overall accuracy = {oa:.3f}\n")
        f.write(f"  Cohen's kappa    = {kappa:.3f}\n")
        f.write(f"  Confusion (predicted LU at survey points):\n")
        for lab, row in conf.iterrows():
            f.write(f"    {lab:<14s}  n = {int(row['n_points']):4d}  "
                    f"({row['share_%']:5.2f} %)\n")
        f.write("\nCrop-resolved suitability KPIs (top by n):\n")
        for _, r in kpis.iterrows():
            f.write(f"  {r['crop']:<16s}  n={int(r['n']):4d}  "
                    f"q50={r['mean_ALES_q50']:5.2f}  "
                    f"%N1+N2={r['pct_N1_or_N2_%']:5.1f}  "
                    f"%HIGH={r['pct_HIGH_pri_%']:5.1f}  "
                    f"meanDKs={r['mean_delta_Ks']:5.2f}\n")
    print(f"\nAll Stage-7b outputs in {OUT}")


if __name__ == "__main__":
    main()
