"""
Stage 15 - Cross-validation reporting and sampling design description.

Purpose (Reviewer 2, comment 3.1): the submitted manuscript does not state the
sampling density, nor give full cross-validation metrics for the interpolator
that generated the published q50 surfaces.

Outputs -> 02_New_Analysis/outputs_rev/
"""
from __future__ import annotations
import pathlib
import numpy as np, pandas as pd, geopandas as gpd, rasterio
from scipy.spatial import cKDTree

import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parent))
from revpaths import REPO as _REPO, OUT as _OUT, WAPOR as _WAPOR

HERE = pathlib.Path(__file__).resolve().parent
REPO = _REPO
OUT = _OUT

ALES_PROPS = ["EC", "ESP", "CaCO3", "SOM", "CEC", "Clay", "Sand_total", "Ks"]


def main():
    cv = pd.read_csv(REPO / "outputs" / "stage3a" / "cv_metrics.csv")
    cv = cv[cv.property.isin(ALES_PROPS)].copy()
    for c in ["RMSE", "MAE", "R2", "CCC", "PICP90", "MPIW"]:
        cv[c] = pd.to_numeric(cv[c], errors="coerce")

    # full table, all learners, both depths
    full = cv.pivot_table(index=["depth", "property"], columns="model",
                          values=["R2", "RMSE", "MAE", "CCC"]).round(3)
    full.to_csv(OUT / "cv_full_metrics.csv")
    print("===== full cross-validation table (R2) =====")
    print(cv.pivot_table(index=["depth", "property"], columns="model", values="R2").round(3).to_string())
    print("\n===== RMSE =====")
    print(cv.pivot_table(index=["depth", "property"], columns="model", values="RMSE").round(3).to_string())

    # best learner per property-depth + its metrics, and the QRF interval metrics
    rows = []
    for (dep, prop), g in cv.groupby(["depth", "property"]):
        b = g.loc[g.R2.idxmax()]
        q = g[g.model == "QRF"].iloc[0]
        rows.append(dict(depth=dep, property=prop, best_learner=b.model,
                         best_R2=round(float(b.R2), 3), best_RMSE=round(float(b.RMSE), 3),
                         best_MAE=round(float(b.MAE), 3), best_CCC=round(float(b.CCC), 3),
                         QRF_R2=round(float(q.R2), 3),
                         QRF_PICP90=round(float(q.PICP90), 3) if np.isfinite(q.PICP90) else np.nan,
                         QRF_MPIW90=round(float(q.MPIW), 3) if np.isfinite(q.MPIW) else np.nan))
    bt = pd.DataFrame(rows).sort_values(["depth", "property"])
    bt.to_csv(OUT / "cv_best_learner_table.csv", index=False)
    print("\n===== best learner per property, with QRF interval calibration =====")
    print(bt.to_string(index=False))
    print("\nmean QRF PICP90 by depth:")
    print(bt.groupby("depth").QRF_PICP90.mean().round(3).to_string())

    # ---- sampling design -------------------------------------------------
    g = gpd.read_file(REPO / "outputs" / "stage1" / "soil_profiles_0_100cm.gpkg")
    xy = np.column_stack([g.geometry.x.values, g.geometry.y.values])
    tree = cKDTree(xy)
    d, _ = tree.query(xy, k=2)
    nn = d[:, 1]
    with rasterio.open(REPO / "outputs" / "stage5" / "ales_index_q50_0_100cm.tif") as s:
        valid = np.isfinite(s.read(1))
    area_km2 = float(valid.sum()) * 0.0009
    water = gpd.read_file(REPO / "outputs" / "stage1" / "water_samples_iwqi.gpkg")
    wxy = np.column_stack([water.geometry.x.values, water.geometry.y.values])
    wd, _ = cKDTree(wxy).query(wxy, k=2)

    summ = dict(
        mapped_area_km2=round(area_km2, 2),
        n_soil_profiles=len(g),
        profile_density_per_100km2=round(100 * len(g) / area_km2, 2),
        mean_area_per_profile_km2=round(area_km2 / len(g), 2),
        profile_nn_min_m=round(float(nn.min()), 0),
        profile_nn_median_m=round(float(np.median(nn)), 0),
        profile_nn_mean_m=round(float(nn.mean()), 0),
        profile_nn_max_m=round(float(nn.max()), 0),
        n_water_samples=len(water),
        water_nn_min_m=round(float(wd[:, 1].min()), 0),
        water_nn_median_m=round(float(np.median(wd[:, 1])), 0),
        water_nn_max_m=round(float(wd[:, 1].max()), 0),
        n_crop_survey_points=993,
        crop_point_density_per_km2_cropland=round(993 / (540958 * 0.0009), 2),
    )
    pd.DataFrame([summ]).T.rename(columns={0: "value"}).to_csv(OUT / "sampling_design.csv")
    print("\n===== sampling design =====")
    for k, v in summ.items():
        print(f"  {k:36s} {v}")

if __name__ == "__main__":
    main()
