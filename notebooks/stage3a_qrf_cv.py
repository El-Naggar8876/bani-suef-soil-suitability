"""
Stage 3a (v2) - QRF + Regression Kriging vs IDW/OK, with trimmed predictors

Improvements vs v1:
  1. Drop climate covariates with negligible AOI variance (CV < 2%) -- they
     act as spurious X/Y proxies at n=60.
  2. Add Regression Kriging (RK): QRF predicts the trend, OK on residuals.
  3. Per-property "best model" reported alongside QRF (kept for uncertainty).
  4. Tighter QRF: max_features='sqrt' to reduce predictor-set bias.

Outputs in outputs/stage3a/:
  cv_metrics.csv             - per (property, depth, model)
  cv_predictions.csv         - per-fold OOB predictions
  feature_importance.csv     - QRF Gini importance per (property, depth)
  best_model_per_property.csv
  predictors_used.txt
  stage3a_summary.txt
"""
from __future__ import annotations
from pathlib import Path
import warnings
import numpy as np
import pandas as pd
import geopandas as gpd
from sklearn.cluster import KMeans
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score
from scipy.spatial import cKDTree
from quantile_forest import RandomForestQuantileRegressor
from pykrige.ok import OrdinaryKriging

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

ROOT   = Path(__file__).resolve().parents[1]
STAGE1 = ROOT / "outputs" / "stage1"
STAGE2 = ROOT / "outputs" / "stage2"
OUT    = ROOT / "outputs" / "stage3a"
OUT.mkdir(parents=True, exist_ok=True)

RANDOM_STATE = 42
N_FOLDS = 5
N_TREES = 500
QUANTILES = [0.05, 0.50, 0.95]

TARGETS = ["EC", "ESP", "CaCO3", "pH", "CEC", "Clay", "SOM",
           "Ks", "AW", "N_avail", "P_avail", "K_avail", "Sand_total"]

DEPTHS = {
    "0_100cm": STAGE1 / "soil_profiles_0_100cm.gpkg",
    "0_30cm":  STAGE1 / "soil_profiles_0_30cm.gpkg",
}

# ---------------------------------------------------------------------------
def ccc(yt, yp):
    yt = np.asarray(yt, float); yp = np.asarray(yp, float)
    mt, mp = yt.mean(), yp.mean()
    cov = ((yt - mt) * (yp - mp)).mean()
    denom = yt.var() + yp.var() + (mt - mp) ** 2
    return float(2 * cov / denom) if denom > 0 else np.nan

def reg_metrics(yt, yp):
    return (float(np.sqrt(mean_squared_error(yt, yp))),
            float(mean_absolute_error(yt, yp)),
            float(r2_score(yt, yp)),
            ccc(yt, yp))

def picp_mpiw(yt, lo, hi):
    yt = np.asarray(yt, float)
    return float(((yt >= lo) & (yt <= hi)).mean()), float(np.mean(hi - lo))

def spatial_kfold(coords, k=N_FOLDS, seed=RANDOM_STATE):
    return KMeans(n_clusters=k, n_init=20, random_state=seed).fit(coords).labels_

def idw_predict(xy_train, y_train, xy_test, power=2, k_neighbours=8):
    tree = cKDTree(xy_train)
    dist, idx = tree.query(xy_test, k=min(k_neighbours, len(xy_train)))
    dist = np.where(dist == 0, 1e-9, dist)
    w = 1.0 / dist ** power
    return (w * y_train[idx]).sum(axis=1) / w.sum(axis=1)

def ok_predict(xy_train, y_train, xy_test, model="spherical"):
    try:
        ok = OrdinaryKriging(
            xy_train[:, 0], xy_train[:, 1], y_train,
            variogram_model=model, verbose=False, enable_plotting=False, nlags=8,
        )
        z, _ = ok.execute("points", xy_test[:, 0], xy_test[:, 1])
        return np.asarray(z, float)
    except Exception:
        return np.full(len(xy_test), float(np.mean(y_train)))

# ---------------------------------------------------------------------------
def load_depth(depth_key):
    soil = gpd.read_file(DEPTHS[depth_key], engine="pyogrio")
    soil["Sand_total"] = soil["Coarse_sand"] + soil["Fine_sand"]
    soil_utm = soil.to_crs(32636)
    soil["X_utm"] = soil_utm.geometry.x.values
    soil["Y_utm"] = soil_utm.geometry.y.values
    cov = pd.read_csv(STAGE2 / "covariates_at_profiles.csv")
    return soil.drop(columns="geometry").merge(cov, on="Profile", how="inner")

def trim_predictors(df, all_bands, cv_threshold=0.02):
    """Drop covariates whose CV across the 60 profiles < threshold."""
    keep, drop = [], []
    for b in all_bands:
        x = df[b].astype(float).to_numpy()
        if np.allclose(x.std(), 0):
            drop.append((b, 0.0)); continue
        cv = x.std() / (abs(x.mean()) + 1e-9)
        if cv < cv_threshold:
            drop.append((b, cv))
        else:
            keep.append(b)
    return keep, drop

# ---------------------------------------------------------------------------
def run_target(df, target, depth_key, predictor_cols):
    sub = df.dropna(subset=[target] + list(predictor_cols)).reset_index(drop=True)
    if len(sub) < 30:
        return [], [], []

    coords = sub[["X_utm", "Y_utm"]].to_numpy()
    fold_id = spatial_kfold(coords, k=N_FOLDS)
    X = sub[predictor_cols].to_numpy()
    y = sub[target].to_numpy()

    qrf_oof = np.full((len(sub), 3), np.nan)
    rk_oof  = np.full(len(sub), np.nan)
    idw_oof = np.full(len(sub), np.nan)
    ok_oof  = np.full(len(sub), np.nan)

    for f in range(N_FOLDS):
        train = fold_id != f; test = fold_id == f
        if test.sum() == 0 or train.sum() < 5:
            continue

        qrf = RandomForestQuantileRegressor(
            n_estimators=N_TREES, min_samples_leaf=3, max_features="sqrt",
            random_state=RANDOM_STATE, n_jobs=-1).fit(X[train], y[train])
        q = qrf.predict(X[test], quantiles=QUANTILES)
        qrf_oof[test] = q

        # Regression Kriging: QRF trend + OK on residuals
        trend_train = qrf.predict(X[train], quantiles=[0.50]).ravel()
        resid_train = y[train] - trend_train
        resid_test  = ok_predict(coords[train], resid_train, coords[test])
        rk_oof[test] = q[:, 1] + resid_test

        idw_oof[test] = idw_predict(coords[train], y[train], coords[test])
        ok_oof[test]  = ok_predict(coords[train], y[train], coords[test])

    valid = ~np.isnan(qrf_oof[:, 1])
    if valid.sum() < 10:
        return [], [], []
    yt = y[valid]

    rows_metric = []
    for name, pred in [("QRF", qrf_oof[valid, 1]),
                       ("RK",  rk_oof[valid]),
                       ("IDW", idw_oof[valid]),
                       ("OK",  ok_oof[valid])]:
        rmse, mae, r2, c = reg_metrics(yt, pred)
        if name == "QRF":
            picp, mpiw = picp_mpiw(yt, qrf_oof[valid, 0], qrf_oof[valid, 2])
        else:
            picp, mpiw = (np.nan, np.nan)
        rows_metric.append(dict(property=target, depth=depth_key, model=name,
                                n=int(valid.sum()), RMSE=rmse, MAE=mae,
                                R2=r2, CCC=c, PICP90=picp, MPIW=mpiw))

    rows_pred = []
    for i in np.where(valid)[0]:
        rows_pred.append(dict(
            Profile=int(sub.loc[i, "Profile"]),
            property=target, depth=depth_key,
            obs=float(y[i]),
            qrf_q05=float(qrf_oof[i, 0]),
            qrf_q50=float(qrf_oof[i, 1]),
            qrf_q95=float(qrf_oof[i, 2]),
            rk=float(rk_oof[i]), idw=float(idw_oof[i]), ok=float(ok_oof[i]),
            fold=int(fold_id[i]),
        ))

    qrf_full = RandomForestQuantileRegressor(
        n_estimators=N_TREES, min_samples_leaf=3, max_features="sqrt",
        random_state=RANDOM_STATE, n_jobs=-1).fit(X, y)
    rows_imp = [dict(property=target, depth=depth_key, predictor=c,
                     importance=float(v))
                for c, v in zip(predictor_cols, qrf_full.feature_importances_)]
    return rows_metric, rows_pred, rows_imp

# ---------------------------------------------------------------------------
def main():
    band_names = (STAGE2 / "stack_band_names.txt").read_text().splitlines()
    print(f"Stage-2 covariates available: {len(band_names)}")

    df_for_trim = load_depth("0_100cm")
    keep_bands, dropped = trim_predictors(df_for_trim, band_names, cv_threshold=0.02)
    print(f"Predictors kept after CV filter: {len(keep_bands)}")
    print("Dropped (low variance):")
    for b, cv in sorted(dropped, key=lambda t: t[1]):
        print(f"  {b:18s} CV={cv:.4f}")
    predictor_cols = keep_bands + ["X_utm", "Y_utm"]
    (OUT / "predictors_used.txt").write_text("\n".join(predictor_cols))

    all_metrics, all_preds, all_imp = [], [], []
    for depth_key in DEPTHS:
        print(f"\n=== Depth {depth_key} ===")
        df = load_depth(depth_key)
        for t in TARGETS:
            print(f"  -> {t}", flush=True)
            m, p, i = run_target(df, t, depth_key, predictor_cols)
            all_metrics.extend(m); all_preds.extend(p); all_imp.extend(i)

    metrics_df = pd.DataFrame(all_metrics)
    metrics_df.to_csv(OUT / "cv_metrics.csv", index=False)
    pd.DataFrame(all_preds).to_csv(OUT / "cv_predictions.csv", index=False)
    pd.DataFrame(all_imp).to_csv(OUT / "feature_importance.csv", index=False)

    pivot_rmse = metrics_df.pivot_table(index=["property", "depth"], columns="model",
                                         values="RMSE").round(3)
    pivot_r2   = metrics_df.pivot_table(index=["property", "depth"], columns="model",
                                         values="R2").round(3)

    best = (metrics_df.loc[metrics_df.groupby(["property", "depth"])["RMSE"].idxmin()]
            [["property", "depth", "model", "RMSE", "R2", "CCC"]])
    best.to_csv(OUT / "best_model_per_property.csv", index=False)
    win_counts = best["model"].value_counts().to_dict()

    qrf_picp = (metrics_df.query("model == 'QRF'")
                .groupby("depth")["PICP90"].mean().round(3))
    mean_metrics = (metrics_df.groupby(["model", "depth"])
                    [["RMSE", "R2", "CCC"]].mean().round(3))

    summary = []
    summary.append("Stage 3a v2 - QRF + Regression Kriging vs IDW/OK")
    summary.append("=" * 60)
    summary.append(f"Folds          : {N_FOLDS} (k-means clusters on UTM x,y)")
    summary.append(f"QRF            : {N_TREES} trees, max_features=sqrt, min_samples_leaf=3")
    summary.append(f"Predictors     : {len(predictor_cols)} (after CV>=2% filter + X,Y)")
    summary.append(f"Properties     : {TARGETS}")
    summary.append(f"Depths         : {list(DEPTHS)}")
    summary.append("")
    summary.append("RMSE per (property, depth):")
    summary.append(pivot_rmse.to_string()); summary.append("")
    summary.append("R2 per (property, depth):")
    summary.append(pivot_r2.to_string()); summary.append("")
    summary.append("Mean across all properties:")
    summary.append(mean_metrics.to_string()); summary.append("")
    summary.append(f"Best-RMSE wins per model: {win_counts}")
    summary.append("")
    summary.append(f"QRF mean 90% PI coverage (target ~0.90):\n{qrf_picp.to_string()}")

    text = "\n".join(summary)
    print("\n" + text)
    (OUT / "stage3a_summary.txt").write_text(text)
    print(f"\nWrote: {OUT}")

if __name__ == "__main__":
    main()
