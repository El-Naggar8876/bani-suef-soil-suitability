"""
Stage 1 — Data audit and harmonization
Beni Suef soil/water suitability paper (Agric. Water Manag. submission)

Tasks:
  1. Load and clean soil profile horizon data (60 profiles, ~4 horizons each).
  2. Aggregate horizons to depth-weighted means over 0-100 cm (primary)
     and 0-30 cm (sensitivity), per profile per property.
  3. Audit: missing values, outliers (IQR), cation-anion ion balance error,
     EC vs sum-of-cations consistency.
  4. Join cleaned per-profile table to Soil_Profiles.shp -> GeoPackage.
  5. Load and audit water samples; verify Meireles-style IWQI computation
     against the published formulae (Meireles et al. 2010).
  6. Write outputs to outputs/stage1/.
"""

from pathlib import Path
import numpy as np
import pandas as pd
import geopandas as gpd

ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = ROOT / "Analysis"
LAYERS = ROOT / "Layers"
OUT = ROOT / "outputs" / "stage1"
OUT.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# 1. Load soil horizon table
# ---------------------------------------------------------------------------
soil_path = ANALYSIS / "Analysis_Banisuef_New_FF__March_2025.xls"
raw = pd.read_excel(soil_path, sheet_name="total", header=0)
# Row 0 holds units, drop it
raw = raw.iloc[1:].reset_index(drop=True)

# Forward-fill profile id (only first horizon row carries it)
raw["profile"] = raw["profile"].ffill()
raw["profile"] = pd.to_numeric(raw["profile"], errors="coerce").astype("Int64")

# Parse depth ranges "0 - 15", "15 - 30", "30 -  70", etc.
def parse_depth(s):
    if pd.isna(s):
        return (np.nan, np.nan)
    s = str(s).replace("–", "-").replace("—", "-")
    parts = [p.strip() for p in s.split("-")]
    parts = [p for p in parts if p != ""]
    try:
        top = float(parts[0])
        bot = float(parts[1])
        return (top, bot)
    except Exception:
        return (np.nan, np.nan)

raw[["depth_top", "depth_bot"]] = raw["depth"].apply(
    lambda s: pd.Series(parse_depth(s))
)
raw["thickness"] = raw["depth_bot"] - raw["depth_top"]

# Numeric columns we will aggregate
NUMERIC_COLS = [
    "PH", "EC", "CO3--", "HCO3-", "Cl-", "SO4-2",
    "Ca++", "Mg++", "Na+", "K+", "SAR", "ESP", "SP",
    "C. sand", "F. Sand", "Silt", "Clay",
    "SOM", "Caco3", "CEC", "AW", "BD", "Ks",
    "N", "P", "K", "Ca", "Mg", "Fe", "Zn", "Mn", "Cu",
]
for c in NUMERIC_COLS:
    raw[c] = pd.to_numeric(raw[c], errors="coerce")

# Texture passes through as the modal class for the profile
raw["Texture"] = raw["Texture"].astype("string")

print(f"Horizon rows after cleaning: {len(raw)}")
print(f"Distinct profiles: {raw['profile'].nunique()}")
print(f"Min/Max profile id: {raw['profile'].min()} / {raw['profile'].max()}")

# Profile depth (deepest bottom per profile)
profile_depth = raw.groupby("profile")["depth_bot"].max().rename("profile_depth_cm")

# ---------------------------------------------------------------------------
# 2. Depth-weighted aggregation
# ---------------------------------------------------------------------------
def depth_weighted(df: pd.DataFrame, cols, top_limit=0, bot_limit=100):
    """Depth-weighted mean of `cols` over [top_limit, bot_limit] cm per profile."""
    df = df.copy()
    # Clip horizons to the slice
    df["t"] = df["depth_top"].clip(lower=top_limit, upper=bot_limit)
    df["b"] = df["depth_bot"].clip(lower=top_limit, upper=bot_limit)
    df["w"] = (df["b"] - df["t"]).clip(lower=0)
    df = df[df["w"] > 0]
    out = {}
    for c in cols:
        sub = df.dropna(subset=[c])
        if sub.empty:
            out[c] = np.nan
        else:
            out[c] = np.average(sub[c], weights=sub["w"])
    out["coverage_cm"] = df["w"].sum()
    return pd.Series(out)

agg_100 = (
    raw.groupby("profile")
    .apply(lambda g: depth_weighted(g, NUMERIC_COLS, 0, 100), include_groups=False)
    .reset_index()
)
agg_30 = (
    raw.groupby("profile")
    .apply(lambda g: depth_weighted(g, NUMERIC_COLS, 0, 30), include_groups=False)
    .reset_index()
)

# Modal texture per profile (most frequent across horizons, weighted by thickness)
def modal_texture(g):
    s = g.dropna(subset=["Texture"])
    if s.empty:
        return pd.NA
    return (
        s.groupby("Texture")["thickness"].sum().sort_values(ascending=False).index[0]
    )

texture_mode = raw.groupby("profile").apply(modal_texture, include_groups=False)
texture_mode.name = "Texture_mode"

agg_100 = agg_100.merge(profile_depth, on="profile").merge(
    texture_mode, on="profile"
)
agg_30 = agg_30.merge(profile_depth, on="profile").merge(
    texture_mode, on="profile"
)

print(f"\nAggregated profiles 0-100 cm: {len(agg_100)}")
print(f"Aggregated profiles 0-30 cm:  {len(agg_30)}")

# ---------------------------------------------------------------------------
# 3. Audit
# ---------------------------------------------------------------------------
audit_lines = []

# 3a. Missing values
miss = agg_100[NUMERIC_COLS].isna().sum()
audit_lines.append("== Missing values per property (0-100 cm aggregated, n=60) ==")
audit_lines.append(miss.to_string())

# 3b. IQR-based outlier flags per property (in the horizon data)
def iqr_flags(s):
    q1, q3 = s.quantile([0.25, 0.75])
    iqr = q3 - q1
    return ((s < q1 - 3 * iqr) | (s > q3 + 3 * iqr)).sum()

out_counts = {c: iqr_flags(raw[c]) for c in NUMERIC_COLS}
audit_lines.append("\n== Extreme outliers per property (>3*IQR, horizon level) ==")
audit_lines.append(pd.Series(out_counts).to_string())

# 3c. Cation-anion balance (meq/L)
cat = raw[["Ca++", "Mg++", "Na+", "K+"]].sum(axis=1)
ani = raw[["CO3--", "HCO3-", "Cl-", "SO4-2"]].sum(axis=1)
bal = (cat - ani) / ((cat + ani) / 2.0) * 100.0  # % difference
raw["ion_balance_pct"] = bal
audit_lines.append("\n== Cation-anion balance error (%) summary ==")
audit_lines.append(bal.describe().to_string())
n_bad_balance = (bal.abs() > 10).sum()
audit_lines.append(f"\nHorizons with |balance error| > 10%: {n_bad_balance} of {len(bal)}")
if bal.abs().max() < 0.01:
    audit_lines.append(
        "  WARNING: cation-anion sums match to machine precision (< 1e-2 %). "
        "This means at least one ion (likely Na+) was back-calculated to force "
        "electroneutrality rather than independently measured. The ion balance "
        "is therefore *not* an independent QC check on this dataset; we will "
        "declare this transparently in the Methods."
    )

# 3d. EC vs TDS-from-cations sanity  (sum cations meq/L * ~64 ~= TDS mg/L; EC dS/m * 640 ~= TDS)
audit_lines.append("\n== EC range (dS/m) and sum-cations (meq/L) ==")
audit_lines.append(f"EC: min {raw['EC'].min():.2f}  max {raw['EC'].max():.2f}  mean {raw['EC'].mean():.2f}")
audit_lines.append(f"Sum cations meq/L: min {cat.min():.2f}  max {cat.max():.2f}  mean {cat.mean():.2f}")

# 3e. Texture vs % sand+silt+clay sums
psum = raw[["C. sand", "F. Sand", "Silt", "Clay"]].sum(axis=1)
audit_lines.append("\n== Particle-size sum (%) summary (should be ~100) ==")
audit_lines.append(psum.describe().to_string())

# Coverage check
audit_lines.append("\n== 0-100 cm coverage per profile (cm) summary ==")
audit_lines.append(agg_100["coverage_cm"].describe().to_string())
short_profiles = agg_100[agg_100["coverage_cm"] < 80]["profile"].tolist()
audit_lines.append(f"Profiles with <80 cm coverage in 0-100 cm slice: {len(short_profiles)}")
if short_profiles:
    audit_lines.append(f"  -> Profile IDs: {short_profiles}")

audit_text = "\n".join(audit_lines)
print("\n" + audit_text)
(OUT / "stage1_audit_report.txt").write_text(audit_text, encoding="utf-8")

# ---------------------------------------------------------------------------
# 4. Join to Soil_Profiles.shp
# ---------------------------------------------------------------------------
sp = gpd.read_file(LAYERS / "Soil_Profiles.shp", engine="pyogrio")
sp["Profile"] = sp["Profile"].astype(int)
agg_100["profile"] = agg_100["profile"].astype(int)
agg_30["profile"] = agg_30["profile"].astype(int)

# Normalize column names for SQLite/GPKG compatibility
RENAME = {
    "PH": "pH", "EC": "EC", "CO3--": "CO3", "HCO3-": "HCO3", "Cl-": "Cl",
    "SO4-2": "SO4", "Ca++": "Ca_meq", "Mg++": "Mg_meq", "Na+": "Na_meq",
    "K+": "K_meq", "C. sand": "Coarse_sand", "F. Sand": "Fine_sand",
    "Caco3": "CaCO3", "N": "N_avail", "P": "P_avail", "K": "K_avail",
    "Ca": "Ca_avail", "Mg": "Mg_avail",
}
agg_100 = agg_100.rename(columns=RENAME)
agg_30 = agg_30.rename(columns=RENAME)
agg_100 = agg_100.rename(columns={"profile": "Profile"})
agg_30 = agg_30.rename(columns={"profile": "Profile"})

sp_join_100 = sp.merge(agg_100, on="Profile", how="left")
sp_join_30 = sp.merge(agg_30, on="Profile", how="left")

n_unmatched = sp_join_100[list(RENAME.values())[0]].isna().sum()
print(f"\nProfiles in shapefile without matched analytics: {n_unmatched}")

sp_join_100.to_file(OUT / "soil_profiles_0_100cm.gpkg", layer="profiles", driver="GPKG", engine="pyogrio")
sp_join_30.to_file(OUT / "soil_profiles_0_30cm.gpkg", layer="profiles", driver="GPKG", engine="pyogrio")
agg_100.to_csv(OUT / "soil_per_profile_0_100cm.csv", index=False)
agg_30.to_csv(OUT / "soil_per_profile_0_30cm.csv", index=False)
raw.to_csv(OUT / "soil_horizons_clean.csv", index=False)

# ---------------------------------------------------------------------------
# 5. Water samples & Meireles IWQI verification
# ---------------------------------------------------------------------------
w_path = ANALYSIS / "water_analyses Benisueif_FFF.xls"
w_total = pd.read_excel(w_path, sheet_name="Total", header=0).iloc[1:].reset_index(drop=True)
w_total = w_total.dropna(subset=["NO"]).copy()
for c in ["pH", "EC ", "Na", "Ca", "Mg", "K", "HCO3", "Cl", "So4", "SAR"]:
    w_total[c] = pd.to_numeric(w_total[c], errors="coerce")
w_total = w_total.rename(columns={"EC ": "EC", "So4": "SO4"})
w_total["NO"] = w_total["NO"].astype(int)

print(f"\nWater samples loaded: {len(w_total)}")
print(w_total[["NO", "Name", "pH", "EC", "Na", "Cl", "HCO3", "SAR"]].head())

# --- Meireles et al. (2010) IWQI implementation ---------------------------
# Five parameters with normalized quality (qi) functions and weights wi.
# qi piecewise-linear functions defined by parameter ranges in the paper.
# Weights (wi): EC=0.211, SAR=0.189, Na=0.204, Cl=0.194, HCO3=0.202.
#   IWQI = sum_i (qi * wi);  range 0-100; classes:
#     85-100 No restriction; 70-85 Low; 55-70 Moderate; 40-55 High; 0-40 Severe.
# Reference: Meireles ACM, Andrade EM, Chaves LCG, Frischkorn H, Crisostomo LA (2010).
#   A new proposal of the classification of irrigation water. Rev. Cienc. Agron. 41(3): 349-357.

def qi_EC(ec):
    # ec in dS/m at 25 C
    if ec < 0.20:
        return _interp(ec, 0.0, 0.20, 0, 35)        # 0-35
    if ec < 0.75:
        return _interp(ec, 0.20, 0.75, 35, 100)     # 35-100 (best)
    if ec < 1.50:
        return _interp(ec, 0.75, 1.50, 100, 60)
    if ec < 3.00:
        return _interp(ec, 1.50, 3.00, 60, 0)
    return 0

def qi_SAR(sar):
    if sar < 2:
        return _interp(sar, 0, 2, 80, 100)
    if sar < 3:
        return _interp(sar, 2, 3, 100, 80)
    if sar < 6:
        return _interp(sar, 3, 6, 80, 60)
    if sar < 12:
        return _interp(sar, 6, 12, 60, 0)
    return 0

def qi_Na(na):  # meq/L
    if na < 2:
        return _interp(na, 0, 2, 80, 100)
    if na < 3:
        return _interp(na, 2, 3, 100, 80)
    if na < 6:
        return _interp(na, 3, 6, 80, 60)
    if na < 9:
        return _interp(na, 6, 9, 60, 0)
    return 0

def qi_Cl(cl):  # meq/L
    if cl < 1:
        return _interp(cl, 0, 1, 80, 100)
    if cl < 4:
        return _interp(cl, 1, 4, 100, 80)
    if cl < 7:
        return _interp(cl, 4, 7, 80, 60)
    if cl < 10:
        return _interp(cl, 7, 10, 60, 0)
    return 0

def qi_HCO3(hco3):  # meq/L
    if hco3 < 1:
        return _interp(hco3, 0, 1, 80, 100)
    if hco3 < 1.5:
        return _interp(hco3, 1, 1.5, 100, 80)
    if hco3 < 4.5:
        return _interp(hco3, 1.5, 4.5, 80, 60)
    if hco3 < 8.5:
        return _interp(hco3, 4.5, 8.5, 60, 0)
    return 0

def _interp(x, x0, x1, y0, y1):
    if x1 == x0:
        return y0
    return y0 + (y1 - y0) * (x - x0) / (x1 - x0)

WEIGHTS = {"EC": 0.211, "SAR": 0.189, "Na": 0.204, "Cl": 0.194, "HCO3": 0.202}

def meireles_iwqi(row):
    qi = {
        "EC": qi_EC(row["EC"]),
        "SAR": qi_SAR(row["SAR"]),
        "Na": qi_Na(row["Na"]),
        "Cl": qi_Cl(row["Cl"]),
        "HCO3": qi_HCO3(row["HCO3"]),
    }
    iwqi = sum(qi[k] * WEIGHTS[k] for k in WEIGHTS)
    return pd.Series({**{f"qi_{k}": v for k, v in qi.items()}, "IWQI_meireles": iwqi})

iwqi_df = w_total.apply(meireles_iwqi, axis=1)
w_total = pd.concat([w_total, iwqi_df], axis=1)

def iwqi_class(v):
    if v >= 85:  return "No restriction"
    if v >= 70:  return "Low restriction"
    if v >= 55:  return "Moderate restriction"
    if v >= 40:  return "High restriction"
    return "Severe restriction"

w_total["IWQI_class"] = w_total["IWQI_meireles"].apply(iwqi_class)

# Compare to the IWQI already in the spreadsheet
iwq_ff = pd.read_excel(w_path, sheet_name="IWQ_FF", header=0).iloc[1:].reset_index(drop=True)
# Last column holds the precomputed WQI
iwq_ff = iwq_ff.dropna(subset=["NO"])
iwq_ff["NO"] = pd.to_numeric(iwq_ff["NO"], errors="coerce").astype("Int64")
iwq_ff = iwq_ff.drop_duplicates(subset=["NO"], keep="first")
precomp = iwq_ff[["NO", "Unnamed: 18"]].rename(columns={"Unnamed: 18": "IWQI_precomputed"})
precomp["IWQI_precomputed"] = pd.to_numeric(precomp["IWQI_precomputed"], errors="coerce")

cmp = w_total[["NO", "Name", "IWQI_meireles", "IWQI_class"]].merge(precomp, on="NO", how="left")
cmp["diff"] = cmp["IWQI_meireles"] - cmp["IWQI_precomputed"]
print("\n== Meireles IWQI: my recomputation vs spreadsheet ==")
print(cmp.to_string(index=False))

# Join water samples to shapefile
ws = gpd.read_file(LAYERS / "Water_Samples.shp", engine="pyogrio")
ws["No"] = ws["No"].astype(int)
ws_join = ws.merge(
    w_total.drop(columns=["Name"]).rename(columns={"NO": "No"}),
    on="No", how="left",
)
ws_join.to_file(OUT / "water_samples_iwqi.gpkg", layer="water", driver="GPKG", engine="pyogrio")
w_total.to_csv(OUT / "water_samples_iwqi.csv", index=False)
cmp.to_csv(OUT / "iwqi_recompute_vs_spreadsheet.csv", index=False)

# ---------------------------------------------------------------------------
# Final summary
# ---------------------------------------------------------------------------
summary = {
    "n_horizons_raw": int(len(raw)),
    "n_profiles": int(raw["profile"].nunique()),
    "n_profiles_in_shp": int(len(sp)),
    "n_water_samples": int(len(w_total)),
    "iwqi_min": float(w_total["IWQI_meireles"].min()),
    "iwqi_max": float(w_total["IWQI_meireles"].max()),
    "iwqi_mean": float(w_total["IWQI_meireles"].mean()),
}
print("\n== Stage 1 summary ==")
for k, v in summary.items():
    print(f"  {k}: {v}")

print("\nIWQI class counts:")
print(w_total["IWQI_class"].value_counts())
print("\nOutputs written to:", OUT)
