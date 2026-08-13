"""
check_setup.py - Pre-flight check for the Bani Suef irrigation-suitability pipeline.

Run this BEFORE running any analysis stage:

    python notebooks/check_setup.py

It reports, in plain language, whether your Python environment and your input
data are complete. It changes nothing and writes nothing. Uses only the Python
standard library, so it works even if the conda environment failed to build.

Exit code 0 = ready to run stages. Exit code 1 = something is missing.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

OK = "  OK      "
BAD = "  MISSING "
WARN = "  ---     "

problems: list[str] = []


def line(status: str, label: str, note: str = "") -> None:
    print(f"{status}{label}" + (f"   ({note})" if note else ""))


def header(title: str) -> None:
    print()
    print(title)
    print("-" * len(title))


def check_file(relpath: str, bundle: str | None = None) -> bool:
    path = ROOT / relpath
    if path.is_file():
        size = path.stat().st_size
        line(OK, relpath, human(size))
        return True
    line(BAD, relpath)
    if bundle:
        problems.append(f"{relpath}  -> {bundle}")
    return False


def check_shapefile(relpath: str, bundle: str | None = None) -> bool:
    """A shapefile is only usable if its .shx, .dbf and .prj siblings exist too."""
    base = ROOT / relpath
    if not base.is_file():
        line(BAD, relpath)
        if bundle:
            problems.append(f"{relpath}  -> {bundle}")
        return False
    missing = [ext for ext in (".shx", ".dbf", ".prj") if not base.with_suffix(ext).is_file()]
    if missing:
        line(BAD, relpath, "found, but companions missing: " + ", ".join(missing))
        problems.append(
            f"{relpath} is incomplete (missing {', '.join(missing)}) "
            "-> ask for the whole folder as a zip"
        )
        return False
    line(OK, relpath, "with .shx .dbf .prj")
    return True


def human(n: int) -> str:
    for unit in ("B", "KB", "MB", "GB"):
        if n < 1024 or unit == "GB":
            return f"{n:.0f} {unit}" if unit == "B" else f"{n:.1f} {unit}"
        n /= 1024.0
    return f"{n:.1f} GB"


# ---------------------------------------------------------------------------
print("=" * 68)
print("Bani Suef irrigation-suitability - setup check")
print("=" * 68)
print(f"Project folder: {ROOT}")

# 1. Python -----------------------------------------------------------------
header("1. Python")
v = sys.version_info
if (v.major, v.minor) >= (3, 10):
    line(OK, f"Python {v.major}.{v.minor}.{v.micro}")
else:
    line(BAD, f"Python {v.major}.{v.minor}.{v.micro}", "need 3.10 or newer; 3.12 recommended")
    problems.append("Python is too old - rebuild the conda environment")

# 2. Packages ---------------------------------------------------------------
header("2. Required Python packages")
REQUIRED = [
    ("numpy", "numerical arrays"),
    ("pandas", "tables"),
    ("scipy", "statistics"),
    ("sklearn", "machine learning (scikit-learn)"),
    ("matplotlib", "plots"),
    ("seaborn", "plots"),
    ("geopandas", "vector GIS"),
    ("shapely", "geometry"),
    ("rasterio", "raster GIS"),
    ("pyogrio", "fast vector reader"),
    ("pyproj", "coordinate systems"),
    ("quantile_forest", "quantile random forest (stages 3a, 3b)"),
    ("pykrige", "ordinary kriging (stage 4)"),
    ("openpyxl", "reads .xlsx"),
    ("xlrd", "reads .xls  (stage 1 needs this)"),
    ("tqdm", "progress bars"),
    ("joblib", "parallelism"),
]
OPTIONAL = [
    ("ee", "Earth Engine - only if you rebuild stage 2 yourself"),
    ("geemap", "Earth Engine helper - only for stage 2"),
    ("rioxarray", "used by some optional paths"),
    ("PIL", "image handling (pillow), stage 9"),
]

missing_pkgs = []
for mod, why in REQUIRED:
    if importlib.util.find_spec(mod) is not None:
        line(OK, mod, why)
    else:
        line(BAD, mod, why)
        missing_pkgs.append(mod)

for mod, why in OPTIONAL:
    line(OK if importlib.util.find_spec(mod) is not None else WARN, mod, why)

if missing_pkgs:
    problems.append(
        "Packages not installed: " + ", ".join(missing_pkgs) + "\n"
        "    Most likely cause: the environment is not active.\n"
        "    Fix:  conda activate bani-suef-suitability"
    )

# 3. Bundle A ---------------------------------------------------------------
header("3. Bundle A - field data (from the corresponding author)")
A = "Bundle A"
check_file("Analysis/Analysis_Banisuef_New_FF__March_2025.xls", A)
check_file("Analysis/water_analyses Benisueif_FFF.xls", A)
check_shapefile("Layers/Study_Area_Last.shp", A)
check_shapefile("Layers/Soil_Profiles.shp", A)
check_shapefile("Layers/Water_Samples.shp", A)
check_shapefile("Crop_March2025_11/Crop_March2025_11.shp", A)

# 4. Bundle B ---------------------------------------------------------------
header("4. Bundle B - pre-computed satellite covariates (skips stage 2)")
B = "Bundle B"
have_b = all([
    check_file("outputs/stage2/covariates_at_profiles.csv", B),
    check_file("outputs/stage2/covariates_at_water_samples.csv", B),
    check_file("outputs/stage2b_local_stack/covariate_stack_30m.tif", B),
])
if not have_b:
    print()
    print("  Bundle B is absent. The alternative is to rebuild the stack with")
    print("  Earth Engine (stages 2 and 2b). Checking whether that is set up:")
    key = ROOT / ".secrets" / "gee_service_account.json"
    line(OK if key.is_file() else WARN, ".secrets/gee_service_account.json",
         "service-account key present" if key.is_file() else "no service-account key")

    cfg = ROOT / "config.py"
    if not cfg.is_file():
        line(BAD, "config.py", "missing from the repository")
    else:
        text = cfg.read_text(encoding="utf-8", errors="replace")
        configured = any(
            ln.strip().startswith("GEE_PROJECT") and "None" not in ln.split("#")[0]
            for ln in text.splitlines()
        )
        line(OK if configured else WARN, "config.py",
             "GEE_PROJECT is set" if configured else "GEE_PROJECT not set yet")
        if not configured and not key.is_file():
            print()
            print("  To rebuild it yourself: set GEE_PROJECT in config.py, run")
            print("  'earthengine authenticate', then test with 'python config.py'.")
            print("  See SETUP.md Appendix A. Otherwise, request Bundle B.")

# 5. Repository files -------------------------------------------------------
header("5. Repository files (should already be present from GitHub)")
check_file("config.py")
check_file("data/countries.geojson")
check_file("data/ne_10m_rivers.geojson")
for s in ["stage1_data_audit", "stage3a_qrf_cv", "stage3b_predict_maps",
          "stage4_iwqi_surface", "stage5_ales_montecarlo", "stage6_counterfactual",
          "stage7_cropping_mismatch", "stage7b_crop_validation", "stage8_figures",
          "stage9_graphical_abstract"]:
    check_file(f"notebooks/{s}.py")

# 6. Progress ---------------------------------------------------------------
header("6. Stages already completed on this machine")
STAGE_MARKERS = [
    ("stage 1", "outputs/stage1/soil_profiles_0_100cm.gpkg"),
    ("stage 3a", "outputs/stage3a"),
    ("stage 3b", "outputs/stage3b/Ks_0_100cm_q50_qrf.tif"),
    ("stage 4", "outputs/stage4/iwqi_q50.tif"),
    ("stage 5", "outputs/stage5/ales_index_q50_0_100cm.tif"),
    ("stage 6", "outputs/stage6/gap_decomposition_0_100cm.csv"),
    ("stage 7", "outputs/stage7/policy_priority_0_100cm.tif"),
    ("stage 7b", "outputs/stage7b/crop_kpis.csv"),
    ("stage 8", "outputs/figures"),
]
for name, marker in STAGE_MARKERS:
    p = ROOT / marker
    done = p.exists() and (any(p.iterdir()) if p.is_dir() else True)
    line(OK if done else WARN, name, "done" if done else "not run yet")

# Verdict -------------------------------------------------------------------
print()
print("=" * 68)
if problems:
    print("NOT READY -", len(problems), "thing(s) need attention:")
    print("=" * 68)
    for i, p in enumerate(problems, 1):
        print(f"\n  {i}. {p}")
    print("\nSee SETUP.md, Step 4 and Step 5, for what to do about each.")
    sys.exit(1)

print("READY. You can start the pipeline:")
print("=" * 68)
print("""
    cd notebooks
    python stage1_data_audit.py
    python stage3a_qrf_cv.py
    python stage3b_predict_maps.py
    python stage4_iwqi_surface.py
    python stage5_ales_montecarlo.py
    python stage6_counterfactual.py
    python stage7_cropping_mismatch.py
    python stage7b_crop_validation.py
    python stage8_figures.py
    python stage9_graphical_abstract.py

Allow 2-4 hours in total. Stages 5 and 6 are the slow ones.
""")
sys.exit(0)
