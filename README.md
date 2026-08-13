# Bani Suef irrigation-suitability — probabilistic counterfactual digital soil mapping

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21279524.svg)](https://doi.org/10.5281/zenodo.21279524)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.12](https://img.shields.io/badge/python-3.12-blue.svg)](https://www.python.org/downloads/release/python-3120/)

Reproducible analysis pipeline accompanying the manuscript:

> **El-Naggar, A. *et al.* (2026)** — *Drainage, not salinity, is the binding irrigation-suitability constraint in a Middle-Egypt Nile floodplain: a probabilistic, counterfactual assessment.* Submitted to *Geoderma Regional*.

> ### 👉 New here? Read [`SETUP.md`](SETUP.md) first.
>
> **This repository contains the analysis code, not the data.** Cloning it and
> running the scripts will fail immediately, by design: the field measurements
> are embargoed until manuscript acceptance, and the 30 m covariate stack is far
> too large for GitHub. [`SETUP.md`](SETUP.md) is a step-by-step, no-experience-
> assumed guide covering what to install, what data to request, and where to put
> it. Run `python notebooks/check_setup.py` at any point to see exactly what is
> still missing.

The repository combines (i) **digital soil mapping** with quantile random forest (QRF) and rigorous spatial cross-validation, (ii) **probabilistic land evaluation** via Monte Carlo propagation through the Sys–Verheye ALES-Arid rating equations, (iii) **counterfactual gap decomposition** that isolates the binding limitation pixel-by-pixel, and (iv) **independent ground-truth validation** against a 993-point Sentinel-2 cropland survey (March 2025) over the 641.98 km² Beni Suef floodplain, Middle Egypt.

---

## Repository layout

```
.
├── SETUP.md                # ← start here: full installation & run guide
├── notebooks/              # Sequential analysis scripts (stages 1–7 + auxiliary 8–9)
│   ├── check_setup.py                         # pre-flight check: what am I missing?
│   ├── stage1_data_audit.{py,ipynb}
│   ├── stage2_covariate_stack.{py,ipynb}     # main covariate stack (53 bands, 30 m)
│   ├── stage2b_download_stack.py              # local Earth Engine download helper
│   ├── stage3a_qrf_cv.py                      # spatial 5-fold CV benchmarking (QRF/RK/IDW/OK)
│   ├── stage3b_predict_maps.py                # QRF prediction at q05/q50/q95
│   ├── stage4_iwqi_surface.py                 # ordinary kriging of IWQI on 30 m grid
│   ├── stage5_ales_montecarlo.py              # 200 MC realisations of Sys–Verheye ALES-Arid
│   ├── stage6_counterfactual.py               # 7 do-style counterfactual scenarios
│   ├── stage7_cropping_mismatch.py            # ALES × Sentinel-2 cropland confrontation
│   ├── stage7b_crop_validation.py             # validation against 993 ground-truth points
│   ├── stage8_figures.py                      # Figures 1–6 (publication-ready)
│   └── stage9_graphical_abstract.py           # Graphical-abstract composition
│   └── README.md
├── data/                   # Lightweight reference vectors only (see data/README.md)
│   ├── countries.geojson
│   └── ne_10m_rivers.geojson
├── outputs/                # Recreated by running the pipeline (see outputs/README.md)
├── requirements.txt        # pip dependencies
├── environment.yml         # conda environment
├── CITATION.cff            # machine-readable citation (used by GitHub & Zenodo)
├── .zenodo.json            # Zenodo-side metadata
├── LICENSE                 # MIT
└── README.md               # this file
```

### Inputs that are *not* in this repository

The scripts additionally expect three folders at the repository root. They are
listed in `.gitignore` and are **not** published here, because the field data is
embargoed until manuscript acceptance:

```
Analysis/
├── Analysis_Banisuef_New_FF__March_2025.xls   # 60 soil profiles (sheet "total")
└── water_analyses Benisueif_FFF.xls           # 20 water samples (sheets "Total", "IWQ_FF")
Layers/
├── Study_Area_Last.shp                        # AOI polygon      (+ .shx .dbf .prj)
├── Soil_Profiles.shp                          # profile points   (+ .shx .dbf .prj)
└── Water_Samples.shp                          # sample points    (+ .shx .dbf .prj)
Crop_March2025_11/
└── Crop_March2025_11.shp                      # 993 survey points (+ .shx .dbf .prj)
```

All field datasets are in **EPSG:32636** (UTM 36N). Until the companion Zenodo
dataset is minted at acceptance, request these from the corresponding author.
See [`data/README.md`](data/README.md) and [`SETUP.md`](SETUP.md).

---

## The seven analysis stages

| # | Stage | Script | Inputs | Outputs |
|---|-------|--------|--------|---------|
| 1 | Data audit & QC | `stage1_data_audit.py` | Raw soil profiles, water samples | `outputs/stage1/` summary tables |
| 2 | Covariate stack assembly | `stage2_covariate_stack.py` **and** `stage2b_download_stack.py` | GEE: Sentinel-1 SAR, Sentinel-2 SR, ERA5-Land, SRTM | Covariates sampled at profiles/water samples in `outputs/stage2/`; 53-band 30 m raster stack in `outputs/stage2b_local_stack/` |
| 3 | Spatial CV + QRF prediction | `stage3a_qrf_cv.py`, `stage3b_predict_maps.py` | Stage 1 + Stage 2 outputs | q05/q50/q95 soil-property rasters, PICP₉₀ = 0.85 |
| 4 | IWQI kriging surface | `stage4_iwqi_surface.py` | 20 surface-water samples | IWQI surface on 30 m grid |
| 5 | Probabilistic ALES-Arid | `stage5_ales_montecarlo.py` | Stage 3 quantile maps | 200 MC realisations of Sys–Verheye index, class probabilities |
| 6 | Counterfactual gap decomposition | `stage6_counterfactual.py` | Stage 5 baseline | 7 counterfactual scenarios, per-pixel Δ rasters |
| 7 | Cropping confrontation + validation | `stage7_cropping_mismatch.py`, `stage7b_crop_validation.py` | Stage 5 + Sentinel-2 cropland mask + 993 survey points | Policy-priority raster, contingency tables, user's accuracy = 0.988 |

Auxiliary stages 8 (publication figures) and 9 (graphical abstract) are *not* part of the seven core analytical stages but are included for full reproducibility of the manuscript's visuals.

---

## Quick start

Full narrative instructions, including what to do if you have never used Python:
[`SETUP.md`](SETUP.md). The condensed version follows.

### Option A — conda (recommended; reliable on Windows)

```powershell
git clone https://github.com/El-Naggar8876/bani-suef-soil-suitability.git
cd bani-suef-soil-suitability
conda env create -f environment.yml
conda activate bani-suef-suitability
```

### Option B — pip + venv

```powershell
python -m venv .venv
.venv\Scripts\Activate.ps1                      # PowerShell
# (or)  source .venv/bin/activate                # bash/zsh
pip install -r requirements.txt
```

`geopandas`, `rasterio` and `fiona` need compiled GDAL bindings and frequently
fail to install this way on Windows. Prefer Option A.

### Verify your setup

```powershell
python notebooks/check_setup.py
```

This reports every required package and input file as `OK` or `MISSING`, and
exits non-zero if anything is absent. Run it before starting the pipeline.

### Google Earth Engine credentials (Stages 2 and 2b only)

Stages 2 and 2b fetch Sentinel-1/2, SRTM, CHIRPS, TerraClimate and ERA5-Land
from Google Earth Engine. They authenticate with a **Google Cloud service
account**, not with the interactive `earthengine authenticate` browser flow. The
key must be saved as:

```
.secrets/gee_service_account.json
```

`.secrets/` is git-ignored. **Never commit or share this file.**

**Most users should skip Stages 2 and 2b entirely** by obtaining the three
pre-computed artefacts from the corresponding author and placing them at:

```
outputs/stage2/covariates_at_profiles.csv
outputs/stage2/covariates_at_water_samples.csv
outputs/stage2b_local_stack/covariate_stack_30m.tif
```

Stages 3–9 then run with no Google account of any kind.

### Running the full pipeline

```powershell
cd notebooks
python stage1_data_audit.py
# Stages 2 and 2b — only if rebuilding the covariate stack yourself (see above).
# They are complementary, not alternatives: stage2 writes the covariate values
# sampled at the 60 profiles (needed by 3a and 3b); stage2b writes the
# wall-to-wall raster stack (needed by 3b, 4, 5, 6, 7 and 8). Run both.
#   python stage2_covariate_stack.py
#   python stage2b_download_stack.py
python stage3a_qrf_cv.py
python stage3b_predict_maps.py
python stage4_iwqi_surface.py
python stage5_ales_montecarlo.py
python stage6_counterfactual.py
python stage7_cropping_mismatch.py
python stage7b_crop_validation.py
python stage8_figures.py                    # optional: regenerate figures
python stage9_graphical_abstract.py         # optional: regenerate graphical abstract
```

Total runtime on a workstation with 16 GB RAM and 8 cores: approximately 2–4 hours, dominated by Stage 5 (Monte Carlo) and Stage 6 (counterfactuals). Stage 2 download time depends on Earth Engine queue.

Each stage reads its inputs from disk and writes its outputs to disk, so the run
can be interrupted after any stage and resumed later.

---

## Reproducing the manuscript results

The published numbers in the manuscript correspond to the seed-fixed run with `RANDOM_STATE = 42` declared at the top of every stage. Key reproducibility hooks:

- Spatial 5-fold k-Means cross-validation seeded for stage 3a.
- 200 Monte Carlo realisations (stage 5) and 100 realisations per counterfactual scenario (stage 6) seeded.
- Sentinel-2 cropland mask date-locked to the **March 2025** composite to align with the 993-point ground-truth survey.

---

## Data availability

| Dataset | Location | License |
|---------|----------|---------|
| 60 soil-profile measurements | Zenodo dataset DOI (companion archive) | CC-BY 4.0 |
| 20 surface-water IWQI samples | Zenodo dataset DOI (companion archive) | CC-BY 4.0 |
| 993-point crop-survey (March 2025) | Zenodo dataset DOI (companion archive) | CC-BY 4.0 |
| Sentinel-1/2, ERA5-Land, SRTM | Google Earth Engine (Copernicus / NASA) | CC-BY-SA / public domain |
| Country & rivers reference vectors | Natural Earth (bundled in `data/`) | Public domain |

Zenodo dataset DOI will be minted at acceptance and inserted here.

---

## Citation

If you use this code, please cite **both** the software and the article:

```bibtex
@software{elnaggar_2026_banisuef_code,
  author       = {El-Naggar, Ahmed;  Abou Alfotoh, M.S.M; Abdel Ghaffar, M. K; El-GamalB. A; Abdellatif D. Abdellatif},
  title        = {Bani Suef irrigation-suitability digital-soil-mapping repository},
  year         = 2026,
  publisher    = {Zenodo},
  version      = {v1.0.0},
  doi          = {10.5281/zenodo.21279524},
  url          = {https://doi.org/10.5281/zenodo.21279524}
}

@article{elnaggar_2026_banisuef_article,
  author  = {El-Naggar, Ahmed;  Abou Alfotoh, M.S.M; Abdel Ghaffar, M. K; El-GamalB. A; Abdellatif D. Abdellatif},
  title   = {Drainage, not salinity, is the binding irrigation-suitability constraint in a Middle-Egypt Nile floodplain: a probabilistic, counterfactual assessment},
  journal = {Geoderma Regional},
  year    = 2026,
  note    = {Submitted}
}
```

A machine-readable citation is provided in [`CITATION.cff`](CITATION.cff).

---

## Licence

Code: [MIT](LICENSE).
Data (when bundled): [CC-BY 4.0](https://creativecommons.org/licenses/by/4.0/).

---

## Acknowledgements

This work was supported by IHE Delft, the Soil, Water and Environment Research Institute (SWERI), and the Egyptian Agricultural Research Centre. Data acquisition relied on the Copernicus Sentinel programme, NASA's SRTM mission, and the ECMWF ERA5-Land reanalysis. We thank the Google Earth Engine team for cloud-resident access to the imagery archives.

---

## Contact

Corresponding author: **Dr Ahmed El-Naggar** — `a.elnaggar@un-ihe.org`

Issue tracker: please open a GitHub issue for bug reports or reproducibility questions.
