# Bani Suef irrigation-suitability — probabilistic counterfactual digital soil mapping

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.PLACEHOLDER.svg)](https://doi.org/10.5281/zenodo.PLACEHOLDER)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.12](https://img.shields.io/badge/python-3.12-blue.svg)](https://www.python.org/downloads/release/python-3120/)

Reproducible analysis pipeline accompanying the manuscript:

> **ElNaggar, A. *et al.* (2026)** — *Drainage, not salinity, is the binding irrigation-suitability constraint in a Middle-Egypt Nile floodplain: a probabilistic, counterfactual assessment.* Submitted to *Geoderma Regional*.

The repository combines (i) **digital soil mapping** with quantile random forest (QRF) and rigorous spatial cross-validation, (ii) **probabilistic land evaluation** via Monte Carlo propagation through the Sys–Verheye ALES-Arid rating equations, (iii) **counterfactual gap decomposition** that isolates the binding limitation pixel-by-pixel, and (iv) **independent ground-truth validation** against a 993-point Sentinel-2 cropland survey (March 2025) over the 641.98 km² Beni Suef floodplain, Middle Egypt.

---

## Repository layout

```
.
├── notebooks/              # Sequential analysis scripts (stages 1–7 + auxiliary 8–9)
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

> **Note** Raw soil-profile measurements, irrigation-water samples and the unprocessed Sentinel-1/2 / SRTM / ERA5-Land tiles are **not** redistributed in this repository. See [`data/README.md`](data/README.md) for download instructions and the Zenodo dataset DOI used for archival snapshots.

---

## The seven analysis stages

| # | Stage | Script | Inputs | Outputs |
|---|-------|--------|--------|---------|
| 1 | Data audit & QC | `stage1_data_audit.py` | Raw soil profiles, water samples | `outputs/stage1/` summary tables |
| 2 | Covariate stack assembly | `stage2_covariate_stack.py` (+ `stage2b_download_stack.py`) | GEE: Sentinel-1 SAR, Sentinel-2 SR, ERA5-Land, SRTM | 53-band 30 m raster stack in `outputs/stage2/` |
| 3 | Spatial CV + QRF prediction | `stage3a_qrf_cv.py`, `stage3b_predict_maps.py` | Stage 1 + Stage 2 outputs | q05/q50/q95 soil-property rasters, PICP₉₀ = 0.85 |
| 4 | IWQI kriging surface | `stage4_iwqi_surface.py` | 20 surface-water samples | IWQI surface on 30 m grid |
| 5 | Probabilistic ALES-Arid | `stage5_ales_montecarlo.py` | Stage 3 quantile maps | 200 MC realisations of Sys–Verheye index, class probabilities |
| 6 | Counterfactual gap decomposition | `stage6_counterfactual.py` | Stage 5 baseline | 7 counterfactual scenarios, per-pixel Δ rasters |
| 7 | Cropping confrontation + validation | `stage7_cropping_mismatch.py`, `stage7b_crop_validation.py` | Stage 5 + Sentinel-2 cropland mask + 993 survey points | Policy-priority raster, contingency tables, user's accuracy = 0.988 |

Auxiliary stages 8 (publication figures) and 9 (graphical abstract) are *not* part of the seven core analytical stages but are included for full reproducibility of the manuscript's visuals.

---

## Quick start

### Option A — pip + venv

```powershell
git clone https://github.com/<USER>/bani-suef-soil-suitability.git
cd bani-suef-soil-suitability
python -m venv .venv
.venv\Scripts\Activate.ps1                      # PowerShell
# (or)  source .venv/bin/activate                # bash/zsh
pip install -r requirements.txt
```

### Option B — conda

```powershell
conda env create -f environment.yml
conda activate bani-suef-suitability
```

### Authenticating to Google Earth Engine (Stage 2 only)

Stage 2 fetches Sentinel-1/2, SRTM and ERA5-Land collections from Google Earth Engine. You must have an approved GEE account.

```powershell
earthengine authenticate
```

If you already have a downloaded covariate stack, you can skip Stage 2 entirely and place the local stack in `outputs/stage2b_local_stack/`.

### Running the full pipeline

```powershell
cd notebooks
python stage1_data_audit.py
python stage2_covariate_stack.py            # OR stage2b_download_stack.py for local
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
  author       = {ElNaggar, Ahmed and co-authors},
  title        = {Bani Suef irrigation-suitability digital-soil-mapping repository},
  year         = 2026,
  publisher    = {Zenodo},
  version      = {v1.0.0},
  doi          = {10.5281/zenodo.PLACEHOLDER},
  url          = {https://doi.org/10.5281/zenodo.PLACEHOLDER}
}

@article{elnaggar_2026_banisuef_article,
  author  = {ElNaggar, Ahmed and co-authors},
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

Corresponding author: **Dr Ahmed ElNaggar** — `a.elnaggar@un-ihe.org`

Issue tracker: please open a GitHub issue for bug reports or reproducibility questions.
