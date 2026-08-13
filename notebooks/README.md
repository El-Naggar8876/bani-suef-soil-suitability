# `notebooks/` — sequential analysis scripts

The seven analysis stages described in §2.11 of the manuscript are implemented as standalone Python scripts (with companion notebooks where applicable) that are designed to be run in order. Each stage reads its inputs from `outputs/<previous-stage>/` (or `data/`) and writes its outputs to `outputs/<this-stage>/`.

Before running anything, verify your installation and input data with
[`check_setup.py`](check_setup.py):

```
python notebooks/check_setup.py
```

| Order | Script | Purpose |
|-------|--------|---------|
| 0 | [`check_setup.py`](check_setup.py) | Pre-flight check — reports missing packages and missing input files |
| 1 | [`stage1_data_audit.py`](stage1_data_audit.py) | Quality-control of soil and water datasets |
| 2 | [`stage2_covariate_stack.py`](stage2_covariate_stack.py) | Assemble 53-band 30 m covariate stack from Earth Engine |
| 2b | [`stage2b_download_stack.py`](stage2b_download_stack.py) | Downloads the wall-to-wall stack locally — **required alongside stage 2**, not an alternative to it (stage 2 writes the covariates sampled at points; stage 2b writes the raster used by stages 3b–8) |
| 3a | [`stage3a_qrf_cv.py`](stage3a_qrf_cv.py) | Spatial 5-fold cross-validation benchmarking (QRF / RK / IDW / OK) |
| 3b | [`stage3b_predict_maps.py`](stage3b_predict_maps.py) | QRF prediction at q05/q50/q95 for 8 soil properties |
| 4 | [`stage4_iwqi_surface.py`](stage4_iwqi_surface.py) | Ordinary kriging of the Irrigation Water Quality Index |
| 5 | [`stage5_ales_montecarlo.py`](stage5_ales_montecarlo.py) | 200 Monte Carlo realisations of Sys–Verheye ALES-Arid |
| 6 | [`stage6_counterfactual.py`](stage6_counterfactual.py) | Seven *do*-style counterfactual scenarios |
| 7 | [`stage7_cropping_mismatch.py`](stage7_cropping_mismatch.py) | ALES × Sentinel-2 cropland confrontation; policy-priority raster |
| 7b | [`stage7b_crop_validation.py`](stage7b_crop_validation.py) | Validation against 993 ground-truth crop-survey points |
| 8 (aux) | [`stage8_figures.py`](stage8_figures.py) | Publication figures (1–6) |
| 9 (aux) | [`stage9_graphical_abstract.py`](stage9_graphical_abstract.py) | Graphical abstract composition |

Stages 8 and 9 produce the manuscript figures only; they are *not* part of the seven analytical stages.

## Notebook companions

`stage1_data_audit.ipynb` and `stage2_covariate_stack.ipynb` mirror the corresponding `.py` scripts and offer interactive exploration of intermediate diagnostics. The `.py` scripts are the canonical execution path used in the published results.
