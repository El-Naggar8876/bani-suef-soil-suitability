# Revision 1 stages (Geoderma Regional GEODRS-D-26-00837)

These scripts were added for the first revision. They depend on the outputs of
the original seven-stage pipeline in `../notebooks/`, which must have been run
first, and they write to `../outputs/stage_rev/`.

Run order:

```
python wapor_to_gtiff.py              # needs ../outputs/wapor/raw/*.npy
python stage10_sensitivity.py
python stage11_copula_mc.py 0_100cm independent 200 150000
python stage11_copula_mc.py 0_100cm correlated  200 150000
python stage11_copula_mc.py 0_30cm  independent 200 150000
python stage11_copula_mc.py 0_30cm  correlated  200 150000
python stage12_spatial_uncertainty.py
python stage13_wapor_check.py
python stage13b_wapor_robustness.py
python stage14_survey_bias.py
python stage15_cv_sampling.py
```

`stage10` imports the rating functions directly from
`../notebooks/stage5_ales_montecarlo.py`, so the sensitivity analysis is
guaranteed to use the same implementation that produced the published surfaces.

## WaPOR inputs

`outputs/wapor/` holds FAO WaPOR v3 level-2 (100 m) annual fields clipped to the
study area and reprojected to the 30 m analysis grid in EPSG:32636:
transpiration for 2025 and as a 2021-2025 mean, and evapotranspiration and total
biomass production as 2021-2025 means. The source rasters are global
cloud-optimised GeoTIFFs published by FAO at
`gs://fao-gismgr-wapor-3-data/DATA/WAPOR-3/MAPSET/`; only the study-area window
is read. Scale factors are 0.1 for T and AETI (mm/year) and 1 for TBP (kg/ha).

## Additional dependencies

Beyond `environment.yml`: `scipy` (copula, bootstrap, rank statistics) and
`python-docx` (used only by the two verification scripts, which read the
manuscript files).
