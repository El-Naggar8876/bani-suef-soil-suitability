# `data/` — input datasets

This directory holds **only the small, redistributable reference vectors** used by the figure-generation stages:

| File | Source | License |
|------|--------|---------|
| `countries.geojson` | Natural Earth (1:10m cultural) | Public domain |
| `ne_10m_rivers.geojson` | Natural Earth (1:10m physical) | Public domain |

## Larger inputs are NOT bundled here

To keep the repository lightweight and respect data-sharing constraints, the following inputs are **not** committed and must be obtained separately:

| Dataset | Used by | Where to get it |
|---------|---------|-----------------|
| 60 soil-profile measurements (chemistry + Ks) | Stages 1, 3, 5, 6 | Companion Zenodo dataset DOI (minted at acceptance) |
| 20 surface-water IWQI samples | Stage 4 | Companion Zenodo dataset DOI |
| 993-point crop-survey (March 2025), 15 crops | Stage 7b | Companion Zenodo dataset DOI |
| Sentinel-1 SAR, Sentinel-2 SR, ERA5-Land, SRTM | Stage 2 | Google Earth Engine (auto-fetched) |

## Where the non-bundled inputs actually go

They do **not** go in `data/`. The stage scripts read them from three folders at
the **repository root**, using these exact names:

```
<repo root>/
├── Analysis/
│   ├── Analysis_Banisuef_New_FF__March_2025.xls   # stage 1, sheet "total"
│   └── water_analyses Benisueif_FFF.xls           # stage 1, sheets "Total" and "IWQ_FF"
├── Layers/
│   ├── Study_Area_Last.shp                        # stages 2, 2b, 9
│   ├── Soil_Profiles.shp                          # stages 1, 2
│   └── Water_Samples.shp                          # stages 1, 2
└── Crop_March2025_11/
    └── Crop_March2025_11.shp                      # stage 7b
```

Every `.shp` must be accompanied by its `.shx`, `.dbf` and `.prj` siblings.
All three folders are listed in `.gitignore`, so they cannot be committed by
accident.

Verify placement with:

```
python notebooks/check_setup.py
```

## Skipping Stage 2

Stages 2 and 2b require a Google Cloud service-account key at
`.secrets/gee_service_account.json`. To avoid this, obtain the three
pre-computed artefacts from the corresponding author and place them at:

```
outputs/stage2/covariates_at_profiles.csv
outputs/stage2/covariates_at_water_samples.csv
outputs/stage2b_local_stack/covariate_stack_30m.tif
```

Stages 3–9 then run with no Google account. See [`SETUP.md`](../SETUP.md).

## Coordinate reference system

All field datasets must be in **EPSG:32636** (UTM 36N) before being passed to Stages 3–7. The pipeline will reproject lat/long inputs but not other CRSs.
