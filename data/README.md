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

Once you have local copies, place them as:

```
data/
├── soil_profiles.csv          # or .xlsx — see stage1_data_audit.py for expected schema
├── water_samples.csv
├── crop_survey_2025_03.gpkg
└── (existing reference vectors)
```

If you only have the Earth Engine credentials, you can skip the manual datasets and run Stage 2 first; intermediate rasters will land in `outputs/stage2/`.

## Coordinate reference system

All field datasets must be in **EPSG:32636** (UTM 36N) before being passed to Stages 3–7. The pipeline will reproject lat/long inputs but not other CRSs.
