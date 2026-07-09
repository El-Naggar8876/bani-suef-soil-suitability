"""
Stage 2 — Climate & remote-sensing covariate stack
Beni Suef paper (target: Agric. Water Manag.)

Pipeline:
  1. Authenticate with Google Earth Engine (service account).
  2. Define AOI from the local Study_Area_Last.shp (UTM-36N -> WGS84 for GEE).
  3. Build covariates at 30-m resolution over AOI:
       Sentinel-2 SR (2024-2025), cloud-masked, dry & wet season medians
         -> bands B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12
         -> indices: NDVI, EVI, NDWI, MNDWI, BSI, SI, NDSI, SAVI, NBR
       Sentinel-1 SAR (VV, VH) annual median (texture context)
       SRTM DEM derivatives: elevation, slope, aspect (sin/cos), TPI, TWI
       CHIRPS rainfall annual mean
       TerraClimate ET0, AI annual mean
       ERA5-Land Tair annual mean
  4. Resample/reproject to AOI 30-m grid (UTM 36N).
  5. Export the stack as a single multi-band GeoTIFF (Drive or local).
  6. Sample the stack at all soil profiles and water samples; save CSVs.

Outputs in outputs/stage2/.
"""
from pathlib import Path
import os
import json
import time
import ee
import geopandas as gpd
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
LAYERS = ROOT / "Layers"
OUT = ROOT / "outputs" / "stage2"
STAGE1 = ROOT / "outputs" / "stage1"
OUT.mkdir(parents=True, exist_ok=True)

SA_KEY = ROOT / ".secrets" / "gee_service_account.json"
SA_EMAIL = json.loads(SA_KEY.read_text())["client_email"]
PROJECT = json.loads(SA_KEY.read_text())["project_id"]

ee.Initialize(
    credentials=ee.ServiceAccountCredentials(SA_EMAIL, str(SA_KEY)),
    project=PROJECT,
)
print("GEE initialised, project:", PROJECT)

# ---------------------------------------------------------------------------
# 1. AOI
# ---------------------------------------------------------------------------
aoi_gdf = gpd.read_file(LAYERS / "Study_Area_Last.shp", engine="pyogrio")
aoi_wgs = aoi_gdf.to_crs(4326)
geom = aoi_wgs.geometry.iloc[0]
# Build EE geometry
if geom.geom_type == "Polygon":
    coords = [list(geom.exterior.coords)]
else:
    coords = [list(p.exterior.coords) for p in geom.geoms]
AOI = ee.Geometry.MultiPolygon(coords)
print("AOI area (km^2):", AOI.area().divide(1e6).getInfo())

# Output projection: UTM 36N, 30 m
TARGET_CRS = "EPSG:32636"
SCALE = 30

# ---------------------------------------------------------------------------
# 2. Sentinel-2 SR seasonal medians (cloud masked)
# ---------------------------------------------------------------------------
def s2_mask(img):
    qa = img.select("QA60")
    cloud_bit_mask = 1 << 10
    cirrus_bit_mask = 1 << 11
    mask = (qa.bitwiseAnd(cloud_bit_mask).eq(0)
              .And(qa.bitwiseAnd(cirrus_bit_mask).eq(0)))
    scl = img.select("SCL")
    # Drop classes: 3 cloud shadow, 8 medium cloud, 9 high cloud, 10 thin cirrus
    scl_mask = (scl.neq(3).And(scl.neq(8)).And(scl.neq(9)).And(scl.neq(10)))
    return (img.updateMask(mask).updateMask(scl_mask)
               .divide(10000)
               .copyProperties(img, ["system:time_start"]))

S2_BANDS = ["B2","B3","B4","B5","B6","B7","B8","B8A","B11","B12"]

def s2_collection(start, end):
    return (ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
            .filterBounds(AOI)
            .filterDate(start, end)
            .filter(ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE", 40))
            .map(s2_mask)
            .select(S2_BANDS))

# Egypt seasons: winter (Nov-Apr) and summer (May-Oct)
s2_winter = s2_collection("2024-11-01", "2025-04-30").median()
s2_summer = s2_collection("2025-05-01", "2025-10-31").median()

def add_indices(img, suffix):
    b = img.select
    NDVI  = b("B8").subtract(b("B4")).divide(b("B8").add(b("B4"))).rename(f"NDVI_{suffix}")
    EVI   = b("B8").subtract(b("B4")).multiply(2.5).divide(
              b("B8").add(b("B4").multiply(6)).subtract(b("B2").multiply(7.5)).add(1)
            ).rename(f"EVI_{suffix}")
    NDWI  = b("B3").subtract(b("B8")).divide(b("B3").add(b("B8"))).rename(f"NDWI_{suffix}")
    MNDWI = b("B3").subtract(b("B11")).divide(b("B3").add(b("B11"))).rename(f"MNDWI_{suffix}")
    BSI   = (b("B11").add(b("B4")).subtract(b("B8").add(b("B2")))).divide(
              b("B11").add(b("B4")).add(b("B8")).add(b("B2"))
            ).rename(f"BSI_{suffix}")
    SI    = b("B11").subtract(b("B8")).divide(b("B11").add(b("B8"))).rename(f"SI_{suffix}")  # salinity index proxy
    NDSI  = b("B11").subtract(b("B12")).divide(b("B11").add(b("B12"))).rename(f"NDSI_{suffix}")
    SAVI  = b("B8").subtract(b("B4")).multiply(1.5).divide(
              b("B8").add(b("B4")).add(0.5)
            ).rename(f"SAVI_{suffix}")
    NBR   = b("B8").subtract(b("B12")).divide(b("B8").add(b("B12"))).rename(f"NBR_{suffix}")
    renamed = img.rename([f"{x}_{suffix}" for x in S2_BANDS])
    return renamed.addBands([NDVI, EVI, NDWI, MNDWI, BSI, SI, NDSI, SAVI, NBR])

s2_w = add_indices(s2_winter, "win")
s2_s = add_indices(s2_summer, "sum")

# ---------------------------------------------------------------------------
# 3. Sentinel-1 annual median VV/VH
# ---------------------------------------------------------------------------
s1 = (ee.ImageCollection("COPERNICUS/S1_GRD")
      .filterBounds(AOI)
      .filterDate("2024-11-01", "2025-10-31")
      .filter(ee.Filter.eq("instrumentMode", "IW"))
      .filter(ee.Filter.listContains("transmitterReceiverPolarisation","VV"))
      .filter(ee.Filter.listContains("transmitterReceiverPolarisation","VH"))
      .select(["VV","VH"])
      .median()
      .rename(["S1_VV","S1_VH"]))

# ---------------------------------------------------------------------------
# 4. DEM derivatives (SRTM 30 m)
# ---------------------------------------------------------------------------
dem = ee.Image("USGS/SRTMGL1_003").clip(AOI).rename("DEM")
slope = ee.Terrain.slope(dem).rename("slope")
aspect = ee.Terrain.aspect(dem)
asp_rad = aspect.multiply(3.14159265 / 180.0)
asp_sin = asp_rad.sin().rename("aspect_sin")
asp_cos = asp_rad.cos().rename("aspect_cos")
# TPI: difference from focal mean
focal_mean = dem.focalMean(radius=5, units="pixels")
tpi = dem.subtract(focal_mean).rename("TPI")
# TWI surrogate: ln(a / tan(slope+0.001)), with a = upslope area proxy via flow accumulation not in EE std.
# Use a simple wetness proxy: inverse slope.
twi_proxy = ee.Image(1).divide(slope.add(0.1)).log().rename("TWI_proxy")
dem_stack = dem.addBands([slope, asp_sin, asp_cos, tpi, twi_proxy])

# ---------------------------------------------------------------------------
# 5. Climate covariates
# ---------------------------------------------------------------------------
chirps = (ee.ImageCollection("UCSB-CHG/CHIRPS/DAILY")
          .filterDate("2015-01-01", "2024-12-31")
          .sum().divide(10).rename("rain_mm_yr_mean"))  # mean annual

terra = ee.ImageCollection("IDAHO_EPSCOR/TERRACLIMATE")\
        .filterDate("2015-01-01", "2024-12-31")
# TerraClimate pet/aet/def are stored in 0.1 mm units; tmmx/tmmn in 0.1 deg C
pet = terra.select("pet").mean().multiply(0.1).rename("PET_mm_mo")
aet = terra.select("aet").mean().multiply(0.1).rename("AET_mm_mo")
def_ = terra.select("def").mean().multiply(0.1).rename("water_deficit")
tmax = terra.select("tmmx").mean().multiply(0.1).rename("Tmax_C")
tmin = terra.select("tmmn").mean().multiply(0.1).rename("Tmin_C")
# Aridity index = annual P / annual PET (UNEP convention)
ai = chirps.divide(pet.multiply(12)).rename("AI")

clim_stack = chirps.addBands([pet, aet, def_, tmax, tmin, ai])

# ---------------------------------------------------------------------------
# 6. Combine and reproject
# ---------------------------------------------------------------------------
# Build stack groups separately. Do NOT call .reproject() before sampling -
# it forces full-stack rasterization at sample time and causes timeouts.
stack_s2 = s2_w.addBands(s2_s).clip(AOI)
stack_s1 = s1.clip(AOI)
stack_dem = dem_stack.clip(AOI)
stack_clim = clim_stack.clip(AOI)

stack = stack_s2.addBands(stack_s1).addBands(stack_dem).addBands(stack_clim)
band_names = stack.bandNames().getInfo()
print(f"Stack bands ({len(band_names)}):")
for b in band_names:
    print("  ", b)

# ---------------------------------------------------------------------------
# 7. Sample at soil profiles and water samples
# ---------------------------------------------------------------------------
sp = gpd.read_file(LAYERS / "Soil_Profiles.shp", engine="pyogrio").to_crs(4326)
ws = gpd.read_file(LAYERS / "Water_Samples.shp", engine="pyogrio").to_crs(4326)

def gdf_to_fc(gdf, id_col):
    feats = []
    for _, row in gdf.iterrows():
        feats.append(ee.Feature(
            ee.Geometry.Point([row.geometry.x, row.geometry.y]),
            {id_col: int(row[id_col])}
        ))
    return ee.FeatureCollection(feats)

sp_fc = gdf_to_fc(sp, "Profile")
ws_fc = gdf_to_fc(ws, "No")

def sample_grouped(points_fc, id_col):
    """Sample each band-group separately to avoid GEE memory/time-out.
    Sentinel-2 at 10 m, Sentinel-1 at 10 m, DEM at 30 m, climate at native.
    """
    out = {}
    print("  - Sentinel-2 ...")
    out['s2'] = stack_s2.sampleRegions(collection=points_fc, scale=10, geometries=False, tileScale=4).getInfo()
    print("  - Sentinel-1 ...")
    out['s1'] = stack_s1.sampleRegions(collection=points_fc, scale=10, geometries=False, tileScale=4).getInfo()
    print("  - DEM derivatives ...")
    out['dem'] = stack_dem.sampleRegions(collection=points_fc, scale=30, geometries=False, tileScale=4).getInfo()
    print("  - Climate ...")
    out['clim'] = stack_clim.sampleRegions(collection=points_fc, scale=1000, geometries=False, tileScale=4).getInfo()
    return out

def merge_groups(grouped, id_col):
    dfs = []
    for k, fc in grouped.items():
        rows = [f['properties'] for f in fc['features']]
        df = pd.DataFrame(rows)
        dfs.append(df.set_index(id_col))
    out = pd.concat(dfs, axis=1)
    out = out.loc[:, ~out.columns.duplicated()].reset_index()
    return out.sort_values(id_col).reset_index(drop=True)

print("Sampling stack at 60 profiles ...")
sp_samples = sample_grouped(sp_fc, 'Profile')
print("Sampling stack at 20 water samples ...")
ws_samples = sample_grouped(ws_fc, 'No')

def fc_to_df(grouped, id_col):
    return merge_groups(grouped, id_col)

sp_df = fc_to_df(sp_samples, "Profile")
ws_df = fc_to_df(ws_samples, "No")
print(f"Soil profile covariate rows: {len(sp_df)}, columns: {len(sp_df.columns)}")
print(f"Water sample covariate rows: {len(ws_df)}, columns: {len(ws_df.columns)}")

sp_df.to_csv(OUT / "covariates_at_profiles.csv", index=False)
ws_df.to_csv(OUT / "covariates_at_water_samples.csv", index=False)

# Quick sanity: print missing-value summary
miss = sp_df.isna().sum().sort_values(ascending=False).head(10)
print("\nTop missing-value columns at profiles:\n", miss)

# ---------------------------------------------------------------------------
# 8. Export the full stack to Google Drive (asynchronous task).
#     Direct download is too big for getDownloadURL (~600 km^2 x 53 bands).
# ---------------------------------------------------------------------------
stack_export = stack.reproject(crs=TARGET_CRS, scale=SCALE)
task = ee.batch.Export.image.toDrive(
    image=stack_export,
    description="benisuef_covariate_stack",
    folder="benisuef_paper",
    fileNamePrefix="covariate_stack_30m",
    region=AOI,
    scale=SCALE,
    crs=TARGET_CRS,
    maxPixels=1e10,
)
task.start()
print("\nDrive export task started:", task.id)
print("Track at https://code.earthengine.google.com/tasks")

# Save metadata
(OUT / "stack_band_names.txt").write_text("\n".join(band_names))
print("\nStage 2 outputs written to:", OUT)
