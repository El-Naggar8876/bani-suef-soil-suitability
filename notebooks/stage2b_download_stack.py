"""
Stage 2b - Download covariate stack to local disk (tiled GeoTIFFs).

Bypasses Drive entirely. Splits AOI into a grid of tiles, fetches each via
getDownloadURL (raw GeoTIFF zip), unzips into a single multi-band TIF per
tile, then mosaics tiles into one VRT (no resampling) and one COG.

Requires: gdal command line (or rasterio.merge if gdal not in PATH).
"""
from __future__ import annotations
from pathlib import Path
import json, io, sys, zipfile, math, time, urllib.request
import ee
import geopandas as gpd
import rasterio
from rasterio.merge import merge as rio_merge

ROOT = Path(__file__).resolve().parents[1]
LAYERS = ROOT / "Layers"
OUT = ROOT / "outputs" / "stage2b_local_stack"
TILES = OUT / "tiles"
OUT.mkdir(parents=True, exist_ok=True)
TILES.mkdir(parents=True, exist_ok=True)

# Earth Engine credentials are configured in config.py at the repository root.
# Edit that file, not this one. See SETUP.md, Appendix A.
sys.path.insert(0, str(ROOT))
from config import init_ee

init_ee()

# ---------------------------------------------------------------------------
# Re-build the same stack as Stage 2 but cast everything to Float32
# to avoid the export dtype-mismatch error.
# ---------------------------------------------------------------------------
aoi_gdf = gpd.read_file(LAYERS / "Study_Area_Last.shp", engine="pyogrio").to_crs(4326)
geom = aoi_gdf.geometry.iloc[0]
coords = ([list(geom.exterior.coords)] if geom.geom_type == "Polygon"
          else [list(p.exterior.coords) for p in geom.geoms])
AOI = ee.Geometry.MultiPolygon(coords)

TARGET_CRS = "EPSG:32636"
SCALE = 30

def s2_mask(img):
    qa = img.select("QA60")
    m = (qa.bitwiseAnd(1 << 10).eq(0).And(qa.bitwiseAnd(1 << 11).eq(0)))
    scl = img.select("SCL")
    sm = scl.neq(3).And(scl.neq(8)).And(scl.neq(9)).And(scl.neq(10))
    return (img.updateMask(m).updateMask(sm).divide(10000)
            .copyProperties(img, ["system:time_start"]))

S2_BANDS = ["B2","B3","B4","B5","B6","B7","B8","B8A","B11","B12"]

def s2_collection(start, end):
    return (ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
            .filterBounds(AOI).filterDate(start, end)
            .filter(ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE", 40))
            .map(s2_mask).select(S2_BANDS))

def add_indices(img, suffix):
    b = img.select
    NDVI  = b("B8").subtract(b("B4")).divide(b("B8").add(b("B4"))).rename(f"NDVI_{suffix}")
    EVI   = b("B8").subtract(b("B4")).multiply(2.5).divide(
              b("B8").add(b("B4").multiply(6)).subtract(b("B2").multiply(7.5)).add(1)
            ).rename(f"EVI_{suffix}")
    NDWI  = b("B3").subtract(b("B8")).divide(b("B3").add(b("B8"))).rename(f"NDWI_{suffix}")
    MNDWI = b("B3").subtract(b("B11")).divide(b("B3").add(b("B11"))).rename(f"MNDWI_{suffix}")
    BSI   = (b("B11").add(b("B4")).subtract(b("B8").add(b("B2")))).divide(
              b("B11").add(b("B4")).add(b("B8")).add(b("B2"))).rename(f"BSI_{suffix}")
    SI    = b("B11").subtract(b("B8")).divide(b("B11").add(b("B8"))).rename(f"SI_{suffix}")
    NDSI  = b("B11").subtract(b("B12")).divide(b("B11").add(b("B12"))).rename(f"NDSI_{suffix}")
    SAVI  = b("B8").subtract(b("B4")).multiply(1.5).divide(
              b("B8").add(b("B4")).add(0.5)).rename(f"SAVI_{suffix}")
    NBR   = b("B8").subtract(b("B12")).divide(b("B8").add(b("B12"))).rename(f"NBR_{suffix}")
    renamed = img.rename([f"{x}_{suffix}" for x in S2_BANDS])
    return renamed.addBands([NDVI, EVI, NDWI, MNDWI, BSI, SI, NDSI, SAVI, NBR])

s2_w = add_indices(s2_collection("2024-11-01", "2025-04-30").median(), "win")
s2_s = add_indices(s2_collection("2025-05-01", "2025-10-31").median(), "sum")

s1 = (ee.ImageCollection("COPERNICUS/S1_GRD")
      .filterBounds(AOI).filterDate("2024-11-01", "2025-10-31")
      .filter(ee.Filter.eq("instrumentMode", "IW"))
      .filter(ee.Filter.listContains("transmitterReceiverPolarisation","VV"))
      .filter(ee.Filter.listContains("transmitterReceiverPolarisation","VH"))
      .select(["VV","VH"]).median().rename(["S1_VV","S1_VH"]))

dem = ee.Image("USGS/SRTMGL1_003").rename("DEM")
slope = ee.Terrain.slope(dem).rename("slope")
asp = ee.Terrain.aspect(dem).multiply(3.14159265 / 180.0)
asp_sin = asp.sin().rename("aspect_sin")
asp_cos = asp.cos().rename("aspect_cos")
tpi = dem.subtract(dem.focalMean(radius=5, units="pixels")).rename("TPI")
twi = ee.Image(1).divide(slope.add(0.1)).log().rename("TWI_proxy")
dem_stack = dem.addBands([slope, asp_sin, asp_cos, tpi, twi])

chirps = (ee.ImageCollection("UCSB-CHG/CHIRPS/DAILY")
          .filterDate("2015-01-01", "2024-12-31").sum().divide(10)
          .rename("rain_mm_yr_mean"))
terra = (ee.ImageCollection("IDAHO_EPSCOR/TERRACLIMATE")
         .filterDate("2015-01-01", "2024-12-31"))
pet = terra.select("pet").mean().multiply(0.1).rename("PET_mm_mo")
aet = terra.select("aet").mean().multiply(0.1).rename("AET_mm_mo")
def_ = terra.select("def").mean().multiply(0.1).rename("water_deficit")
tmax = terra.select("tmmx").mean().multiply(0.1).rename("Tmax_C")
tmin = terra.select("tmmn").mean().multiply(0.1).rename("Tmin_C")
ai = chirps.divide(pet.multiply(12)).rename("AI")
clim = chirps.addBands([pet, aet, def_, tmax, tmin, ai])

stack = (s2_w.addBands(s2_s).addBands(s1).addBands(dem_stack).addBands(clim)
         .clip(AOI).toFloat())  # cast everything to Float32 for export
band_names = stack.bandNames().getInfo()
print(f"Stack bands ({len(band_names)})")

# ---------------------------------------------------------------------------
# Tile AOI in EPSG:32636 and download each tile via getDownloadURL.
# ---------------------------------------------------------------------------
aoi_utm = aoi_gdf.to_crs(32636)
minx, miny, maxx, maxy = aoi_utm.total_bounds
print(f"AOI UTM bbox: {minx:.0f} {miny:.0f} -> {maxx:.0f} {maxy:.0f}")

TILE_SIZE_M = 8000   # 8 km tiles -> ~267x267 px @ 30 m -> ~71k px x 53 bands x 4 B = ~15 MB raw
nx = math.ceil((maxx - minx) / TILE_SIZE_M)
ny = math.ceil((maxy - miny) / TILE_SIZE_M)
print(f"Grid: {nx} x {ny} = {nx*ny} tiles")

aoi_geom_utm = aoi_utm.geometry.iloc[0]

tile_files = []
for j in range(ny):
    for i in range(nx):
        x0 = minx + i * TILE_SIZE_M
        y0 = miny + j * TILE_SIZE_M
        x1 = min(x0 + TILE_SIZE_M, maxx)
        y1 = min(y0 + TILE_SIZE_M, maxy)
        # Skip tile if it doesn't intersect AOI polygon
        from shapely.geometry import box
        tile_box = box(x0, y0, x1, y1)
        if not tile_box.intersects(aoi_geom_utm):
            continue
        out_tif = TILES / f"tile_{j:02d}_{i:02d}.tif"
        if out_tif.exists() and out_tif.stat().st_size > 1000:
            tile_files.append(out_tif)
            continue
        region = ee.Geometry.Rectangle([x0, y0, x1, y1], proj=TARGET_CRS, geodesic=False)
        try:
            url = stack.getDownloadURL({
                "scale": SCALE, "region": region, "crs": TARGET_CRS,
                "fileFormat": "GeoTIFF", "filePerBand": False,
            })
        except Exception as e:
            print(f"  tile {j},{i} URL failed: {e}")
            continue
        # Download zip
        try:
            with urllib.request.urlopen(url, timeout=300) as r:
                payload = r.read()
        except Exception as e:
            print(f"  tile {j},{i} download failed: {e}")
            continue
        z = zipfile.ZipFile(io.BytesIO(payload))
        names = [n for n in z.namelist() if n.endswith(".tif")]
        if not names:
            print(f"  tile {j},{i} no tif in zip")
            continue
        out_tif.write_bytes(z.read(names[0]))
        size_mb = out_tif.stat().st_size / 1e6
        print(f"  tile {j:02d}_{i:02d} -> {size_mb:.1f} MB")
        tile_files.append(out_tif)
        time.sleep(0.5)  # be polite to GEE

print(f"\nDownloaded {len(tile_files)} tiles. Mosaicking ...")

# ---------------------------------------------------------------------------
# Mosaic tiles into a single multi-band GeoTIFF
# ---------------------------------------------------------------------------
srcs = [rasterio.open(p) for p in tile_files]
mosaic, transform = rio_merge(srcs)
profile = srcs[0].profile.copy()
profile.update({
    "height": mosaic.shape[1], "width": mosaic.shape[2],
    "transform": transform,
    "compress": "deflate", "predictor": 3, "tiled": True,
    "blockxsize": 512, "blockysize": 512, "BIGTIFF": "IF_NEEDED",
})
out_tif = OUT / "covariate_stack_30m.tif"
with rasterio.open(out_tif, "w", **profile) as dst:
    dst.write(mosaic)
    dst.descriptions = tuple(band_names)
for s in srcs:
    s.close()

print(f"\nWrote: {out_tif}  ({out_tif.stat().st_size/1e6:.1f} MB, {len(band_names)} bands)")
(OUT / "band_names.txt").write_text("\n".join(band_names))
print("Band names list:", OUT / "band_names.txt")
