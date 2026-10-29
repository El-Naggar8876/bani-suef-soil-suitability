"""Convert the WaPOR v3 windows (fetched as Int16 .npy) into EPSG:4326 GeoTIFFs,
then reproject onto the 30 m EPSG:32636 analysis grid used by Stages 5-7."""
import json, numpy as np, rasterio, pathlib
from rasterio.transform import from_origin
from rasterio.warp import reproject, Resampling

import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parent))
from revpaths import REPO as _REPO, OUT as _OUT, WAPOR as _WAPOR

HERE = pathlib.Path(__file__).resolve().parent
RAW = _WAPOR / 'raw'
OUT  = HERE / "wapor"; OUT.mkdir(exist_ok=True)
REPO = _REPO
GRID = REPO / 'outputs' / 'stage5' / 'ales_index_q50_0_100cm.tif'

SCALE = {"T_2025": 0.1, "T_mean_2021_2025": 0.1, "AETI_mean_i16": 0.1, "TBP_mean_i16": 1.0}
UNITS = {"T_2025": "mm/year", "T_mean_2021_2025": "mm/year",
         "AETI_mean_i16": "mm/year", "TBP_mean_i16": "kg/ha"}
RENAME = {"AETI_mean_i16": "AETI_mean_2021_2025", "TBP_mean_i16": "TBP_mean_2021_2025"}

meta = {m["name"]: m for m in json.load(open(RAW / "meta.json"))}

with rasterio.open(GRID) as g:
    dst_prof = g.profile.copy()
    dst_crs, dst_tr, dH, dW = g.crs, g.transform, g.height, g.width
    valid_mask = np.isfinite(g.read(1))

for name, m in meta.items():
    a = np.load(RAW / f"{name}.npy").astype("float32")
    nod = m["nodata"]
    vals = np.where(a == nod, np.nan, a * m.get("div", 1) * SCALE[name]).astype("float32")
    tr = from_origin(m["west"], m["north"], m["res"], m["res"])
    out_name = RENAME.get(name, name)

    # 1. native EPSG:4326 GeoTIFF
    p4326 = OUT / f"wapor_{out_name}_4326.tif"
    with rasterio.open(p4326, "w", driver="GTiff", height=m["h"], width=m["w"], count=1,
                       dtype="float32", crs="EPSG:4326", transform=tr,
                       nodata=np.float32(np.nan), compress="deflate") as d:
        d.write(vals, 1); d.update_tags(1, UNIT=UNITS[name], SOURCE="FAO WaPOR v3 L2 100 m")

    # 2. reprojected onto the 30 m analysis grid (bilinear: continuous variable)
    dst = np.full((dH, dW), np.nan, "float32")
    reproject(vals, dst, src_transform=tr, src_crs="EPSG:4326",
              dst_transform=dst_tr, dst_crs=dst_crs,
              src_nodata=np.nan, dst_nodata=np.nan, resampling=Resampling.bilinear)
    dst[~valid_mask] = np.nan
    p32636 = OUT / f"wapor_{out_name}_30m.tif"
    prof = dict(dst_prof, dtype="float32", count=1, nodata=np.float32(np.nan), compress="deflate")
    with rasterio.open(p32636, "w", **prof) as d:
        d.write(dst, 1); d.update_tags(1, UNIT=UNITS[name], SOURCE="FAO WaPOR v3 L2 100 m")

    v = dst[np.isfinite(dst)]
    print(f"{out_name:22s} {UNITS[name]:8s} n={v.size:7d} min={v.min():8.1f} "
          f"med={np.median(v):8.1f} mean={v.mean():8.1f} max={v.max():8.1f}")
print("\nwrote:", *[p.name for p in sorted(OUT.glob('*.tif'))], sep="\n  ")
