"""
Stage 9 v2: Advanced graphical abstract for AWM submission.

3-panel landscape (16 x 6 in @ 300 dpi) with:
  Panel 1 INPUT  : real Egypt outline + Nile + AOI bounding box; data-source list
  Panel 2 METHOD : 4-step flow + horizontal bar of CF deltas (Ks highlighted)
  Panel 3 OUTPUT : real priority raster overlaid with AOI outline + KPI strip
"""
from pathlib import Path
import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Polygon as MplPoly
import rasterio
import geopandas as gpd

ROOT = Path(__file__).resolve().parents[1]
FIG_DIR = ROOT / "outputs" / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

# palette
C_BG        = "white"
C_PANEL_BG  = "#FBFAF7"
C_TEXT      = "#222222"
C_MUTED     = "#6E6E6E"
C_LAND      = "#EFE7D6"
C_LAND_EDGE = "#9C9180"
C_NILE      = "#3F7AA6"
C_HL        = "#C0392B"
C_SOIL      = "#B6915C"
C_GREEN     = "#5E8C4A"

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 10,
    "axes.edgecolor": C_TEXT,
})

# assets
COUNTRIES = ROOT / "data" / "countries.geojson"
RIVERS    = ROOT / "data" / "ne_10m_rivers.geojson"
AOI_SHP   = ROOT / "Layers" / "Study_Area_Last.shp"
PRIO_TIF  = ROOT / "outputs" / "stage7" / "policy_priority_0_100cm.tif"

egypt_poly = None
if COUNTRIES.exists():
    data = json.load(open(COUNTRIES, encoding="utf-8"))
    feats = [f for f in data["features"] if f["properties"].get("name") == "Egypt"]
    if feats:
        geom = feats[0]["geometry"]
        if geom["type"] == "Polygon":
            egypt_poly = [np.array(geom["coordinates"][0])]
        else:
            egypt_poly = [np.array(p[0]) for p in geom["coordinates"]]

# Real Nile geometry from Natural Earth (main stem + Rosetta + Damietta delta
# branches) so that the Y-shaped bifurcation north of Cairo reads correctly.
nile_lines = []
nile_names = {"Nile", "Rosetta Branch", "Damietta Branch"}
if RIVERS.exists():
    rdata = json.load(open(RIVERS, encoding="utf-8"))
    for feat in rdata["features"]:
        name = (feat["properties"].get("name") or "").strip()
        if name in nile_names:
            geom = feat["geometry"]
            if geom["type"] == "LineString":
                nile_lines.append(np.array(geom["coordinates"]))
            elif geom["type"] == "MultiLineString":
                for seg in geom["coordinates"]:
                    nile_lines.append(np.array(seg))

aoi_bbox_ll = dict(lon_min=30.80, lon_max=31.42, lat_min=28.70, lat_max=29.42)
aoi_centre_ll = (0.5*(aoi_bbox_ll["lon_min"]+aoi_bbox_ll["lon_max"]),
                 0.5*(aoi_bbox_ll["lat_min"]+aoi_bbox_ll["lat_max"]))

aoi_gdf = gpd.read_file(AOI_SHP, engine="pyogrio") if AOI_SHP.exists() else None

# figure
fig = plt.figure(figsize=(16, 6.0), facecolor=C_BG)
gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.25, 1.05], wspace=0.04,
                      left=0.015, right=0.985, top=0.90, bottom=0.10)

fig.text(0.5, 0.965,
         "Counterfactual decomposition of irrigated land suitability — Beni Suef, Nile floodplain",
         ha="center", va="top", fontsize=12.5, color=C_TEXT, fontweight="bold")

# Panel 1
ax1 = fig.add_subplot(gs[0, 0]); ax1.set_xlim(0, 1); ax1.set_ylim(0, 1); ax1.axis("off")
ax1.add_patch(FancyBboxPatch((0.005, 0.005), 0.99, 0.99, boxstyle="round,pad=0.005",
                              facecolor=C_PANEL_BG, edgecolor="#DCD6C8", lw=1,
                              transform=ax1.transAxes))
ax1.text(0.5, 0.96, "1. INPUT  ·  multi-source DSM stack", ha="center", va="top",
         fontsize=11, fontweight="bold", color=C_TEXT)

egypt_ax = fig.add_axes([0.060, 0.45, 0.23, 0.34])
egypt_ax.set_aspect("equal")
egypt_ax.set_facecolor("#EAF1F6")
if egypt_poly is not None:
    for poly in egypt_poly:
        egypt_ax.add_patch(MplPoly(poly, closed=True, facecolor=C_LAND,
                                    edgecolor=C_LAND_EDGE, lw=0.8, zorder=2))
# Real Nile (main stem + Rosetta + Damietta) from Natural Earth
if nile_lines:
    for seg in nile_lines:
        egypt_ax.plot(seg[:, 0], seg[:, 1], color=C_NILE, lw=1.6, zorder=3,
                      solid_capstyle="round")
else:
    nile_lon = [32.90, 32.85, 32.50, 31.65, 31.20, 31.20, 31.40, 30.45, 31.55, 32.35]
    nile_lat = [24.10, 25.70, 27.20, 28.50, 29.80, 30.05, 30.40, 31.30, 31.55, 31.40]
    egypt_ax.plot(nile_lon, nile_lat, color=C_NILE, lw=1.6, zorder=3,
                   solid_capstyle="round")
egypt_ax.add_patch(MplPoly([[32.5,22.5],[33.0,22.0],[33.2,22.5],[32.7,23.5],[32.4,23.0]],
                            closed=True, facecolor=C_NILE, edgecolor=C_NILE, alpha=0.7, zorder=3))
egypt_ax.plot(31.24, 30.05, "o", color="black", ms=3.5, zorder=5)
egypt_ax.text(31.55, 30.20, "Cairo", fontsize=7.5, color=C_TEXT, zorder=5)

ax_w = aoi_bbox_ll["lon_max"] - aoi_bbox_ll["lon_min"]
ax_h = aoi_bbox_ll["lat_max"] - aoi_bbox_ll["lat_min"]
egypt_ax.add_patch(mpatches.Rectangle(
    (aoi_bbox_ll["lon_min"], aoi_bbox_ll["lat_min"]), ax_w, ax_h,
    facecolor=C_HL, edgecolor=C_HL, alpha=0.85, lw=0.5, zorder=6))
egypt_ax.annotate("AOI\n642 km²", xy=aoi_centre_ll, xytext=(28.0, 27.5),
                   fontsize=8, color=C_HL, fontweight="bold", ha="center",
                   arrowprops=dict(arrowstyle="-", color=C_HL, lw=0.8), zorder=7)

egypt_ax.set_xlim(24, 37); egypt_ax.set_ylim(21.5, 32.2)
egypt_ax.set_xticks([]); egypt_ax.set_yticks([])
for s in egypt_ax.spines.values():
    s.set_edgecolor(C_LAND_EDGE); s.set_linewidth(0.8)
egypt_ax.text(25.3, 31.6, "N", fontsize=8, fontweight="bold", color=C_TEXT, ha="center")
egypt_ax.annotate("", xy=(25.3, 31.5), xytext=(25.3, 30.7),
                   arrowprops=dict(arrowstyle="->", color=C_TEXT, lw=0.8))
egypt_ax.text(35.5, 22.0, "~600 km", fontsize=6.5, color=C_MUTED, ha="right")
egypt_ax.plot([34.5, 35.5], [22.3, 22.3], color=C_TEXT, lw=1)

items = [
    ("60",  "soil profiles · 0-30 / 30-60 / 60-100 cm", C_SOIL),
    ("20",  "irrigation-water samples · IWQI", C_NILE),
    ("53",  "covariates @ 30 m · Sentinel-2, SRTM, climate", C_GREEN),
    ("993", "crop-survey points · March 2025 · 15 crops", C_HL),
]
for i, (n, label, col) in enumerate(items):
    y = 0.36 - i * 0.072
    ax1.add_patch(FancyBboxPatch((0.06, y - 0.028), 0.88, 0.056,
                                  boxstyle="round,pad=0.005",
                                  facecolor=col, alpha=0.18,
                                  edgecolor=col, lw=1, transform=ax1.transAxes))
    ax1.text(0.11, y, n, fontsize=13, fontweight="bold", color=col,
             va="center", transform=ax1.transAxes)
    ax1.text(0.22, y, label, fontsize=8.4, color=C_TEXT,
             va="center", transform=ax1.transAxes)

# Panel 2
ax2 = fig.add_subplot(gs[0, 1]); ax2.set_xlim(0, 1); ax2.set_ylim(0, 1); ax2.axis("off")
ax2.add_patch(FancyBboxPatch((0.005, 0.005), 0.99, 0.99, boxstyle="round,pad=0.005",
                              facecolor=C_PANEL_BG, edgecolor="#DCD6C8", lw=1,
                              transform=ax2.transAxes))
ax2.text(0.5, 0.96, "2. METHOD  ·  probabilistic ALES + counterfactual engine",
         ha="center", va="top", fontsize=11, fontweight="bold", color=C_TEXT)

steps = [
    ("Quantile RF\n8 soil props", "#E6E0D2", C_TEXT),
    ("FAO–Sys ALES\nS1 → N1",     "#E6E0D2", C_TEXT),
    ("Counterfactual\ndo(X = x*)", "#FCE0DA", C_HL),
    ("Δ-suitability\nper limitation", "#FCE0DA", C_HL),
]
n = len(steps); box_w = 0.18; box_h = 0.11; y_flow = 0.80
xs = np.linspace(0.04, 0.96 - box_w, n)
for i, ((s, fc, ec), x) in enumerate(zip(steps, xs)):
    ax2.add_patch(FancyBboxPatch((x, y_flow - box_h/2), box_w, box_h,
                                  boxstyle="round,pad=0.006",
                                  facecolor=fc, edgecolor=ec, lw=1.3,
                                  transform=ax2.transAxes))
    ax2.text(x + box_w/2, y_flow, s, ha="center", va="center",
             fontsize=8.3, color=ec, transform=ax2.transAxes)
    if i < n - 1:
        ax2.annotate("", xy=(xs[i+1], y_flow), xytext=(x + box_w, y_flow),
                      xycoords=ax2.transAxes,
                      arrowprops=dict(arrowstyle="->", color=C_MUTED, lw=1.2))

ax_bar = fig.add_axes([0.395, 0.24, 0.27, 0.44])
levers = ["CF_full", "CF_K$_s$", "CF_chem", "CF_EC", "CF_SOM", "CF_ESP", "CF_CaCO$_3$"]
deltas = [24.85, 21.10, 2.10, 2.10, 0.31, 0.0, 0.0]
colours = ["#555555", C_HL, C_MUTED, C_MUTED, C_MUTED, C_MUTED, C_MUTED]
y_pos = np.arange(len(levers))[::-1]
ax_bar.barh(y_pos, deltas, color=colours, edgecolor="black", lw=0.4, height=0.72)
ax_bar.set_yticks(y_pos); ax_bar.set_yticklabels(levers, fontsize=8.5)
ax_bar.set_xlabel("Mean Δ ALES-index  (0–100 cm)", fontsize=8.5, color=C_TEXT)
ax_bar.set_xlim(0, 30)
ax_bar.tick_params(labelsize=8)
ax_bar.spines["top"].set_visible(False)
ax_bar.spines["right"].set_visible(False)
ax_bar.spines["left"].set_color(C_MUTED)
ax_bar.spines["bottom"].set_color(C_MUTED)
for yi, val in zip(y_pos, deltas):
    if val > 0.5:
        ax_bar.text(val + 0.4, yi, f"+{val:.1f}", va="center", fontsize=8, color=C_TEXT)
    else:
        ax_bar.text(0.4, yi, "≈ 0", va="center", fontsize=8, color=C_MUTED)
ks_y = y_pos[1]
ax_bar.annotate("85 % of\nclosable gap",
                 xy=(21.1, ks_y), xytext=(15.5, ks_y - 1.6),
                 fontsize=8.5, fontweight="bold", color=C_HL, ha="center",
                 arrowprops=dict(arrowstyle="->", color=C_HL, lw=1.0))

ax2.text(0.50, 0.06,
         "K$_s$ (a physical limit) — not salinity — drives the gap",
         ha="center", va="bottom", fontsize=9.2, style="italic",
         fontweight="bold", color=C_HL, transform=ax2.transAxes)

# Panel 3
ax3 = fig.add_subplot(gs[0, 2]); ax3.set_xlim(0, 1); ax3.set_ylim(0, 1); ax3.axis("off")
ax3.add_patch(FancyBboxPatch((0.005, 0.005), 0.99, 0.99, boxstyle="round,pad=0.005",
                              facecolor=C_PANEL_BG, edgecolor="#DCD6C8", lw=1,
                              transform=ax3.transAxes))
ax3.text(0.5, 0.96, "3. OUTPUT  ·  policy-priority surface + KPIs",
         ha="center", va="top", fontsize=11, fontweight="bold", color=C_TEXT)

map_ax = fig.add_axes([0.730, 0.45, 0.225, 0.34])
if PRIO_TIF.exists():
    with rasterio.open(PRIO_TIF) as src:
        arr = src.read(1)
        bounds = src.bounds
        crs = src.crs
        rgba = np.zeros((*arr.shape, 4), dtype=float)
        rgba[arr >= 0]   = [0.86, 0.86, 0.86, 1.0]
        rgba[arr == 2]   = [0.97, 0.62, 0.50, 1.0]
        rgba[arr == 1]   = [0.75, 0.18, 0.16, 1.0]
        rgba[arr <  0]   = [1, 1, 1, 0]
        rgba[arr == 255] = [1, 1, 1, 0]
        map_ax.imshow(rgba, extent=[bounds.left, bounds.right,
                                     bounds.bottom, bounds.top],
                       interpolation="nearest")
        map_ax.set_xlim(bounds.left, bounds.right)
        map_ax.set_ylim(bounds.bottom, bounds.top)
        if aoi_gdf is not None:
            aoi_proj = aoi_gdf.to_crs(crs)
            for geom in aoi_proj.geometry:
                geoms = [geom] if geom.geom_type == "Polygon" else list(geom.geoms)
                for g in geoms:
                    xs_, ys_ = g.exterior.xy
                    map_ax.plot(xs_, ys_, color=C_TEXT, lw=0.9, zorder=10)
        x0 = bounds.left + (bounds.right - bounds.left) * 0.05
        y0 = bounds.bottom + (bounds.top - bounds.bottom) * 0.05
        map_ax.plot([x0, x0 + 5000], [y0, y0], color="black", lw=2)
        map_ax.text(x0 + 2500, y0 + (bounds.top - bounds.bottom) * 0.02,
                     "5 km", fontsize=7, ha="center", color="black")
        nx_ = bounds.right - (bounds.right - bounds.left) * 0.06
        ny_ = bounds.bottom + (bounds.top - bounds.bottom) * 0.10
        map_ax.text(nx_, ny_ + (bounds.top-bounds.bottom)*0.06, "N",
                     fontsize=8, fontweight="bold", ha="center", color="black")
        map_ax.annotate("", xy=(nx_, ny_ + (bounds.top-bounds.bottom)*0.04),
                         xytext=(nx_, ny_ - (bounds.top-bounds.bottom)*0.02),
                         arrowprops=dict(arrowstyle="->", color="black", lw=1))
map_ax.set_xticks([]); map_ax.set_yticks([])
for s in map_ax.spines.values():
    s.set_edgecolor(C_MUTED); s.set_linewidth(0.6)
map_ax.set_aspect("equal")

leg_handles = [
    mpatches.Patch(facecolor="#C0392B", edgecolor="black", lw=0.4, label="HIGH"),
    mpatches.Patch(facecolor="#F79F80", edgecolor="black", lw=0.4, label="Medium"),
    mpatches.Patch(facecolor="#DCDCDC", edgecolor="black", lw=0.4, label="AOI"),
]
map_ax.legend(handles=leg_handles, loc="upper left", fontsize=6.8,
              frameon=True, framealpha=0.85, borderpad=0.3,
              handlelength=1.2, handletextpad=0.4)

kpis = [("85 %",       "of closable gap → K$_s$"),
        ("159 km²",    "HIGH-priority cropped land"),
        ("0.19 → 0.78", "P(≥ S2): baseline → full-CF"),
        ("OA 0.988",   "cropland mask vs. 993 ground-truth pts")]
for i, (val, label) in enumerate(kpis):
    y = 0.36 - i * 0.072
    ax3.add_patch(FancyBboxPatch((0.05, y - 0.030), 0.90, 0.060,
                                  boxstyle="round,pad=0.005",
                                  facecolor="white", edgecolor=C_HL, lw=0.8,
                                  transform=ax3.transAxes))
    ax3.text(0.085, y, val, fontsize=12, fontweight="bold",
             color=C_HL, va="center", transform=ax3.transAxes)
    ax3.text(0.46, y, label, fontsize=8.0, color=C_TEXT,
             va="center", ha="left", transform=ax3.transAxes)

for x_arrow in (0.346, 0.692):
    fig.patches.append(FancyArrowPatch(
        (x_arrow, 0.62), (x_arrow + 0.013, 0.62),
        transform=fig.transFigure, arrowstyle="-|>",
        mutation_scale=22, color=C_MUTED, lw=2.0))

fig.text(0.5, 0.020,
         "Fix infiltration first — not salinity — to unlock irrigated productivity in the Nile floodplain",
         ha="center", va="bottom", fontsize=10.5, style="italic", color=C_TEXT)

out_png = FIG_DIR / "graphical_abstract.png"
out_pdf = FIG_DIR / "graphical_abstract.pdf"
fig.savefig(out_png, dpi=300, facecolor=C_BG)
fig.savefig(out_pdf, facecolor=C_BG)
plt.close(fig)
print(f"Wrote {out_png.name}  ({out_png.stat().st_size/1024:.0f} KB)")
print(f"Wrote {out_pdf.name}  ({out_pdf.stat().st_size/1024:.0f} KB)")
