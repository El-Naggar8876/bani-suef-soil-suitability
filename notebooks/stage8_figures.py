"""
Stage 8 - Publication-quality figures for the AWM manuscript.

Produces:
    fig1_study_area.png            - AOI on Sentinel-2 RGB w/ samples (placeholder if no S2 RGB)
    fig2_soil_uncertainty.png      - 6 panels: q05/q50/q95 of Ks and EC at 0-100 cm
    fig3_baseline_ales.png         - 2x3: q50 / class / P(>=S2), both depths
    fig4_counterfactuals.png       - 2x3 delta_q50 maps + bar chart
    fig5_cropping_priority.png     - landuse + contingency heatmap + priority
    fig6_workflow.png              - schematic (boxes + arrows, matplotlib)

Notes:
- Uses matplotlib only; no contextily/cartopy dependencies.
- All maps overlay AOI extent boundary; nodata masked.
- Output dir: outputs/figures/  (300 DPI PNG + PDF for vector versions of f4/f6).
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle
from matplotlib.colors import TwoSlopeNorm, LogNorm
import matplotlib.patches as mpatches
import rasterio
import geopandas as gpd

ROOT     = Path(__file__).resolve().parents[1]
S1       = ROOT / "outputs" / "stage1"
S2B      = ROOT / "outputs" / "stage2b_local_stack"
S3B      = ROOT / "outputs" / "stage3b"
S5       = ROOT / "outputs" / "stage5"
S6       = ROOT / "outputs" / "stage6"
S7       = ROOT / "outputs" / "stage7"
OUT      = ROOT / "outputs" / "figures"
OUT.mkdir(parents=True, exist_ok=True)

mpl.rcParams.update({
    "figure.dpi": 120,
    "savefig.dpi": 300,
    "font.size": 9,
    "axes.titlesize": 10,
    "axes.labelsize": 9,
    "axes.linewidth": 0.8,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8,
    "font.family": "DejaVu Sans",
})

ALES_COLORS = ["#1a9850", "#a6d96a", "#fee08b", "#fdae61", "#d73027"]
ALES_LABELS = ["S1", "S2", "S3", "N1", "N2"]
ALES_CMAP = ListedColormap(ALES_COLORS)
ALES_NORM = BoundaryNorm([0.5, 1.5, 2.5, 3.5, 4.5, 5.5], ALES_CMAP.N)

LU_COLORS = ["#2b83ba", "#1a9850", "#a6d96a", "#d9d9d9", "#fee08b"]
LU_LABELS = ["water", "double crop", "single crop", "bare/desert", "sparse veg."]
LU_CMAP = ListedColormap(LU_COLORS)
LU_NORM = BoundaryNorm([0.5, 1.5, 2.5, 3.5, 4.5, 5.5], LU_CMAP.N)


def read_masked(path: Path):
    with rasterio.open(path) as r:
        a = r.read(1).astype(np.float32)
        nd = r.nodata
        if nd is not None:
            a[a == nd] = np.nan
        return np.ma.masked_invalid(a), r.bounds, r.transform


def add_scale_bar(ax, length_m=10000, label="10 km"):
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    x0 = xmin + 0.05 * (xmax - xmin)
    y0 = ymin + 0.06 * (ymax - ymin)
    ax.plot([x0, x0 + length_m], [y0, y0], color="k", lw=2.5, solid_capstyle="butt")
    ax.text(x0 + length_m / 2, y0 + 0.012 * (ymax - ymin), label,
            ha="center", va="bottom", fontsize=7)


def add_north_arrow(ax):
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    x = xmin + 0.95 * (xmax - xmin)
    y = ymin + 0.92 * (ymax - ymin)
    dy = 0.06 * (ymax - ymin)
    ax.annotate("N", xy=(x, y), xytext=(x, y - dy),
                arrowprops=dict(arrowstyle="-|>", color="k", lw=1.2),
                ha="center", va="center", fontsize=9, fontweight="bold")


# ============================================================
# Figure 1 - Study area
# ============================================================
def fig1_study_area():
    """Two-panel cartographic figure:
      Left  : Egypt locator with real country outline + Nile + AOI bbox
      Right : AOI elevation map with soil/water sampling points,
              scale bar, north arrow, lat/lon graticule.
    """
    import json
    from matplotlib.patches import Polygon as MplPoly

    soil = gpd.read_file(S1 / "soil_profiles_0_100cm.gpkg", engine="pyogrio")
    water = gpd.read_file(S1 / "water_samples_iwqi.gpkg", engine="pyogrio") \
        if (S1 / "water_samples_iwqi.gpkg").exists() else None

    # ---- elevation backdrop
    elev_path = S2B / "covariate_stack_30m.tif"
    bands = (S2B / "band_names.txt").read_text().splitlines()
    with rasterio.open(elev_path) as r:
        elev = r.read(bands.index("DEM") + 1).astype(np.float32)
        if r.nodata is not None:
            elev[elev == r.nodata] = np.nan
        ext = (r.bounds.left, r.bounds.right, r.bounds.bottom, r.bounds.top)
        crs_proj = r.crs
        bounds = r.bounds
    elev = np.ma.masked_invalid(elev)
    elev_vals = elev.compressed()

    # ---- Egypt outline
    countries_geojson = ROOT / "data" / "countries.geojson"
    egypt_polys = None
    if countries_geojson.exists():
        with open(countries_geojson, encoding="utf-8") as f:
            data = json.load(f)
        feats = [f for f in data["features"]
                 if f["properties"].get("name") == "Egypt"]
        if feats:
            geom = feats[0]["geometry"]
            if geom["type"] == "Polygon":
                egypt_polys = [np.array(geom["coordinates"][0])]
            else:
                egypt_polys = [np.array(p[0]) for p in geom["coordinates"]]

    # ---- Real Nile geometry from Natural Earth (main stem + Delta branches).
    # Natural Earth labels the two delta distributaries separately:
    # "Rosetta Branch" (western) and "Damietta Branch" (eastern); we plot
    # them in the same blue as the main "Nile" feature so the Y-shaped
    # bifurcation north of Cairo reads as part of the river system.
    nile_names = {"Nile", "Rosetta Branch", "Damietta Branch"}
    nile_lines = []
    rivers_geojson = ROOT / "data" / "ne_10m_rivers.geojson"
    if rivers_geojson.exists():
        with open(rivers_geojson, encoding="utf-8") as f:
            rdata = json.load(f)
        for feat in rdata["features"]:
            name = (feat["properties"].get("name") or "").strip()
            if name in nile_names:
                geom = feat["geometry"]
                if geom["type"] == "LineString":
                    nile_lines.append(np.array(geom["coordinates"]))
                elif geom["type"] == "MultiLineString":
                    for seg in geom["coordinates"]:
                        nile_lines.append(np.array(seg))

    # ---- AOI bbox in lat/lon (transform raster bounds to EPSG:4326)
    try:
        from pyproj import Transformer
        tr = Transformer.from_crs(crs_proj.to_epsg() or 32636, 4326,
                                   always_xy=True)
        lon_min, lat_min = tr.transform(bounds.left,  bounds.bottom)
        lon_max, lat_max = tr.transform(bounds.right, bounds.top)
    except Exception:
        # fallback to known AOI envelope
        lon_min, lat_min, lon_max, lat_max = 30.80, 28.70, 31.42, 29.42
    aoi_centre_ll = (0.5 * (lon_min + lon_max), 0.5 * (lat_min + lat_max))

    # ============================================================ figure
    fig = plt.figure(figsize=(13, 7.2), facecolor="white")
    gs = fig.add_gridspec(1, 2, width_ratios=[0.85, 1.55], wspace=0.18,
                          left=0.06, right=0.94, top=0.92, bottom=0.10)

    # ---- LEFT: Egypt locator
    ax_loc = fig.add_subplot(gs[0, 0])
    ax_loc.set_aspect("equal")
    ax_loc.set_facecolor("#EAF1F6")  # sea/desert backdrop
    if egypt_polys is not None:
        for poly in egypt_polys:
            ax_loc.add_patch(MplPoly(poly, closed=True,
                                      facecolor="#EFE7D6",
                                      edgecolor="#9C9180", lw=1.0, zorder=2))
    # Real Nile course (Natural Earth)
    if nile_lines:
        for seg in nile_lines:
            ax_loc.plot(seg[:, 0], seg[:, 1], color="#3F7AA6", lw=1.7,
                         zorder=3, solid_capstyle="round")
    else:
        # Fallback schematic if data missing
        nile_lon = [32.90, 32.85, 32.50, 31.65, 31.20, 31.20, 31.40, 30.45, 31.55, 32.35]
        nile_lat = [24.10, 25.70, 27.20, 28.50, 29.80, 30.05, 30.40, 31.30, 31.55, 31.40]
        ax_loc.plot(nile_lon, nile_lat, color="#3F7AA6", lw=1.7, zorder=3,
                    solid_capstyle="round")
    # Lake Nasser
    ax_loc.add_patch(MplPoly([[32.5, 22.5], [33.0, 22.0], [33.2, 22.5],
                              [32.7, 23.5], [32.4, 23.0]],
                             closed=True, facecolor="#3F7AA6",
                             edgecolor="#3F7AA6", alpha=0.7, zorder=3))
    # Cairo dot
    ax_loc.plot(31.24, 30.05, "o", color="black", ms=4, zorder=5)
    ax_loc.text(31.55, 30.20, "Cairo", fontsize=9, color="#222222", zorder=5)
    # AOI red bbox (drawn at true geographic location)
    ax_loc.add_patch(plt.Rectangle(
        (lon_min, lat_min), lon_max - lon_min, lat_max - lat_min,
        facecolor="#C0392B", edgecolor="#C0392B", alpha=0.85,
        lw=0.8, zorder=6))
    ax_loc.annotate("Beni Suef\nAOI",
                    xy=aoi_centre_ll, xytext=(28.0, 27.4),
                    fontsize=10, color="#C0392B", fontweight="bold",
                    ha="center",
                    arrowprops=dict(arrowstyle="-", color="#C0392B", lw=1),
                    zorder=7)
    ax_loc.set_xlim(24, 37); ax_loc.set_ylim(21.5, 32.2)
    ax_loc.set_xlabel("Longitude (°E)", fontsize=9)
    ax_loc.set_ylabel("Latitude (°N)", fontsize=9)
    ax_loc.tick_params(labelsize=8)
    ax_loc.set_title("(a) Locator — Egypt", fontsize=11, fontweight="bold",
                     loc="left", pad=6)
    # north arrow
    ax_loc.annotate("N", xy=(25.3, 31.7), fontsize=10, fontweight="bold",
                    color="black", ha="center")
    ax_loc.annotate("", xy=(25.3, 31.55), xytext=(25.3, 30.7),
                    arrowprops=dict(arrowstyle="->", color="black", lw=1))
    # scale (degrees → ~111 km @ this latitude; 5° ≈ 480 km)
    ax_loc.plot([34.5, 35.5], [22.2, 22.2], color="black", lw=2)
    ax_loc.text(35.0, 22.5, "~100 km", fontsize=7.5, ha="center",
                color="#333333")

    # ---- RIGHT: AOI map
    ax = fig.add_subplot(gs[0, 1])
    im = ax.imshow(elev, extent=ext, cmap="terrain", origin="upper",
                   vmin=np.percentile(elev_vals, 2),
                   vmax=np.percentile(elev_vals, 98))
    soil = soil.to_crs(crs_proj)
    ax.scatter(soil.geometry.x, soil.geometry.y, s=26, c="#C0392B",
               edgecolor="white", lw=0.6,
               label=f"Soil profiles (n={len(soil)})", zorder=5)
    if water is not None:
        water = water.to_crs(crs_proj)
        ax.scatter(water.geometry.x, water.geometry.y, s=34, c="#1F6FA3",
                   marker="^", edgecolor="white", lw=0.6,
                   label=f"Water samples (n={len(water)})", zorder=5)
    ax.set_xlabel("Easting (m, UTM 36 N)", fontsize=10)
    ax.set_ylabel("Northing (m, UTM 36 N)", fontsize=10)
    ax.set_title("(b) Beni Suef floodplain — sampling design",
                 fontsize=11, fontweight="bold", loc="left", pad=6)
    leg = ax.legend(loc="upper left", framealpha=0.92, fontsize=9,
                    edgecolor="#888888")
    leg.set_zorder(20)

    # colorbar (vertical, attached to right)
    cbar = plt.colorbar(im, ax=ax, fraction=0.045, pad=0.02)
    cbar.set_label("Elevation (m a.s.l.)", fontsize=10)
    cbar.ax.tick_params(labelsize=8)
    add_scale_bar(ax)
    add_north_arrow(ax)

    plt.savefig(OUT / "fig1_study_area.png", dpi=300, bbox_inches="tight",
                facecolor="white")
    plt.savefig(OUT / "fig1_study_area.pdf", bbox_inches="tight",
                facecolor="white")
    plt.close()
    print("Wrote fig1_study_area.{png,pdf}")


# ============================================================
# Figure 2 - Soil property uncertainty (Ks and EC, 0-100cm)
# ============================================================
def fig2_soil_uncertainty():
    panels = [
        ("Ks q05",  S3B / "Ks_0_100cm_q05_qrf.tif", "viridis", "Ks (cm h$^{-1}$)"),
        ("Ks q50",  S3B / "Ks_0_100cm_q50_idw.tif", "viridis", "Ks (cm h$^{-1}$)"),
        ("Ks q95",  S3B / "Ks_0_100cm_q95_qrf.tif", "viridis", "Ks (cm h$^{-1}$)"),
        ("EC q05",  S3B / "EC_0_100cm_q05_qrf.tif", "magma",   "EC (dS m$^{-1}$)"),
        ("EC q50",  S3B / "EC_0_100cm_q50_idw.tif", "magma",   "EC (dS m$^{-1}$)"),
        ("EC q95",  S3B / "EC_0_100cm_q95_qrf.tif", "magma",   "EC (dS m$^{-1}$)"),
    ]
    fig, axes = plt.subplots(2, 3, figsize=(11, 7.2), constrained_layout=True)
    for ax, (title, path, cmap, label) in zip(axes.flat, panels):
        a, b, _ = read_masked(path)
        ext = (b.left, b.right, b.bottom, b.top)
        vals = a.compressed()
        vmin = np.percentile(vals, 2)
        vmax = np.percentile(vals, 98)
        im = ax.imshow(a, extent=ext, origin="upper", cmap=cmap, vmin=vmin, vmax=vmax)
        ax.set_title(title)
        ax.set_xticks([]); ax.set_yticks([])
        cbar = plt.colorbar(im, ax=ax, fraction=0.045, pad=0.02)
        cbar.set_label(label, fontsize=8)
    plt.savefig(OUT / "fig2_soil_uncertainty.png", bbox_inches="tight")
    plt.close()
    print("Wrote fig2_soil_uncertainty.png")


# ============================================================
# Figure 3 - Baseline ALES suitability
# ============================================================
def fig3_baseline_ales():
    fig, axes = plt.subplots(2, 3, figsize=(11, 7.2), constrained_layout=True)
    for row, depth in enumerate(["0_30cm", "0_100cm"]):
        idx, b, _ = read_masked(S5 / f"ales_index_q50_{depth}.tif")
        cls, _, _ = read_masked(S5 / f"ales_class_q50_{depth}.tif")
        ps2, _, _ = read_masked(S5 / f"p_at_least_S2_{depth}.tif")
        ext = (b.left, b.right, b.bottom, b.top)

        ax = axes[row, 0]
        idx_vals = idx.compressed()
        vmax = float(np.percentile(idx_vals, 98)) if idx_vals.size else 75
        im = ax.imshow(idx, extent=ext, origin="upper", cmap="RdYlGn",
                       vmin=0, vmax=max(vmax, 30))
        ax.set_title(f"ALES index q50 — {depth.replace('_', '–')}")
        ax.set_xticks([]); ax.set_yticks([])
        plt.colorbar(im, ax=ax, fraction=0.045, pad=0.02, label="index")

        ax = axes[row, 1]
        im = ax.imshow(cls, extent=ext, origin="upper", cmap=ALES_CMAP, norm=ALES_NORM)
        ax.set_title(f"Modal class — {depth.replace('_', '–')}")
        ax.set_xticks([]); ax.set_yticks([])
        cbar = plt.colorbar(im, ax=ax, fraction=0.045, pad=0.02,
                            ticks=[1, 2, 3, 4, 5])
        cbar.ax.set_yticklabels(ALES_LABELS)

        ax = axes[row, 2]
        im = ax.imshow(ps2, extent=ext, origin="upper", cmap="Blues",
                       vmin=0, vmax=1)
        ax.set_title(f"P(class ≥ S2) — {depth.replace('_', '–')}")
        ax.set_xticks([]); ax.set_yticks([])
        plt.colorbar(im, ax=ax, fraction=0.045, pad=0.02, label="probability")

    plt.savefig(OUT / "fig3_baseline_ales.png", bbox_inches="tight")
    plt.close()
    print("Wrote fig3_baseline_ales.png")


# ============================================================
# Figure 4 - Counterfactual decomposition
# ============================================================
def fig4_counterfactuals():
    scenarios = ["CF_EC", "CF_ESP", "CF_CaCO3", "CF_Ks", "CF_SOM"]
    fig = plt.figure(figsize=(13, 8.5))
    gs = fig.add_gridspec(2, 4, height_ratios=[1, 0.9], hspace=0.25, wspace=0.2)

    # 5 delta maps at 0-100 cm
    vmax_global = 35
    for i, scen in enumerate(scenarios):
        ax = fig.add_subplot(gs[i // 4 * 0 + (0 if i < 4 else 1),  # row
                                 i % 4 if i < 4 else 0])
        # Simpler: explicit positions
    # rebuild with explicit grid:
    plt.close(fig)

    fig, axes = plt.subplots(2, 3, figsize=(13, 8.4), constrained_layout=True)
    axes_list = list(axes.flat)
    # All counterfactual gains are >= 0 (relaxing a constraint cannot decrease
    # suitability), so a sequential Reds ramp gives a faithful legend mapping
    # 0 -> white, vmax -> dark red, with linearly distributed ticks.
    for ax, scen in zip(axes_list[:5], scenarios):
        a, b, _ = read_masked(S6 / f"cf_{scen}_0_100cm_delta_q50.tif")
        ext = (b.left, b.right, b.bottom, b.top)
        im = ax.imshow(a, extent=ext, origin="upper", cmap="Reds",
                       vmin=0, vmax=vmax_global)
        nice = scen.replace("CF_", "Δ ")
        ax.set_title(f"{nice}  (0–100 cm)")
        ax.set_xticks([]); ax.set_yticks([])
        plt.colorbar(im, ax=ax, fraction=0.045, pad=0.02, label="Δ index q50")

    # Bar chart in 6th panel
    ax = axes_list[5]
    g30  = pd.read_csv(S6 / "gap_decomposition_0_30cm.csv")
    g100 = pd.read_csv(S6 / "gap_decomposition_0_100cm.csv")
    order = ["CF_EC", "CF_ESP", "CF_CaCO3", "CF_SOM", "CF_chem", "CF_Ks", "CF_full"]
    d30  = [g30 .set_index("scenario").loc[s, "mean_delta_index"] for s in order]
    d100 = [g100.set_index("scenario").loc[s, "mean_delta_index"] for s in order]
    x = np.arange(len(order))
    w = 0.4
    ax.bar(x - w / 2, d30,  w, label="0–30 cm",  color="#fdae61")
    ax.bar(x + w / 2, d100, w, label="0–100 cm", color="#d73027")
    ax.set_xticks(x)
    ax.set_xticklabels([s.replace("CF_", "") for s in order], rotation=30, ha="right")
    ax.set_ylabel("Mean Δ ALES index q50")
    ax.legend(frameon=False)
    ax.grid(axis="y", alpha=0.3)

    plt.savefig(OUT / "fig4_counterfactuals.png", bbox_inches="tight")
    plt.savefig(OUT / "fig4_counterfactuals.pdf", bbox_inches="tight")
    plt.close()
    print("Wrote fig4_counterfactuals.{png,pdf}")


# ============================================================
# Figure 5 - Cropping mismatch & priority
# ============================================================
def fig5_cropping_priority():
    fig, axes = plt.subplots(1, 3, figsize=(14, 5.2), constrained_layout=True,
                             gridspec_kw={"width_ratios": [1, 0.85, 1]})

    lu, b, _ = read_masked(S7 / "landuse_class.tif")
    ext = (b.left, b.right, b.bottom, b.top)
    ax = axes[0]
    im = ax.imshow(lu, extent=ext, origin="upper", cmap=LU_CMAP, norm=LU_NORM)
    ax.set_title("Sentinel-2 land-use classification")
    ax.set_xticks([]); ax.set_yticks([])
    cbar = plt.colorbar(im, ax=ax, fraction=0.045, pad=0.02, ticks=[1, 2, 3, 4, 5])
    cbar.ax.set_yticklabels(LU_LABELS)

    # Crop x ALES contingency (independent ground-truth, n = 993, March 2025).
    # Replaces the previous LU x ALES pixel-count heatmap with the named-crop
    # version that uses survey points as the row dimension.
    s7b = ROOT / "outputs" / "stage7b" / "crop_x_ales_class_0_100cm.csv"
    ct = pd.read_csv(s7b, index_col=0)
    # Sort crops by total survey n (descending) so dominant rotations sit on top
    ct = ct.loc[ct.sum(axis=1).sort_values(ascending=False).index]
    ax = axes[1]
    arr = ct.values.astype(float)
    arr_log = np.log10(arr + 1)
    im = ax.imshow(arr_log, cmap="Greys", aspect="auto")
    ax.set_xticks(range(len(ct.columns))); ax.set_xticklabels(ct.columns)
    ax.set_yticks(range(len(ct.index)));   ax.set_yticklabels(ct.index)
    ax.set_xlabel("ALES class"); ax.set_ylabel("Surveyed crop (March 2025)")
    ax.set_title("Crop × ALES contingency\n(n = 993 survey points)")
    for i in range(arr.shape[0]):
        for j in range(arr.shape[1]):
            v = int(arr[i, j])
            if v > 0:
                ax.text(j, i, f"{v}", ha="center", va="center",
                        fontsize=7,
                        color="white" if arr_log[i, j] > arr_log.max() * 0.55 else "black")

    pr, b2, _ = read_masked(S7 / "policy_priority_0_100cm.tif")
    ext2 = (b2.left, b2.right, b2.bottom, b2.top)
    ax = axes[2]
    pri_cmap = ListedColormap(["#d73027", "#fdae61", "#1a9850"])
    pri_norm = BoundaryNorm([0.5, 1.5, 2.5, 3.5], pri_cmap.N)
    im = ax.imshow(pr, extent=ext2, origin="upper", cmap=pri_cmap, norm=pri_norm)
    ax.set_title("Policy priority\n(cropped, N1, top-quartile CF_Ks Δ)")
    ax.set_xticks([]); ax.set_yticks([])
    cbar = plt.colorbar(im, ax=ax, fraction=0.045, pad=0.02, ticks=[1, 2, 3])
    cbar.ax.set_yticklabels(["HIGH", "MEDIUM", "WATCH"])

    plt.savefig(OUT / "fig5_cropping_priority.png", bbox_inches="tight")
    plt.close()
    print("Wrote fig5_cropping_priority.png")


# ============================================================
# Figure 6 - Workflow schematic
# ============================================================
def fig6_workflow():
    fig, ax = plt.subplots(figsize=(11, 7.6))
    ax.set_xlim(0, 12); ax.set_ylim(0.4, 8.2); ax.set_axis_off()

    def box(x, y, w, h, text, color="#deebf7", ec="#2b6cb5"):
        b = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.15",
                           facecolor=color, edgecolor=ec, lw=1.2)
        ax.add_patch(b)
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center",
                fontsize=8.5, wrap=True)

    def arrow(x1, y1, x2, y2):
        ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2),
                                     arrowstyle="-|>", mutation_scale=12,
                                     color="#444", lw=1.0))

    # Inputs
    box(0.2, 6.4, 2.4, 1.0, "60 soil profiles\n(EC, ESP, CaCO₃,\nCEC, SOM, Ks…)", "#fff5b3", "#cc9900")
    box(0.2, 4.9, 2.4, 1.0, "20 water samples\n(IWQI inputs)", "#fff5b3", "#cc9900")
    box(0.2, 3.4, 2.4, 1.0, "53-band covariate\nstack 30 m\n(S1, S2, ERA5, DEM)", "#fff5b3", "#cc9900")
    # New input: independent crop-survey ground truth
    box(0.2, 1.7, 2.4, 1.1, "993 crop-survey points\n(March 2025, 15 crops)\nindependent ground truth",
        "#fff5b3", "#cc9900")

    # Stage 3 - prediction
    box(3.4, 5.9, 2.6, 1.5, "Spatial CV benchmark\n(QRF / RK / IDW / OK)\n→ q05 / q50 / q95\nPICP₉₀ ≈ 0.85", "#deebf7")
    # Stage 4
    box(3.4, 4.6, 2.6, 1.0, "Ordinary kriging\nIWQI → 30 m grid", "#deebf7")
    # Stage 5
    box(6.6, 6.0, 2.8, 1.4, "ALES-Arid index\n(Sys-Verheye)\n+ 200 MC draws/pixel", "#bbe3bb", "#2b8a3e")
    # Stage 6
    box(6.6, 4.3, 2.8, 1.3, "Counterfactual\ngap decomposition\n(7 scenarios × 100 MC)", "#bbe3bb", "#2b8a3e")
    # Stage 7
    box(6.6, 2.6, 2.8, 1.3, "Cropping confrontation\n(Sentinel-2 NDVI/MNDWI)\n→ priority raster", "#bbe3bb", "#2b8a3e")
    # Stage 7b - independent validation
    box(6.6, 0.9, 2.8, 1.3, "Independent validation\nvs. 993 crop points\n(User's acc. = 0.988)", "#bbe3bb", "#2b8a3e")

    # Outputs
    box(9.9, 5.9, 1.9, 1.2, "Probabilistic\nsuitability maps\n+ P(≥S2)", "#fbd1d1", "#b94343")
    box(9.9, 4.3, 1.9, 1.2, "Δ-index maps\nper limitation", "#fbd1d1", "#b94343")
    box(9.9, 2.6, 1.9, 1.2, "≈159 km²\nHIGH-priority zone", "#fbd1d1", "#b94343")
    box(9.9, 0.9, 1.9, 1.2, "Validated cropland\nmask + crop-resolved\nKPIs", "#fbd1d1", "#b94343")

    # Arrows
    arrow(2.6, 6.9, 3.4, 6.7)
    arrow(2.6, 5.4, 3.4, 5.1)
    arrow(2.6, 3.9, 3.4, 6.4)
    arrow(2.6, 3.9, 3.4, 4.9)
    arrow(2.6, 2.2, 6.6, 1.6)        # crop survey → validation block
    arrow(6.0, 6.6, 6.6, 6.7)
    arrow(6.0, 5.1, 6.6, 6.0)        # IWQI feeds ALES
    arrow(8.0, 6.0, 8.0, 5.6)        # ALES → CF
    arrow(8.0, 4.3, 8.0, 3.9)        # CF → cropping
    arrow(8.0, 2.6, 8.0, 2.2)        # cropping → validation
    arrow(9.4, 6.5, 9.9, 6.5)
    arrow(9.4, 4.9, 9.9, 4.9)
    arrow(9.4, 3.2, 9.9, 3.2)
    arrow(9.4, 1.5, 9.9, 1.5)

    # Section headers above each column of boxes (column-banner style).
    header_y = 7.7
    bar_y = 7.55
    headers = [
        (0.2, 2.4, "INPUTS",                       "#cc9900"),
        (3.4, 2.6, "SPATIAL PREDICTION",           "#2b6cb5"),
        (6.6, 2.8, "SUITABILITY + CAUSAL ANALYSIS", "#2b8a3e"),
        (9.9, 1.9, "OUTPUTS",                      "#b94343"),
    ]
    for x, w, label, col in headers:
        ax.add_patch(Rectangle((x, bar_y), w, 0.04, facecolor=col,
                               edgecolor="none", zorder=2))
        ax.text(x + w / 2, header_y, label, color=col,
                ha="center", va="bottom",
                fontsize=10, fontweight="bold")

    plt.savefig(OUT / "fig6_workflow.png", bbox_inches="tight")
    plt.savefig(OUT / "fig6_workflow.pdf", bbox_inches="tight")
    plt.close()
    print("Wrote fig6_workflow.{png,pdf}")


# ============================================================
# Figure 7 - Independent crop-survey validation
# ============================================================
def fig7_crop_validation():
    """Three-panel diagnostic: (a) survey points on ALES class map,
    (b) per-crop %N1+N2 + %HIGH-priority bar, (c) per-crop mean DKs bar."""
    s7b = ROOT / "outputs" / "stage7b"
    kpis = pd.read_csv(s7b / "crop_kpis.csv")
    pts  = pd.read_csv(s7b / "crop_points_extracted.csv")

    # Drop rare crops (n < 5) from the bar charts to keep the panels readable;
    # they remain in the contingency table inside Figure 5.
    kpi_main = kpis[kpis["n"] >= 5].copy()
    kpi_main = kpi_main.sort_values("n", ascending=False)

    fig = plt.figure(figsize=(14, 5.6), constrained_layout=True)
    gs  = fig.add_gridspec(1, 3, width_ratios=[1.05, 1.0, 1.0])

    # ---- (a) ALES class map + survey points ---------------------------------
    ax = fig.add_subplot(gs[0, 0])
    cls, b, _ = read_masked(S5 / "ales_class_q50_0_100cm.tif")
    ext = (b.left, b.right, b.bottom, b.top)
    ax.imshow(cls, extent=ext, origin="upper", cmap=ALES_CMAP, norm=ALES_NORM,
              alpha=0.85)
    # Plot dominant crops with distinct markers (top 5 by n)
    top5 = kpi_main.head(5)["crop"].tolist()
    palette = {"Wheat": "#1f77b4", "Clover": "#2ca02c",
               "Citrus": "#ff7f0e", "Garlic": "#9467bd",
               "Medical Plants": "#8c564b"}
    for crop in top5:
        sub = pts[pts["CROP"] == crop]
        ax.scatter(sub["x"], sub["y"], s=10, c=palette.get(crop, "#000000"),
                   edgecolor="white", linewidth=0.25, label=f"{crop} (n={len(sub)})",
                   zorder=4)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title("(a) Crop-survey points on ALES class map (0–100 cm)")
    crop_leg = ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.02),
              ncol=3, fontsize=7, frameon=False, handletextpad=0.4,
              columnspacing=0.9, title="Top-5 surveyed crops",
              title_fontsize=7.5)
    ax.add_artist(crop_leg)
    # Second legend: ALES suitability classes (background colours)
    ales_handles = [mpatches.Patch(facecolor=c, edgecolor="black", lw=0.3,
                                    label=l)
                    for c, l in zip(ALES_COLORS, ALES_LABELS)]
    ax.legend(handles=ales_handles, loc="upper left", fontsize=6.5,
              frameon=True, framealpha=0.85, borderpad=0.3,
              handlelength=1.0, handletextpad=0.3,
              title="ALES class", title_fontsize=7,
              ncol=5, columnspacing=0.6)
    add_scale_bar(ax, length_m=5000, label="5 km")

    # ---- (b) Per-crop %N1+N2 and %HIGH-priority -----------------------------
    ax = fig.add_subplot(gs[0, 1])
    x = np.arange(len(kpi_main))
    w = 0.4
    ax.bar(x - w/2, kpi_main["pct_N1_or_N2_%"], w,
           label="% in N1 + N2", color="#d73027")
    ax.bar(x + w/2, kpi_main["pct_HIGH_pri_%"], w,
           label="% in HIGH-priority", color="#fdae61")
    ax.set_xticks(x)
    ax.set_xticklabels([f"{c}\n(n={n})" for c, n in
                        zip(kpi_main["crop"], kpi_main["n"])],
                       rotation=30, ha="right", fontsize=7.5)
    ax.set_ylabel("% of survey points")
    ax.set_ylim(0, 105)
    ax.axhline(50, color="grey", lw=0.6, ls="--", alpha=0.5)
    ax.legend(frameon=False, fontsize=7.5)
    ax.grid(axis="y", alpha=0.3)
    ax.set_title("(b) Crop exposure to Ks-limited soils")

    # ---- (c) Per-crop mean CF_Ks gain ---------------------------------------
    ax = fig.add_subplot(gs[0, 2])
    ax.barh(np.arange(len(kpi_main))[::-1],
            kpi_main["mean_delta_Ks"], color="#3F7AA6")
    ax.set_yticks(np.arange(len(kpi_main))[::-1])
    ax.set_yticklabels([f"{c} (n={n})" for c, n in
                        zip(kpi_main["crop"], kpi_main["n"])],
                       fontsize=7.5)
    ax.set_xlabel("Mean Δ ALES q50 if Ks set to S2")
    ax.axvline(np.nanmean(kpi_main["mean_delta_Ks"]),
               color="red", ls="--", lw=0.7, alpha=0.7,
               label=f"survey mean = {np.nanmean(kpi_main['mean_delta_Ks']):.1f}")
    ax.legend(frameon=False, fontsize=7)
    ax.grid(axis="x", alpha=0.3)
    ax.set_title("(c) Crop-resolved CF_Ks suitability gain")

    plt.savefig(OUT / "fig7_crop_validation.png", bbox_inches="tight")
    plt.savefig(OUT / "fig7_crop_validation.pdf", bbox_inches="tight")
    plt.close()
    print("Wrote fig7_crop_validation.{png,pdf}")


def main():
    fig1_study_area()
    fig2_soil_uncertainty()
    fig3_baseline_ales()
    fig4_counterfactuals()
    fig5_cropping_priority()
    fig6_workflow()
    fig7_crop_validation()
    print(f"\nAll figures in {OUT}")


if __name__ == "__main__":
    main()
