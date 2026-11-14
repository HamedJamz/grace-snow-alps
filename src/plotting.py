"""Figures for the GRACE / snow project (matplotlib, simple lon/lat maps)."""

from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
NE = ROOT / "data" / "raw" / "natural_earth"
FIG = ROOT / "figures"

C_TWS, C_TWS2, C_SWE, C_SOIL, C_REST = "#1f4e9c", "#6a9fd8", "#2a9d8f", "#c08a3e", "#888888"
SOURCE = ("Data: GravIS GFZ RL06 / COST-G RL02 TWS (GFZ, CC BY 4.0); ERA5-Land monthly "
          "(Copernicus C3S/ECMWF). Boundaries: Natural Earth.")

plt.rcParams.update({"font.size": 9, "axes.titlesize": 10.5, "axes.titleweight": "bold",
                     "axes.spines.top": False, "axes.spines.right": False,
                     "savefig.dpi": 300, "savefig.bbox": "tight"})


def save(fig, name):
    FIG.mkdir(exist_ok=True)
    fig.savefig(FIG / f"{name}.png", facecolor="white")
    plt.close(fig)


def _borders(ax, extent):
    c = gpd.read_file(NE / "ne_10m_admin_0_countries_europe.gpkg")
    lakes = gpd.read_file(NE / "ne_10m_lakes_europe.gpkg")
    c.boundary.plot(ax=ax, color="#333333", linewidth=0.5, zorder=3)
    lakes.plot(ax=ax, color="#d9e6f0", edgecolor="#9fb8cc", linewidth=0.2, zorder=2)
    ax.set_xlim(extent[0], extent[1]); ax.set_ylim(extent[2], extent[3])
    ax.set_aspect(1 / np.cos(np.deg2rad(np.mean(extent[2:]))))
    ax.set_xlabel("Longitude (°E)"); ax.set_ylabel("Latitude (°N)")


def _outline_cells(ax, mask, color="black", lw=1.2):
    """Draw the outer boundary of the selected 1° cells."""
    lat, lon = mask.lat.values, mask.lon.values
    m = mask.values
    for i in range(len(lat)):
        for j in range(len(lon)):
            if not m[i, j]:
                continue
            x0, x1, y0, y1 = lon[j] - .5, lon[j] + .5, lat[i] - .5, lat[i] + .5
            if i == 0 or not m[i - 1, j]: ax.plot([x0, x1], [y0, y0], color=color, lw=lw, zorder=4)
            if i == len(lat) - 1 or not m[i + 1, j]: ax.plot([x0, x1], [y1, y1], color=color, lw=lw, zorder=4)
            if j == 0 or not m[i, j - 1]: ax.plot([x0, x0], [y0, y1], color=color, lw=lw, zorder=4)
            if j == len(lon) - 1 or not m[i, j + 1]: ax.plot([x1, x1], [y0, y1], color=color, lw=lw, zorder=4)


def map_field(da, title, cbar_label, cmap="viridis", vmin=None, vmax=None, extent=(-10, 30, 36, 60),
              region_mask=None, point=None, name=None, grid_1deg=False):
    """Map a 2-D lat/lon field, optionally with the study-region outline and a site marker."""
    fig, ax = plt.subplots(figsize=(7.2, 5.6))
    pm = ax.pcolormesh(da.lon, da.lat, da, cmap=cmap, vmin=vmin, vmax=vmax, shading="nearest", zorder=1)
    _borders(ax, extent)
    if grid_1deg:
        for x in np.arange(np.floor(extent[0]), extent[1] + 1):
            ax.axvline(x, color="#7f8c99", lw=0.4, alpha=0.7, zorder=2.5)
        for y in np.arange(np.floor(extent[2]), extent[3] + 1):
            ax.axhline(y, color="#7f8c99", lw=0.4, alpha=0.7, zorder=2.5)
    if region_mask is not None:
        _outline_cells(ax, region_mask)
    if point is not None:
        ax.plot(point["lon"], point["lat"], marker="^", color="#c0392b", ms=8, mec="white", zorder=5)
        ax.annotate("Zugspitze", (point["lon"], point["lat"]), xytext=(7, -13),
                    textcoords="offset points", fontsize=8.5, fontweight="bold", color="#c0392b", zorder=5)
    cb = fig.colorbar(pm, ax=ax, shrink=0.8, pad=0.02)
    cb.set_label(cbar_label)
    ax.set_title(title, loc="left")
    fig.text(0.01, 0.01, SOURCE, fontsize=6.5, color="#555555")
    if name:
        save(fig, name)
    return fig


def plot_series(series_dict, title, ylabel="Storage anomaly (mm w.e.)", name=None, fits=None,
                gap=("2017-07", "2018-05")):
    """Plot several monthly series; optionally overlay least-squares fits."""
    fig, ax = plt.subplots(figsize=(10, 3.8))
    for label, (s, color) in series_dict.items():
        x = s.index.to_timestamp() if isinstance(s.index, pd.PeriodIndex) else s.index
        ax.plot(x, s.values, color=color, lw=1.2, marker=".", ms=3, label=label)
    if fits:
        for label, (f, color) in fits.items():
            ax.plot(f.index, f.values, color=color, lw=1, ls="--", alpha=0.8, label=label)
    if gap:
        ax.axvspan(pd.Timestamp(gap[0]), pd.Timestamp(gap[1]), color="#eeeeee", zorder=0)
        ax.text(pd.Timestamp("2017-12-01"), ax.get_ylim()[1] * 0.92, "GRACE /\nGRACE-FO\ngap",
                ha="center", va="top", fontsize=7, color="#666666")
    ax.axhline(0, color="#999999", lw=0.6)
    ax.set_ylabel(ylabel); ax.set_title(title, loc="left")
    ax.legend(frameon=False, ncol=4, fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.1))
    ax.grid(axis="y", color="#eeeeee")
    fig.text(0.01, -0.1, SOURCE, fontsize=6.5, color="#555555")
    if name:
        save(fig, name)
    return fig


def plot_climatology(clims, title, ylabel="Mean anomaly (mm w.e.)", name=None):
    """Mean seasonal cycle (calendar-month means) for several series."""
    fig, ax = plt.subplots(figsize=(6.8, 4))
    months = ["J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"]
    for label, (c, color, style) in clims.items():
        ax.plot(c.index, c.values, color=color, ls=style, lw=2, marker="o", ms=4, label=label)
    ax.set_xticks(range(1, 13)); ax.set_xticklabels(months)
    ax.axhline(0, color="#999999", lw=0.6)
    ax.set_ylabel(ylabel); ax.set_title(title, loc="left")
    ax.legend(frameon=False, fontsize=8)
    ax.grid(axis="y", color="#eeeeee")
    fig.text(0.01, -0.03, SOURCE, fontsize=6.5, color="#555555")
    if name:
        save(fig, name)
    return fig
