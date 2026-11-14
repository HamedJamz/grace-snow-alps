"""
Read GravIS Level-3 terrestrial water storage (TWS) grids and build regional time series.

Data: GFZ GravIS RL06 (doi:10.5880/GFZ.GRAVIS_06_L3_TWS) and COST-G GravIS RL02
(doi:10.5880/COST-G.GRAVIS_02_L3_TWS), both CC BY 4.0. 1 degree grid, land only,
TWS anomalies in cm of equivalent water height relative to the mean 2002-04 to 2020-03.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"

PRODUCTS = {
    "GFZ RL06": "GRAVIS-3_GFZOP_0600_TWS_europe_subset.nc",
    "COST-G RL02": "GRAVIS-3_COSTG_0200_TWS_europe_subset.nc",
}

# Superconducting gravimeter site named in the SnowShifts project:
# Environmental Research Station Schneefernerhaus, Zugspitze (approx. 2650 m a.s.l.)
ZUGSPITZE = {"lat": 47.4167, "lon": 10.9796}


def load_tws(product="GFZ RL06"):
    """Open a GravIS TWS subset and convert TWS-related layers from cm to mm."""
    ds = xr.open_dataset(RAW / PRODUCTS[product])
    ds = ds.sortby("lat")                       # ascending latitude is easier to slice
    for v in ["tws", "std_tws", "leakage", "model_atmosphere"]:
        ds[v] = ds[v] * 10.0
        ds[v].attrs["units"] = "mm"
    return ds


def area_weights(lat):
    """Weights proportional to the area of a regular lat/lon cell (cos of latitude)."""
    return np.cos(np.deg2rad(lat))


def regional_mean(da, mask):
    """Area-weighted mean of a (time, lat, lon) field over cells where mask is True."""
    w = area_weights(da.lat) * mask
    return (da * w).sum(("lat", "lon")) / w.sum(("lat", "lon"))


def to_series(da):
    """Convert a 1-D DataArray over time to a pandas Series."""
    return pd.Series(da.values, index=pd.DatetimeIndex(da.time.values), name=da.name)


def decimal_year(index):
    """Decimal year for a DatetimeIndex (e.g. 2010-07-02 -> 2010.5)."""
    start = pd.to_datetime(index.year.astype(str) + "-01-01")
    end = pd.to_datetime((index.year + 1).astype(str) + "-01-01")
    return index.year + (index - start) / (end - start)
