"""
ERA5-Land monthly snow water equivalent (SWE) and soil moisture.

Data: Copernicus Climate Change Service (2019): ERA5-Land monthly averaged data from 1950 to
present. C3S Climate Data Store. doi:10.24381/cds.68d2bb30. Model description: Muñoz-Sabater et al.
(2021), Earth Syst. Sci. Data 13, 4349-4383.

Variables used:
  sd     snow depth water equivalent       [m of water equivalent]
  swvl1-4 volumetric soil water, layers 0-7, 7-28, 28-100, 100-289 cm   [m3 m-3]
"""

from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
ERA5_FILE = "era5_land_monthly_alps_2002_2026.nc"

SOIL_LAYER_DEPTH_M = {"swvl1": 0.07, "swvl2": 0.21, "swvl3": 0.72, "swvl4": 1.89}

# ERA5-Land treats glacier points with a fixed, very large snow pack (about 10 m w.e.).
# These are not real seasonal snow and are removed before averaging.
GLACIER_SWE_THRESHOLD_MM = 2000.0


def load_era5(path=None):
    """Open the ERA5-Land file and return SWE and soil water storage in mm."""
    ds = xr.open_dataset(path or RAW / ERA5_FILE)
    ren = {k: v for k, v in {"valid_time": "time", "latitude": "lat", "longitude": "lon"}.items()
           if k in ds.dims or k in ds.coords}
    ds = ds.rename(ren).sortby("lat")
    out = xr.Dataset(coords={"time": ds.time, "lat": ds.lat, "lon": ds.lon})
    out["swe"] = ds["sd"] * 1000.0
    out["swe"].attrs["units"] = "mm"
    if all(k in ds for k in SOIL_LAYER_DEPTH_M):
        out["soil"] = sum(ds[k] * d for k, d in SOIL_LAYER_DEPTH_M.items()) * 1000.0
        out["soil"].attrs["units"] = "mm"
    return out


def glacier_mask(swe):
    """True for cells whose SWE never drops below the glacier threshold (permanent snow)."""
    return swe.min("time") > GLACIER_SWE_THRESHOLD_MM


def to_grace_grid(da, grace_lat, grace_lon):
    """Average 0.1 degree ERA5-Land cells into the 1 degree GRACE cells.

    A GRACE cell centred at (47.5, 10.5) covers 47-48 N and 10-11 E; every ERA5-Land
    point inside that square contributes equally (land points only, NaN over sea).
    """
    # small offset so that values like 47.99999 (float rounding) land in the right cell
    lat_c = np.floor(da.lat.values.astype(float) + 1e-6) + 0.5
    lon_c = np.floor(da.lon.values.astype(float) + 1e-6) + 0.5
    la, lo = np.meshgrid(lat_c, lon_c, indexing="ij")
    vals = da.transpose("time", "lat", "lon").values.reshape(da.sizes["time"], -1)
    df = pd.DataFrame(vals.T)
    df["lat"], df["lon"] = la.ravel(), lo.ravel()
    cell_mean = df.groupby(["lat", "lon"]).mean()          # NaN (sea, glacier) is skipped
    out = xr.DataArray(cell_mean.values.T, dims=("time", "cell"),
                       coords={"time": da.time, "cell": cell_mean.index}).unstack("cell")
    return out.reindex(lat=np.asarray(grace_lat, float), lon=np.asarray(grace_lon, float))
