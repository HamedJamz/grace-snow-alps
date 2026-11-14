"""
Cut a GravIS global TWS file to a Europe window.

    python data/subset_gravis.py GRAVIS-3_GFZOP_0600_TWS_GRID_GFZ_0006.nc data/raw/GRAVIS-3_GFZOP_0600_TWS_europe_subset.nc
"""
import sys

import xarray as xr

src, dst = sys.argv[1], sys.argv[2]
ds = xr.open_dataset(src)
ds = ds.assign_coords(lon=(((ds.lon + 180) % 360) - 180)).sortby("lon")
sub = ds.sel(lon=slice(-12, 32), lat=slice(34, 62) if ds.lat[0] < ds.lat[-1] else slice(62, 34))
sub.attrs.update(ds.attrs)
sub.to_netcdf(dst, encoding={v: {"zlib": True, "complevel": 4} for v in sub.data_vars})
print(dst, dict(sub.sizes))
