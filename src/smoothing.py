"""
Gaussian smoothing on the sphere, used as a rough stand-in for the spatial filtering of
GRACE data.

GRACE cannot resolve features smaller than a few hundred km, and the filters needed to
remove noise also spread ("leak") a concentrated signal such as mountain snow into the
surroundings. Smoothing a high-resolution field with a similar filter shows how much of a
regional signal survives at GRACE resolution. GravIS uses anisotropic VDK filters, so the
isotropic Gaussian here is only an approximation.
"""

import numpy as np
import xarray as xr

EARTH_RADIUS_KM = 6371.0


def gaussian_weights(lat, lon, radius_km):
    """Weight matrix (n x n) for a spherical Gaussian with half-weight radius radius_km.

    Jekeli-type kernel: w = exp(-b (1 - cos psi)), b = ln 2 / (1 - cos(r / R)),
    multiplied by cell area (cos latitude) and normalised per row.
    """
    la, lo = np.meshgrid(np.deg2rad(lat), np.deg2rad(lon), indexing="ij")
    la, lo = la.ravel(), lo.ravel()
    cos_psi = (np.sin(la)[:, None] * np.sin(la)[None, :] +
               np.cos(la)[:, None] * np.cos(la)[None, :] * np.cos(lo[:, None] - lo[None, :]))
    b = np.log(2) / (1 - np.cos(radius_km / EARTH_RADIUS_KM))
    w = np.exp(-b * (1 - np.clip(cos_psi, -1, 1))) * np.cos(la)[None, :]
    return w / w.sum(axis=1, keepdims=True)


def smooth(da, radius_km=300.0):
    """Smooth a (time, lat, lon) field; NaN is treated as zero signal (e.g. no snow)."""
    W = gaussian_weights(da.lat.values, da.lon.values, radius_km)
    vals = np.nan_to_num(da.transpose("time", "lat", "lon").values).reshape(da.sizes["time"], -1)
    out = (W @ vals.T).T.reshape(da.sizes["time"], da.sizes["lat"], da.sizes["lon"])
    return xr.DataArray(out, coords={"time": da.time, "lat": da.lat, "lon": da.lon},
                        dims=("time", "lat", "lon"), name=f"{da.name}_smoothed")
