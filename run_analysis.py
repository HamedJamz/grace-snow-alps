"""
Snow-related mass variations in the Alps from GRACE/GRACE-FO and ERA5-Land.

Steps
  1. read GravIS TWS (GFZ RL06 and COST-G RL02) and ERA5-Land SWE / soil moisture
  2. put ERA5-Land on the 1 degree GRACE grid and define the Alpine study region
  3. regional time series, seasonal/trend fit, mean seasonal cycle
  4. how much of the seasonal TWS signal is snow (and soil moisture)?
  5. Zugspitze: local snow signal vs. what GRACE can resolve
  6. figures, tables and a few sanity checks

Usage:  python run_analysis.py
"""

import json
import warnings

import numpy as np
import pandas as pd
import xarray as xr

from src.grace import PRODUCTS, ZUGSPITZE, ROOT, load_tws, regional_mean, to_series, area_weights
from src.era5 import load_era5, glacier_mask, to_grace_grid
from src.timeseries import (fit_seasonal_trend, monthly_climatology, to_monthly, anomaly,
                            missing_months, design_matrix)
from src.smoothing import smooth
from src import plotting as P

warnings.filterwarnings("ignore", category=RuntimeWarning)
OUT = ROOT / "outputs"
OUT.mkdir(exist_ok=True)
REF = ("2002-04", "2020-03")          # GravIS reference period, used for all anomalies
SWE_REGION_THRESHOLD_MM = 50.0        # mean seasonal SWE maximum needed to call a 1° cell "snowy"
ALPS_BOX = dict(lon=(5.0, 16.5), lat=(43.5, 48.5))
checks = []


def check(name, ok, detail=""):
    checks.append({"check": name, "passed": bool(ok), "detail": detail})
    print(("PASS " if ok else "FAIL ") + name + (f"  ({detail})" if detail else ""))


def monthly_series(da):
    """Regional DataArray (time) -> monthly pandas Series with a PeriodIndex."""
    return to_monthly(to_series(da))


# ------------------------------------------------------------------------------------
# 1. data
# ------------------------------------------------------------------------------------
grace = {p: load_tws(p) for p in PRODUCTS}
era = load_era5()
check("ERA5-Land SWE is non-negative (up to float rounding)", float(era.swe.min()) > -1e-6)

glac = glacier_mask(era.swe)
swe = era.swe.where(~glac)
soil = era.soil.where(~glac)
print(f"glacier points removed: {int(glac.sum())} of {int(era.swe.isel(time=0).notnull().sum())} land points")

g0 = grace["GFZ RL06"]
swe_1deg = to_grace_grid(swe, g0.lat, g0.lon)
soil_1deg = to_grace_grid(soil, g0.lat, g0.lon)

# ------------------------------------------------------------------------------------
# 2. study region: snowy 1° cells in the Alpine box with complete GRACE data
# ------------------------------------------------------------------------------------
yearly_max = swe_1deg.groupby("time.year").max("time")
mean_max = yearly_max.where(yearly_max.year < 2026).mean("year")
in_box = ((g0.lon >= ALPS_BOX["lon"][0]) & (g0.lon <= ALPS_BOX["lon"][1]) &
          (g0.lat >= ALPS_BOX["lat"][0]) & (g0.lat <= ALPS_BOX["lat"][1]))
grace_ok = g0.tws.notnull().all("time") & grace["COST-G RL02"].tws.notnull().all("time")
region = (in_box & grace_ok & (mean_max >= SWE_REGION_THRESHOLD_MM)).fillna(False).transpose("lat", "lon")
n_cells = int(region.sum())
area_km2 = float((area_weights(g0.lat) * region).sum() * 111.32 ** 2)
print(f"study region: {n_cells} cells, about {area_km2:,.0f} km2")
check("Study region has at least 10 cells", n_cells >= 10, f"{n_cells} cells")

# ------------------------------------------------------------------------------------
# 3. regional series
# ------------------------------------------------------------------------------------
tws = {p: monthly_series(regional_mean(grace[p].tws, region)) for p in PRODUCTS}
leak = monthly_series(regional_mean(g0.leakage, region))
std_tws = monthly_series(regional_mean(g0.std_tws, region))
swe_r = anomaly(monthly_series(regional_mean(swe_1deg, region)), *REF)
soil_r = anomaly(monthly_series(regional_mean(soil_1deg, region)), *REF)

# SWE smoothed with a 300 km Gaussian: how much of the snow signal would survive GRACE-like filtering?
SMOOTH_RADIUS_KM = 300.0
swe_smooth = smooth(swe_1deg.rename("swe"), SMOOTH_RADIUS_KM)
swe_s_r = anomaly(monthly_series(regional_mean(swe_smooth, region)), *REF)

# compare on GRACE months only
months = tws["GFZ RL06"].index.intersection(swe_r.index)
df = pd.DataFrame({"tws_gfz": tws["GFZ RL06"], "tws_costg": tws["COST-G RL02"],
                   "swe": swe_r, "swe_smoothed_300km": swe_s_r, "soil": soil_r}).loc[months]
df["rest"] = df["tws_gfz"] - df["swe"] - df["soil"]
df.index.name = "month"
df.round(2).to_csv(OUT / "alps_monthly_storage_anomalies_mm.csv")

missing = missing_months(tws["GFZ RL06"].to_timestamp())
check("GRACE months handled: duplicates averaged, gaps kept as gaps",
      tws["GFZ RL06"].index.is_unique, f"{len(missing)} missing months incl. 2017-07..2018-05 gap")

# seasonal / trend fits (on the GRACE months, time stamp = middle of month)
def fit(s):
    s2 = s.copy(); s2.index = s2.index.to_timestamp() + pd.Timedelta(days=14)
    return fit_seasonal_trend(s2)

fits = {k: fit(df[k]) for k in ["tws_gfz", "tws_costg", "swe", "swe_smoothed_300km", "soil", "rest"]}
fit_tab = pd.DataFrame({k: {kk: vv for kk, vv in v.items() if not isinstance(vv, pd.Series)}
                        for k, v in fits.items()}).T
fit_tab.round(2).to_csv(OUT / "alps_seasonal_trend_fit.csv")
print(fit_tab.round(1))

clim = {k: monthly_climatology(df[k].set_axis(df.index.to_timestamp())) for k in df.columns}
pd.DataFrame(clim).round(1).to_csv(OUT / "alps_mean_seasonal_cycle_mm.csv")

# share of the seasonal TWS signal explained by snow
amp = {k: fits[k]["annual_amp"] for k in fits}
snow_share = amp["swe"] / amp["tws_gfz"]
r_swe = np.corrcoef(df["tws_gfz"], df["swe"])[0, 1]
r_swe_soil = np.corrcoef(df["tws_gfz"], df["swe"] + df["soil"])[0, 1]
both = df[["tws_gfz", "tws_costg"]].dropna()
r_products = np.corrcoef(both["tws_gfz"], both["tws_costg"])[0, 1]
rms_products = float(np.sqrt(((both["tws_gfz"] - both["tws_costg"]) ** 2).mean()))
print(f"annual amplitude: TWS {amp['tws_gfz']:.0f} mm, SWE {amp['swe']:.0f} mm, "
      f"soil {amp['soil']:.0f} mm -> SWE/TWS = {snow_share:.2f}")

# ------------------------------------------------------------------------------------
# 4. Zugspitze: point vs. 1° cell vs. GRACE
# ------------------------------------------------------------------------------------
pt_all = era.swe.sel(lat=ZUGSPITZE["lat"], lon=ZUGSPITZE["lon"], method="nearest")
pt_is_glacier = bool(glac.sel(lat=pt_all.lat, lon=pt_all.lon))
# if the nearest ERA5-Land point is a glacier point, use the nearest non-glacier land point
cand = swe.sel(lat=slice(47.2, 47.7), lon=slice(10.7, 11.3))
ok = cand.notnull().all("time")
dist = np.hypot((cand.lat - ZUGSPITZE["lat"]), (cand.lon - ZUGSPITZE["lon"]) * np.cos(np.deg2rad(47.4)))
dist = dist.where(ok)
i_lat, i_lon = np.unravel_index(int(np.nanargmin(dist.transpose("lat", "lon").values)), dist.transpose("lat", "lon").shape)
pt = cand.isel(lat=i_lat, lon=i_lon)
cell = dict(lat=np.floor(ZUGSPITZE["lat"]) + 0.5, lon=np.floor(ZUGSPITZE["lon"]) + 0.5)

z = pd.DataFrame({
    "swe_point": anomaly(to_monthly(to_series(pt)), *REF),
    "swe_1deg_cell": anomaly(monthly_series(swe_1deg.sel(**cell)), *REF),
    "tws_grace_1deg_cell": monthly_series(g0.tws.sel(**cell)),
}).loc[months]
zclim = {k: monthly_climatology(z[k].set_axis(z.index.to_timestamp())) for k in z.columns}
zamp = {k: fit(z[k])["annual_amp"] for k in z.columns}
pd.DataFrame(zclim).round(1).to_csv(OUT / "zugspitze_mean_seasonal_cycle_mm.csv")
print("Zugspitze annual amplitudes (mm):", {k: round(v) for k, v in zamp.items()})

# ------------------------------------------------------------------------------------
# 5. maps
# ------------------------------------------------------------------------------------
def amplitude_map(da):
    """Annual amplitude per cell (same harmonic model, cells with complete series)."""
    s = da.dropna("time", how="all")
    s = s.groupby(s.time.dt.strftime("%Y-%m")).mean("time")
    t = pd.PeriodIndex(s.strftime.values, freq="M").to_timestamp() + pd.Timedelta(days=14)
    from src.grace import decimal_year
    tt = np.asarray(decimal_year(t), float)
    A = design_matrix(tt, 2010.0)
    Y = s.values.reshape(len(tt), -1)
    good = np.isfinite(Y).all(0)
    X = np.full((A.shape[1], Y.shape[1]), np.nan)
    X[:, good] = np.linalg.lstsq(A, Y[:, good], rcond=None)[0]
    amp_ = np.hypot(X[2], X[3]).reshape(s.shape[1:])
    tr = X[1].reshape(s.shape[1:])
    return (xr.DataArray(amp_, coords={"lat": s.lat, "lon": s.lon}, dims=("lat", "lon")),
            xr.DataArray(tr, coords={"lat": s.lat, "lon": s.lon}, dims=("lat", "lon")))

amp_map, trend_map = amplitude_map(g0.tws)
P.map_field(amp_map, "Annual amplitude of GRACE TWS (GravIS GFZ RL06, 2002–2026)",
            "Annual amplitude (mm water equivalent)", cmap="viridis", vmin=0, vmax=120,
            region_mask=region, point=ZUGSPITZE, name="fig1_grace_tws_annual_amplitude")

march = swe.where(swe.time.dt.month == 3, drop=True).mean("time")
P.map_field(march, "ERA5-Land mean March SWE (0.1°) and the 1° GRACE grid",
            "SWE (mm water equivalent)", cmap="Blues", vmin=0, vmax=600, extent=(4, 18, 43, 49),
            region_mask=region, point=ZUGSPITZE, grid_1deg=True, name="fig2_era5_march_swe_alps")

# reindex to every calendar month so that lines break at missing GRACE months
dfp = df.reindex(pd.period_range(df.index.min(), df.index.max(), freq="M"))
P.plot_series({"TWS – GravIS GFZ RL06": (dfp["tws_gfz"], P.C_TWS),
               "TWS – GravIS COST-G RL02": (dfp["tws_costg"], P.C_TWS2),
               "SWE – ERA5-Land": (dfp["swe"], P.C_SWE),
               "Soil moisture (0–2.89 m) – ERA5-Land": (dfp["soil"], P.C_SOIL)},
              f"Alpine snow region ({n_cells} cells × 1°): monthly storage anomalies",
              name="fig3_alps_time_series")

P.plot_climatology({"TWS (GRACE, GFZ)": (clim["tws_gfz"], P.C_TWS, "-"),
                    "SWE (ERA5-Land)": (clim["swe"], P.C_SWE, "-"),
                    "SWE smoothed (Gaussian 300 km)": (clim["swe_smoothed_300km"], P.C_SWE, ":"),
                    "Soil moisture (ERA5-Land)": (clim["soil"], P.C_SOIL, "-"),
                    "TWS − SWE − soil": (clim["rest"], P.C_REST, "--")},
                   "Mean seasonal cycle in the Alpine snow region", name="fig4_alps_seasonal_cycle")

P.plot_climatology({"SWE at ERA5-Land point near Zugspitze": (zclim["swe_point"], "#c0392b", "-"),
                    "SWE averaged over the 1° cell": (zclim["swe_1deg_cell"], P.C_SWE, "-"),
                    "TWS from GRACE for the same 1° cell": (zclim["tws_grace_1deg_cell"], P.C_TWS, "-")},
                   "Zugspitze: local snow signal vs. the 1° scale", name="fig5_zugspitze_scale")

# ------------------------------------------------------------------------------------
# 5b. sensitivity tests: region threshold and smoothing radius
# ------------------------------------------------------------------------------------
def fit_region(field, mask, is_anomaly=False):
    """Seasonal/trend fit of the regional mean of field over mask (GRACE months only)."""
    s = monthly_series(regional_mean(field, mask))
    if not is_anomaly:
        s = anomaly(s, *REF)
    return fit(s.loc[s.index.intersection(months)])

def mm_to_gt(mm, km2):
    """1 mm of water over 1 km2 is 1e-6 Gt (1 Gt = 1 km3 of water)."""
    return mm * km2 * 1e-6

rows = []
for thr in [30, 50, 100]:
    m = (in_box & grace_ok & (mean_max >= thr)).fillna(False).transpose("lat", "lon")
    a_km2 = float((area_weights(g0.lat) * m).sum() * 111.32 ** 2)
    f_g = fit_region(g0.tws, m, True)
    f_c = fit_region(grace["COST-G RL02"].tws, m, True)
    f_s = fit_region(swe_1deg, m)
    rows.append({"swe_threshold_mm": thr, "cells": int(m.sum()), "area_km2": round(a_km2),
                 "tws_gfz_amp_mm": f_g["annual_amp"], "tws_costg_amp_mm": f_c["annual_amp"],
                 "swe_amp_mm": f_s["annual_amp"], "swe_to_tws_ratio": f_s["annual_amp"] / f_g["annual_amp"],
                 "tws_gfz_trend_mm_yr": f_g["trend_per_yr"], "tws_gfz_trend_err_ar1": f_g["trend_err_ar1"],
                 "tws_gfz_trend_gt_yr": mm_to_gt(f_g["trend_per_yr"], a_km2)})
sens_region = pd.DataFrame(rows).round(2)
sens_region.to_csv(OUT / "sensitivity_region_threshold.csv", index=False)
print(sens_region.to_string(index=False))

rows = []
for radius in [150, 200, 300, 400]:
    f = fit_region(smooth(swe_1deg.rename("swe"), radius), region)
    rows.append({"gaussian_radius_km": radius, "smoothed_swe_amp_mm": f["annual_amp"],
                 "ratio_to_grace_tws": f["annual_amp"] / amp["tws_gfz"]})
sens_radius = pd.DataFrame(rows).round(2)
sens_radius.to_csv(OUT / "sensitivity_smoothing_radius.csv", index=False)
print(sens_radius.to_string(index=False))

# ------------------------------------------------------------------------------------
# 6. checks and summary
# ------------------------------------------------------------------------------------
check("Two GRACE products agree (r > 0.9)", r_products > 0.9, f"r = {r_products:.3f}, RMS diff = {rms_products:.1f} mm")
check("SWE peak is in late winter / spring (Feb–Apr)", clim["swe"].idxmax() in (2, 3, 4),
      f"maximum in month {clim['swe'].idxmax()}")
src_mean = float(swe.sel(lat=slice(46.99, 47.95), lon=slice(9.99, 10.95)).isel(time=2).mean())
check("1° SWE cell value equals the mean of its 0.1° source points",
      abs(src_mean - float(swe_1deg.sel(lat=47.5, lon=10.5).isel(time=2))) < 1e-3,
      f"{src_mean:.2f} mm")
check("Harmonic fit reproduces a synthetic annual signal",
      abs(fit_seasonal_trend(pd.Series(100 * np.cos(2 * np.pi * (np.arange(240) / 12 + 0.04)),
          index=pd.date_range("2005-01-15", periods=240, freq="MS") + pd.Timedelta(days=14)))["annual_amp"] - 100) < 1)

summary = {
    "region_cells": n_cells, "region_area_km2": round(area_km2),
    "months_compared": len(df), "first_month": str(df.index[0]), "last_month": str(df.index[-1]),
    "annual_amplitude_mm": {k: round(float(v), 1) for k, v in amp.items()},
    "annual_max_day_of_year": {k: round(float(fits[k]["annual_max_doy"])) for k in fits},
    "trend_mm_per_yr": {k: {"trend": round(float(fits[k]["trend_per_yr"]), 2),
                            "err_formal": round(float(fits[k]["trend_err"]), 2),
                            "err_ar1": round(float(fits[k]["trend_err_ar1"]), 2),
                            "lag1_autocorr": round(float(fits[k]["lag1_autocorr"]), 2)} for k in fits},
    "tws_gfz_trend_gt_per_yr": round(mm_to_gt(float(fits["tws_gfz"]["trend_per_yr"]), area_km2), 2),
    "tws_gfz_trend_err_ar1_gt_per_yr": round(mm_to_gt(float(fits["tws_gfz"]["trend_err_ar1"]), area_km2), 2),
    "swe_to_tws_annual_amplitude_ratio": round(float(snow_share), 2),
    "smoothed_swe_to_tws_annual_amplitude_ratio": round(float(amp["swe_smoothed_300km"] / amp["tws_gfz"]), 2),
    "smoothing_radius_km": SMOOTH_RADIUS_KM,
    "corr_tws_swe": round(float(r_swe), 3), "corr_tws_swe_plus_soil": round(float(r_swe_soil), 3),
    "corr_gfz_costg": round(float(r_products), 3), "rms_gfz_minus_costg_mm": round(rms_products, 1),
    "mean_regional_std_tws_mm": round(float(std_tws.mean()), 1),
    "rms_regional_leakage_mm": round(float(np.sqrt((leak ** 2).mean())), 1),
    "zugspitze": {"era5_point_used": [float(pt.lat), float(pt.lon)], "nearest_point_was_glacier": pt_is_glacier,
                  "annual_amplitude_mm": {k: round(float(v), 1) for k, v in zamp.items()}},
    "glacier_points_removed": int(glac.sum()),
    "checks": checks,
}
(OUT / "results_summary.json").write_text(json.dumps(summary, indent=2, default=str))
print(json.dumps({k: v for k, v in summary.items() if k != "checks"}, indent=1, default=str))
print(f"{sum(c['passed'] for c in checks)}/{len(checks)} checks passed")
