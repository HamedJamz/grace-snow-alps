"""
Simple time-series tools for monthly mass/storage series.

The standard model for GRACE-type series is fitted by least squares:

    y(t) = a + b*(t - t0) + c1*cos(2*pi*t) + s1*sin(2*pi*t)
                          + c2*cos(4*pi*t) + s2*sin(4*pi*t)

t is in decimal years. a = offset, b = linear trend, (c1, s1) = annual and
(c2, s2) = semi-annual harmonics.
"""

import numpy as np
import pandas as pd

from .grace import decimal_year


def design_matrix(t, t0):
    """Columns: offset, trend, annual cos/sin, semi-annual cos/sin."""
    w = 2 * np.pi * t
    return np.column_stack([np.ones_like(t), t - t0, np.cos(w), np.sin(w),
                            np.cos(2 * w), np.sin(2 * w)])


def fit_seasonal_trend(series, t0=2010.0):
    """Least-squares fit of offset, trend, annual and semi-annual terms.

    Returns a dict with the fitted parameters, the annual amplitude and the month of
    the annual maximum, formal 1-sigma errors, and the fitted / residual series.
    Formal errors assume uncorrelated residuals and are optimistic; trend_err_ar1 is the
    trend error corrected for lag-1 autocorrelation (Santer et al., 2000).
    """
    s = series.dropna()
    t = np.asarray(decimal_year(s.index), dtype=float)
    A = design_matrix(t, t0)
    x, *_ = np.linalg.lstsq(A, s.values, rcond=None)
    fit = A @ x
    res = s.values - fit
    dof = len(s) - A.shape[1]
    sigma2 = res @ res / dof
    cov = sigma2 * np.linalg.inv(A.T @ A)
    err = np.sqrt(np.diag(cov))

    amp = np.hypot(x[2], x[3])
    # phase of the annual cosine: maximum at t = atan2(s1, c1) / (2 pi) years
    t_max = (np.arctan2(x[3], x[2]) / (2 * np.pi)) % 1.0
    amp_err = np.sqrt((x[2] * err[2]) ** 2 + (x[3] * err[3]) ** 2) / amp

    # Residuals of monthly storage series are autocorrelated, which makes the formal
    # errors too small. Following Santer et al. (2000), the trend error is inflated using
    # an effective sample size n_eff = n (1 - r1) / (1 + r1), with r1 the lag-1
    # autocorrelation of the residuals (only pairs of consecutive months are used).
    gap_days = np.diff(s.index.values).astype("timedelta64[D]").astype(float)
    consecutive = (gap_days > 20) & (gap_days < 45)
    r1 = np.corrcoef(res[:-1][consecutive], res[1:][consecutive])[0, 1]
    r1 = max(r1, 0.0)
    n_eff = len(s) * (1 - r1) / (1 + r1)
    p = A.shape[1]
    inflation = np.sqrt(dof / (n_eff - p)) if n_eff > p else np.nan
    return {
        "n": len(s),
        "offset": x[0], "trend_per_yr": x[1], "trend_err": err[1],
        "lag1_autocorr": r1, "n_eff": n_eff, "trend_err_ar1": err[1] * inflation,
        "annual_amp": amp, "annual_amp_err": amp_err,
        "annual_max_doy": t_max * 365.25,
        "semiannual_amp": np.hypot(x[4], x[5]),
        "residual_std": np.sqrt(sigma2),
        "fitted": pd.Series(fit, index=s.index),
        "residual": pd.Series(res, index=s.index),
    }


def monthly_climatology(series):
    """Mean value for each calendar month (1-12)."""
    s = series.dropna()
    return s.groupby(s.index.month).mean()


def to_monthly(series):
    """Put a series on a calendar-month index (YYYY-MM); duplicate months are averaged.

    GRACE solutions are centred on irregular dates, and two solutions occasionally fall
    into the same calendar month. Averaging them keeps one value per month.
    """
    s = series.dropna()
    return s.groupby(s.index.to_period("M")).mean()


def anomaly(series, start="2002-04", end="2020-03"):
    """Subtract the mean over a reference period (default: GravIS reference period)."""
    ref = series[(series.index >= pd.Period(start, "M")) & (series.index <= pd.Period(end, "M"))] \
        if isinstance(series.index, pd.PeriodIndex) else series.loc[start:end]
    return series - ref.mean()


def missing_months(series, start=None, end=None):
    """List calendar months without a value between the first and last month."""
    s = to_monthly(series)
    full = pd.period_range(start or s.index.min(), end or s.index.max(), freq="M")
    return [str(p) for p in full if p not in s.index]
