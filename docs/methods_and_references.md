# Methods and the literature behind them

Each processing step and the reference it is based on. Full citations are in the README and in
`references.bib`.

| Step | What is done | Based on |
|---|---|---|
| Satellite gravimetry | Monthly gravity fields from GRACE / GRACE-FO | Tapley et al. (2004); Landerer et al. (2020); Tapley et al. (2019) |
| Gravity → water storage | Stokes coefficients to equivalent water height (thin-layer approximation) | Wahr et al. (1998) |
| Filtering | Anisotropic decorrelation filters (DDK / VDK) – done by GravIS | Swenson & Wahr (2006); Kusche (2007); Kusche et al. (2009); Horvath et al. (2018) |
| Background corrections | C20/C30 from SLR, degree-1, GIA (ICE-6G_D), AOD1B de-aliasing – done by GravIS | Dahle et al. (2025); Swenson et al. (2008); Peltier et al. (2018); Dobslaw et al. (2017) |
| Leakage | Not corrected; discussed and tested with smoothed model output | Landerer & Swenson (2012); Dobslaw et al. (2020) |
| Uncertainty of GravIS TWS | `std_tws` layer and spatial covariance model | Boergens et al. (2020, GEM) |
| Snow and soil moisture | ERA5-Land SWE and soil water layers | Muñoz-Sabater et al. (2021); C3S (2019) |
| Known limits of reanalysis snow | Biases in complex terrain | Monteiro & Morin (2023); Kouki et al. (2023) |
| Seasonal / trend model | Least squares: offset, trend, annual, semi-annual | Humphrey et al. (2016); Tapley et al. (2019) |
| Trend uncertainty | Effective sample size from lag-1 autocorrelation | Santer et al. (2000) |
| Spatial smoothing test | Isotropic Gaussian kernel on the sphere | Jekeli (1981); Wahr et al. (1998) |
| Context for the TWS trend | Alpine glacier mass loss; Central European droughts | Sommer et al. (2020); Boergens et al. (2020, GRL) |
| Local scale / Zugspitze | Superconducting gravimeter as hydrological sensor | Voigt et al. (2021) |

## Formulas used

**Seasonal and trend model** (t in decimal years, t₀ = 2010):

y(t) = a + b (t − t₀) + c₁ cos 2πt + s₁ sin 2πt + c₂ cos 4πt + s₂ sin 4πt

Annual amplitude A = √(c₁² + s₁²); time of the annual maximum t_max = atan2(s₁, c₁) / 2π.

**Trend error with autocorrelation** (Santer et al., 2000): with r₁ the lag-1 autocorrelation of
the residuals, n_eff = n (1 − r₁)/(1 + r₁) and σ_b,AR1 = σ_b,formal · √((n − p)/(n_eff − p)),
p = 6 fitted parameters.

**Gaussian kernel** (Jekeli, 1981): W(ψ) ∝ exp(−b (1 − cos ψ)), b = ln 2 / (1 − cos(r/R)),
with ψ the spherical distance, r the half-weight radius and R = 6371 km.

**Area weights:** w ∝ cos φ for a regular latitude–longitude grid.

**Mass from water height:** 1 mm over 1 km² = 10⁻⁶ Gt (1 Gt = 1 km³ of water).

**Soil water storage:** S = Σᵢ θᵢ dᵢ with layer thicknesses 0.07, 0.21, 0.72 and 1.89 m
(ERA5-Land; Muñoz-Sabater et al., 2021).
