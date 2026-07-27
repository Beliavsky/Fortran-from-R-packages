# yieldcurves

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/yieldcurves)](https://CRAN.R-project.org/package=yieldcurves)
[![CRAN downloads](https://cranlogs.r-pkg.org/badges/yieldcurves)](https://CRAN.R-project.org/package=yieldcurves)
[![Total Downloads](https://cranlogs.r-pkg.org/badges/grand-total/yieldcurves)](https://CRAN.R-project.org/package=yieldcurves)
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

A technical working paper for this package can be found [here](https://charlescoverdale.github.io/files/coverdale_yieldcurves_2026.pdf).

**An R package for working with yield curves: fit smooth curves to bond yields, extract forward rates and discount factors, compute duration and convexity, and decompose curve movements.**

## Installation

```r
install.packages("yieldcurves")

# Or install the development version from GitHub
# install.packages("devtools")
devtools::install_github("charlescoverdale/yieldcurves")
```

```r
library(yieldcurves)

# US Treasury yields: maturities in years, rates as decimals (5% = 0.05)
maturities <- c(0.25, 0.5, 1, 2, 5, 10, 30)
rates <- c(0.052, 0.050, 0.048, 0.045, 0.042, 0.040, 0.043)

fit <- yc_nelson_siegel(maturities, rates)
fit
#> -- Yield Curve (Nelson-Siegel) --
#> * Type: "zero"
#> * Maturities: 7 (0.25Y to 30Y)
#> * Rate range: 4% to 5.2%
#> * RMSE: 8.6 bps
#> * Parameters: beta0=0.04127, beta1=0.01347, beta2=-0.0091, tau=1

plot(fit)
```


## Why yieldcurves?

Government bond yields are published daily at a handful of maturities: 1-month, 3-month, 1-year, 2-year, 5-year, 10-year, 30-year. But most questions in fixed income require yields at maturities that are not published, or require derived quantities (forward rates, discount factors, risk measures) that cannot be read off a table. If you need a discount factor at 7.5 years, or the implied rate for a 5-year loan starting in 5 years, you need a fitted model.

The Nelson-Siegel (1987) and Svensson (1994) models are the standard solution. They fit the entire curve with 4 to 6 parameters that have economic meaning: level (long-run rate), slope (term premium), and curvature (medium-term humps). Over 20 central banks use these models for their official yield curve estimates. Once you have a fitted curve, you can extract forward rates, discount factors, carry and roll-down, duration and convexity, and Z-spreads.

There are two existing R packages in this space. `YieldCurve` (last updated 2022) fits Nelson-Siegel and Svensson models but requires xts/zoo objects as input, does not compute duration, convexity, Z-spreads, or PCA, and has no carry/roll-down analysis. `termstrc` was more comprehensive (cubic splines, forward rates) but was archived from CRAN in 2018 due to failing checks and is no longer maintained. Neither package covers the full workflow from curve fitting through to risk measures and portfolio analytics. `yieldcurves` fills this gap: it works with plain numeric vectors, has no heavy dependencies, and covers fitting, extraction, risk measures, and decomposition in one package.


## How does this compare to existing packages?

| Feature | yieldcurves | YieldCurve | termstrc |
|---------|-------------|------------|----------|
| Nelson-Siegel fitting | Yes | Yes | Yes |
| Svensson fitting | Yes | Yes | Yes |
| Cubic spline | Yes | No | Yes |
| Weighted fitting | Yes | No | No |
| Forward rates (analytical) | Yes | No | Yes |
| Discount factors | Yes | No | Yes |
| Duration and convexity | Yes | No | No |
| Z-spread | Yes | No | No |
| Key rate durations | Yes | No | No |
| Par/zero conversions | Yes | No | Yes |
| PCA decomposition | Yes | No | No |
| Carry and roll-down | Yes | No | No |
| Slope measures | Yes | No | No |
| Works with plain vectors | Yes | No (needs xts/zoo) | No |
| Last updated | 2026 | 2022 | 2015 |


## Examples

### Fit a yield curve

The Nelson-Siegel model takes two inputs: maturities (in years) and rates (as decimals, so 5% = 0.05). It returns a fitted curve you can query at any maturity.

```r
library(yieldcurves)

maturities <- c(0.25, 0.5, 1, 2, 3, 5, 7, 10, 20, 30)
rates <- c(0.052, 0.050, 0.048, 0.045, 0.043, 0.042, 0.041,
           0.040, 0.042, 0.043)

fit <- yc_nelson_siegel(maturities, rates)
fit
#> -- Yield Curve (Nelson-Siegel) --
#> * Type: "zero"
#> * Maturities: 10 (0.25Y to 30Y)
#> * Rate range: 4% to 5.2%
#> * RMSE: 8.6 bps
#> * Parameters: beta0=0.04127, beta1=0.01347, beta2=-0.0091, tau=1

plot(fit)
```

### Extract forward rates and discount factors

Forward rates tell you what the market implies about future interest rates. Discount factors tell you what a future cash flow is worth today.

```r
# Forward rates at 1, 5, and 10 years
yc_forward(fit, maturities = c(1, 5, 10))
#>   maturity forward_rate
#> 1        1       0.0468
#> 2        5       0.0393
#> 3       10       0.0413

# Discount factors
yc_discount(fit, maturities = c(1, 5, 10, 30))
#>   maturity discount_factor
#> 1        1          0.9530
#> 2        5          0.8103
#> 3       10          0.6619
#> 4       30          0.2817   # $1 in 30 years is worth $0.28 today
```

### Compute carry and roll-down

How much do you earn from holding a bond, assuming the curve does not move? Carry is the yield income minus the funding cost. Roll-down is the capital gain as the bond's remaining maturity shortens and it slides to a lower-rate part of the curve.

```r
yc_carry(fit, maturities = c(2, 5, 10, 30))
#>   maturity     carry   rolldown      total
#> 1        2  0.000170  0.005143  0.005313   # 53 bps total on the 2Y
#> 2        5  0.000013  0.004478  0.004491   # 45 bps on the 5Y
#> 3       10 -0.000059  0.003260  0.003201   # 32 bps on the 10Y
#> 4       30  0.000033 -0.001507 -0.001474   # Negative on the 30Y
```

### Compute duration and convexity for a coupon bond

Duration measures how sensitive a bond's price is to rate changes. A modified duration of 7.87 means a 1% rate rise causes roughly a 7.87% price drop.

```r
# 10-year bond with 5% coupon at 4.5% yield, semi-annual coupons
yc_bond_duration(face = 100, coupon_rate = 0.05, maturity = 10,
                 yield = 0.045, frequency = 2)
#> $macaulay_duration
#> [1] 8.05
#>
#> $modified_duration
#> [1] 7.87
#>
#> $convexity
#> [1] 73.7
#>
#> $price
#> [1] 104.01
```

### Compute the Z-spread on a bond

The Z-spread is the constant spread over the benchmark zero curve that reprices a bond to its market price. A positive Z-spread means the bond yields more than the benchmark.

```r
# Benchmark zero curve
curve <- yc_curve(c(0.5, 1, 2, 5, 10), c(0.03, 0.035, 0.04, 0.042, 0.045))

# Bond trading below par
result <- yc_zspread(price = 95, coupon_rate = 0.04, maturity = 5,
                     curve = curve, frequency = 2)
result$zspread
#> [1] 0.0148   # 148 bps over the benchmark curve
```

### Decompose yield curve movements with PCA

Principal component analysis decomposes a time series of yield curves into orthogonal factors. Litterman and Scheinkman (1991) showed that three factors (level, slope, curvature) explain over 95% of yield curve movements.

```r
# Simulate a time series of yield curves (200 days, 5 tenors)
set.seed(42)
n_days <- 200
tenors <- c(1, 2, 5, 10, 30)
base <- c(0.045, 0.043, 0.042, 0.040, 0.043)
level <- cumsum(rnorm(n_days, 0, 0.001))
curves <- matrix(NA, n_days, length(tenors))
for (i in seq_len(n_days)) curves[i, ] <- base + level[i]
colnames(curves) <- paste0(tenors, "Y")

pca <- yc_pca(curves)
pca
#> -- Yield Curve PCA --
#> * Components: 3
#> * PC1 (Level): 97.2% variance
#> * PC2 (Slope): 1.2% variance
#> * PC3 (Curvature): 0.8% variance

plot(pca)
```


## Key terms

- **Zero rates** (spot rates): yields on zero-coupon bonds. These are the building blocks of the yield curve.
- **Par rates**: the coupon rate at which a bond prices at par (100). This is what the news means by "the 10-year yield."
- **Forward rates**: rates implied by the curve for a future period. The "5y5y forward" is the rate the market implies for a 5-year loan starting in 5 years.
- **Basis points (bps)**: 1 bp = 0.01%. All rates in this package are decimals: 5% = `0.05`, 50 bps = `0.005`.


## What data do I need?

You supply two numeric vectors: **maturities** (in years) and **rates** (as decimals, so 5% = 0.05). You can type them in directly or pull them from a data source.

### Option 1: Type rates in directly

```r
library(yieldcurves)

maturities <- c(0.25, 0.5, 1, 2, 3, 5, 7, 10, 20, 30)
rates <- c(0.052, 0.050, 0.048, 0.045, 0.043, 0.042, 0.041,
           0.040, 0.042, 0.043)

fit <- yc_nelson_siegel(maturities, rates)
plot(fit)
```

### Option 2: Pull US Treasury yields from FRED

```r
# 1. Install the fred package (one time)
install.packages("fred")

# 2. Get a free API key from https://fred.stlouisfed.org/docs/api/api_key.html
fred::fred_set_key("your_api_key_here")

# 3. Download Treasury constant maturity rates
library(fred)
treasury <- fred_series(c("DGS1", "DGS2", "DGS5", "DGS10", "DGS30"))

# 4. Extract the most recent observation
latest <- treasury[nrow(treasury), ]
maturities <- c(1, 2, 5, 10, 30)
rates <- as.numeric(latest[, -1]) / 100  # FRED reports percent, convert to decimal

# 5. Fit a yield curve
library(yieldcurves)
fit <- yc_nelson_siegel(maturities, rates)
plot(fit)
```

### Where to find yield data

| Country | Source | R package | Series / function |
|---------|--------|-----------|-------------------|
| US | US Treasury / FRED | [fred](https://cran.r-project.org/package=fred) | DGS1MO, DGS3MO, DGS1, DGS2, DGS5, DGS10, DGS30 |
| UK | Bank of England | [boe](https://cran.r-project.org/package=boe) | `boe_yield_curve()` |
| Euro area | European Central Bank | [readecb](https://cran.r-project.org/package=readecb) | `ecb_yield_curve()` |
| Japan | Ministry of Finance | Download CSV | mof.go.jp |
| Australia | RBA | Download CSV | rba.gov.au |
| Canada | Bank of Canada | [fred](https://cran.r-project.org/package=fred) | DGS series on FRED |


## Functions

### Curve construction and fitting

| Function | Description |
|----------|-------------|
| `yc_curve()` | Create a yield curve object from maturity-rate pairs |
| `yc_nelson_siegel()` | Fit a Nelson-Siegel (1987) model |
| `yc_svensson()` | Fit a Svensson (1994) model |
| `yc_cubic_spline()` | Fit a cubic spline |
| `yc_fit()` | Unified fitting interface (dispatches to the above) |

### Evaluation and extraction

| Function | Description |
|----------|-------------|
| `yc_predict()` | Evaluate a fitted curve at new maturities |
| `yc_forward()` | Extract instantaneous or forward-forward rates |
| `yc_discount()` | Compute discount factors (continuous, annual, or semi-annual) |
| `yc_interpolate()` | Interpolate between observed rates (linear, log-linear, or cubic) |

### Rate conversions

| Function | Description |
|----------|-------------|
| `yc_par_to_zero()` | Convert par rates to zero rates via bootstrap stripping |
| `yc_zero_to_par()` | Convert zero rates to par rates |

### Risk measures

| Function | Description |
|----------|-------------|
| `yc_duration()` | Macaulay/modified duration and convexity for zero-coupon bonds |
| `yc_bond_duration()` | Duration and convexity for coupon-bearing bonds |
| `yc_zspread()` | Z-spread (zero-volatility spread) computation |
| `yc_key_rate_duration()` | Key rate durations with triangular bump profiles |

### Analysis and decomposition

| Function | Description |
|----------|-------------|
| `yc_carry()` | Carry and roll-down decomposition |
| `yc_pca()` | Principal component analysis of yield curve time series |
| `yc_slope()` | Spread measures (2s10s, 2s30s, butterfly, etc.) |
| `yc_level_slope_curvature()` | Level, slope, and curvature factor extraction |

All `yc_curve` and `yc_pca` objects have `print()`, `summary()`, and `plot()` methods.


## Academic references

- Nelson, C.R. and Siegel, A.F. (1987). Parsimonious Modeling of Yield Curves. *The Journal of Business*, 60(4), 473-489. [doi:10.1086/296409](https://doi.org/10.1086/296409)
- Svensson, L.E.O. (1994). Estimating and Interpreting Forward Interest Rates: Sweden 1992-1994. *NBER Working Paper*, 4871. [doi:10.3386/w4871](https://doi.org/10.3386/w4871)
- Litterman, R. and Scheinkman, J. (1991). Common Factors Affecting Bond Returns. *The Journal of Fixed Income*, 1(1), 54-61. [doi:10.3905/jfi.1991.692347](https://doi.org/10.3905/jfi.1991.692347)


## Related packages

| Package | Description |
|---------|-------------|
| [boe](https://github.com/charlescoverdale/boe) | Bank of England data (includes official yield curves) |
| [fred](https://github.com/charlescoverdale/fred) | Federal Reserve Economic Data (includes Treasury rates) |
| [readecb](https://github.com/charlescoverdale/readecb) | European Central Bank data (includes euro area yield curves) |
| [debtkit](https://github.com/charlescoverdale/debtkit) | Debt sustainability analysis and fiscal projections |
| [mpshock](https://github.com/charlescoverdale/mpshock) | Monetary policy shock series (US/UK/AU) |


## Issues

Found a bug or have a feature request? Open an issue at [github.com/charlescoverdale/yieldcurves/issues](https://github.com/charlescoverdale/yieldcurves/issues).


## Keywords

r, r-package, yield-curve, fixed-income, nelson-siegel, svensson, term-structure, interest-rates, bond-math, finance, duration, convexity, z-spread
