# Computational coverage

## Exported R routines

| Upstream routine | Fortran implementation | Status |
|---|---|---|
| `BS_European_Greeks` | `bs_european_greeks` | Complete numerical translation |
| `BS_Geometric_Asian_Greeks` | `bs_geometric_asian_greeks` | Complete numerical translation |
| `BS_Implied_Volatility` | `bs_implied_volatility` | Complete, with safeguarded iteration |
| `BS_Malliavin_Asian_Greeks` | `bs_malliavin_asian_greeks` | Complete scalar-parameter translation |
| `Binomial_American_Greeks` | `binomial_american_greeks` | Complete numerical translation |
| `Greeks` | `calculate_greeks` plus the explicit Malliavin routines | Computational dispatch translated |
| `Implied_Volatility` | `implied_volatility` and model-specific wrappers | Complete callback-based numerical API |
| `Malliavin_Asian_Greeks` | `malliavin_asian_greeks` | Complete numerical translation |
| `Malliavin_European_Greeks` | `malliavin_european_greeks` | Complete numerical translation |
| `Malliavin_Geometric_Asian_Greeks` | `malliavin_geometric_asian_greeks` | Complete numerical translation |
| `Greeks_UI` | Not compiled | Shiny/plotly user interface |

## Internal C++ helpers

The behavior of the following Rcpp helpers is incorporated into
`greeks_paths` and `greeks_malliavin`:

- `make_BM`
- `calc_X`
- `calc_log_X`
- `calc_I`, `calc_I_1`, `calc_I_2`, and `calc_I_3`
- `calc_XW` and `calc_tXW`
- `rowCumsums`

## Models and payoffs

### European Black-Scholes

- Call and put
- Cash-or-nothing call and put
- Asset-or-nothing call and put
- All 18 upstream sensitivities

### Geometric Asian Black-Scholes

- Call and put
- Fair value, delta, rho, vega, theta, gamma, and vomma

### American Black-Scholes

- Call and put
- Corrected binomial fair value
- Delta, gamma, vega, theta, rho, and dividend-yield sensitivity

### Malliavin Monte Carlo

- European call, put, cash binary, and asset binary payoffs
- Arithmetic and geometric Asian payoffs
- User-defined payoff callback
- User-defined payoff-derivative callback
- Black-Scholes diffusion
- Compound-Poisson jump diffusion
- User-defined jump sampler
- Antithetic Brownian paths where supported upstream
- Pathwise standard errors

## Excluded infrastructure

- Shiny application and reactive controls
- `ggplot2`, `plotly`, and other graphics
- `tibble`, `tidyr`, and `magrittr` presentation pipelines
- R S3/Roxygen documentation machinery
- Rcpp registration and R package-loading glue
- Vignette rendering

The original source, vignettes, images, and R metadata are retained under
`original/greeks-1.5.6` for provenance.
