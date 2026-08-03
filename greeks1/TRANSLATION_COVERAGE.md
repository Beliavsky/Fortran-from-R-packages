# Translation coverage

## Exported computational R functions

| R export | Fortran equivalent | Status |
|---|---|---|
| `BS_European_Greeks` | `bs_european_greeks` | Translated |
| `BS_Geometric_Asian_Greeks` | `bs_geometric_asian_greeks` | Translated |
| `BS_Implied_Volatility` | `bs_implied_volatility` | Translated |
| `BS_Malliavin_Asian_Greeks` | `bs_malliavin_asian_greeks` | Translated |
| `Binomial_American_Greeks` | `binomial_american_greeks` | Translated |
| `Greeks` | `option_greeks` | Translated |
| `Implied_Volatility` | `implied_volatility` | Translated |
| `Malliavin_Asian_Greeks` | `malliavin_asian_greeks` | Translated |
| `Malliavin_European_Greeks` | `malliavin_european_greeks` | Translated |
| `Malliavin_Geometric_Asian_Greeks` | `malliavin_geometric_asian_greeks` | Translated |

## Internal C++ helpers

`make_BM`, `rowCumsums`, `calc_I`, `calc_I_1`, `calc_I_2`, `calc_I_3`,
`calc_X`, `calc_log_X`, `calc_XW`, and `calc_tXW` are translated.

## Deliberately omitted

- `Greeks_UI`, Shiny server/UI code, plotly, ggplot2, and image assets.
- The magrittr pipe export.
- Rcpp registration and R package glue.
- R documentation rendering, vignettes, and web-only interactivity.

These omissions do not remove pricing or Greek-estimation algorithms.
