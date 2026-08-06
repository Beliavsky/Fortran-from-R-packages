# Translation coverage

## Translated computational exports

- `bond_future_price`
- `stir_future_price`
- `stir_rate`
- `ctd_bond_yield`
- `proba_ctd`
- `proba_ctd_opt`
- `bond_yield_spread`

## Supporting numerical code

- Two- and three-component lognormal-mixture PDFs, CDFs, quantiles, and moments
- European, American, and futures-style-margin option pricing
- Multi-start bounded Nelder-Mead calibration
- Gregorian dates and four day-count conventions
- Coupon schedules, accrued interest, dirty prices, yield solving, and carry
- CTD self-consistency, net-basis aggregation, and probability integration
- Gaussian-copula simulation and Gaussian KDE

## Intentionally omitted exports

The following functions only retrieve Bloomberg data through `Rblpapi` and do
not contain portable numerical algorithms:

- `bond_future_charac_bbg`
- `stir_future_charac_bbg`
- `deliv_bonds_charac_bbg`
- `option_prices_bbg`

All `ggplot2` plotting code and R container/interface plumbing are also omitted.
