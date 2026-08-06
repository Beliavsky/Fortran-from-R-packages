# Porting notes

## Formula parity

The default calculations preserve the upstream R formulas:

- duration and convexity are evaluated at the current yield;
- yield income is evaluated from the lagged yield;
- the first return is missing;
- a user-supplied duration or convexity series overrides the internally
  calculated values;
- the default annualization scale is 261.

## Numerically stable mode

The direct formulas lose precision near zero and are undefined at exactly zero.
The optional corrected mode uses series expansions near zero and stable
logarithm/exponential calculations elsewhere. This does not change ordinary
positive-yield results materially.

## Missing observations

`carry_forward` reproduces `zoo::na.locf0` behavior for interior and trailing
NaNs. Leading NaNs remain missing. `total_return` itself does not fill missing
values; it propagates them through the lagged-return calculation, as the
upstream function does when called directly.

## Omitted R infrastructure

The following are not numerical algorithms and were not translated:

- FRED network access through `quantmod::getSymbols`;
- `xts`, `tibble`, and date-index conversions;
- dplyr pipelines and package namespace hooks;
- vignette plotting and third-party performance-reporting examples.

A calling program can load yields from CSV, a database, or an HTTP client and
then pass the numeric observations to `prepare_yields` and `total_return`.
