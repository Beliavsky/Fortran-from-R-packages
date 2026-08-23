# Translation coverage

## Implemented in v0.1.0

- Generic parameter object and family dispatch
- PDF/CDF/quantile/RNG for the main lmomco distribution catalog listed in README.md
- Finite-parameter Rice, Eta-Mu, and Kappa-Mu numerical evaluation
- Sample PWMs and L-moments
- L-moment ratios
- Quantile-integrated theoretical L-moments
- Direct Normal, Exponential, Gumbel, and GEV L-moment fitting
- Named compatibility wrappers for principal distributions
- Return-period, plotting-position, harmonic-mean, and Gini helpers

## Not yet translated

High-value remaining numerical targets include:

- right-censored L-moments and BFR workflows
- analytical/bootstrap L-moment covariance routines
- asymmetric/trimmed L-moment catalog beyond ordinary L-moments
- full L-comoment and regional-frequency-analysis workflow parity
- all family-specific `par*`, `lmom*`, and validation wrappers
- MLE/MPS fitting framework across every distribution
- specialized Eta-Mu/Kappa-Mu parameter-table inversion workflows
- `kappa = Inf` Kappa-Mu Dirac-mass limit
- specialized order-statistic/bootstrap confidence procedures

Plotting and R S3/list/formula presentation infrastructure are intentionally outside the Fortran scope.
