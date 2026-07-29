# Porting notes

## Translation policy

The port retains numerical algorithms and examples that make sense as a
Fortran library. R-only plotting, graphics annotations, S3-style argument
routing through `...`, matrix orientation workarounds, and console text are not
translated. Scalar and array dimensions are explicit in the Fortran API.

The original package has 147 exports. Fifty-two are plot/figure/QQ/mean-excess
presentation functions. Their underlying numerical calculations are exposed by
the Fortran routines, but the graphics wrappers themselves are omitted.

## Main mappings

| R family | Fortran API |
|---|---|
| `NormalVaR`, `NormalES` | `normal_var`, `normal_es` |
| `tVaR`, `tES` | `student_t_var`, `student_t_es` |
| `LogNormalVaR`, `LogNormalES` | `lognormal_var`, `lognormal_es` |
| `LogtVaR`, `LogtES` | `log_student_t_var`, `log_student_t_es` |
| `HSVaR`, `HSES` | `historical_var`, `historical_es` |
| kernel VaR/ES variants | `kernel_var`, `kernel_es` with a kernel constant |
| bootstrap VaR/ES and intervals | `bootstrap_var_es`, `bootstrap_confidence_interval` |
| `BoxCoxVaR`, `BoxCoxES` | `boxcox_var`, `boxcox_es` |
| Cornish-Fisher functions | `cornish_fisher_var`, `cornish_fisher_es` |
| Gumbel/Frechet/GPD functions | corresponding routines in `dowd_risk` |
| Hill/Pickands functions | `hill_estimator`, `hill_quantile_estimator`, `pickands_estimator` |
| variance-covariance functions | routines in `dowd_portfolio` |
| normal/adjusted hotspots | leave-one-position-out hotspot subroutines |
| `PCAPrelim`, `PCAVaR`, `PCAES` | `pca_prelim`, `pca_var`, `pca_es` |
| backtests | routines in `dowd_backtests` |
| Black-Scholes functions | routines in `dowd_options` |
| American put functions | binomial price, VaR, ES, and simulation routines |
| copula functions | routines in `dowd_copulas` |
| insurance/pension/bond/filter examples | routines in `dowd_simulations` |
| spectral normal risk | `spectral_risk_normal` |
| quantile standard errors | normal and Student-t standard-error routines |

## Deliberate numerical corrections

The original package description warns that inherited toolbox errors remain.
This port does not reproduce obvious defects when a standard formula is clear.
Notable corrections include:

1. Portfolio variance is converted to standard deviation with `sqrt(p' C p)`.
   Several R routines multiply by `p' C p` directly as though it were a
   standard deviation.
2. Christoffersen independence uses the full binary exception sequence, the
   supplied VaR vector, and the standard likelihood-ratio expression. The R
   source replaces the VaR input with the P/L input, removes nonexceptions
   before transition counting, and has malformed likelihood arithmetic.
3. Historical ES uses an integral of interpolated empirical quantiles. This
   avoids the R source's repeated `ceiling(index)` lower-index error.
4. The Hill estimator uses the standard average log-excess formula. The R
   implementation repeatedly adds the same order statistic and has shifted
   indexing.
5. The Gumbel copula is evaluated as
   `exp(-(((-log u)^beta + (-log v)^beta)^(1/beta)))`.
6. Gaussian-copula normal-margin sums use their exact normal distribution and
   preserve the requested correlation. The R sum routine passes zero
   correlation internally.
7. American-put node outcomes are weighted by their binomial probabilities for
   VaR and ES. The R code treats the `m+1` nodes as equally likely.
8. European option routines consistently convert day inputs to years before
   Black-Scholes valuation.
9. Student-t ES uses the closed-form tail mean rather than a coarse average of
   1,000 quantiles.
10. Box-Cox lambda estimation is implemented internally by profile likelihood,
    replacing the R `forecast` dependency.

These corrections mean exact agreement with erroneous R outputs is not always
a design goal. Standard textbook identities and independent reference values
are tested instead.

## Omitted functionality

- All plotting and figure routines
- R graphics side effects in statistical simulation routines
- R package datasets and help-page examples that only create plots
- The `*DFPerc` presentation family, whose source contains several argument
  assignment defects; the reusable normal/Student-t quantile standard errors
  and bootstrap confidence intervals are provided instead
- Dynamic `...` argument dispatch; use explicit Fortran arguments

## Naming and units

Fortran names use lower-case snake case. Holding periods are in days only for
functions whose names or documentation explicitly say `*_days`; ordinary risk
functions accept a holding-period multiplier in the same units as the input
mean and volatility. Black-Scholes pricing accepts maturity in years, while the
Dowd-style option risk wrappers accept days and convert internally.
