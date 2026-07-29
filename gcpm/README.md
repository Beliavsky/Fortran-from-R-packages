# GCPM modern Fortran translation

This is a modern Fortran/FPM translation of the computational parts of the R
package **GCPM 1.2.2**, Generalized Credit Portfolio Model.

The original package was written by Kevin Jakob, with the startup copyright
notice naming Kevin Jakob and Dr. Matthias Fischer. The original package is
licensed under GNU GPL version 2. This translation is distributed under the
same GPL version 2 terms; see `LICENSE` and `ORIGIN.md`.

## Implemented scope

The library implements:

- portfolio validation and loss discretization;
- analytical expected loss and standard-deviation decomposition;
- analytical CreditRisk+ loss-distribution recursion;
- CRP Monte Carlo with Bernoulli or Poisson defaults;
- CreditMetrics-style Monte Carlo with Bernoulli or Poisson defaults;
- scenario likelihood ratios and scenario recycling;
- empirical loss-distribution construction;
- expected loss, standard deviation, VaR, economic capital, and expected
  shortfall;
- analytical CreditRisk+ contributions to EC, VaR, ES, and standard deviation;
- simulation-based contributions to EC, VaR, and ES;
- standard normal CDF and inverse CDF, Poisson, gamma, and normal random-number
  support;
- a reader for the original GCPM semicolon-delimited portfolio format,
  including decimal-comma numeric fields.

R-specific S4 classes, plotting, progress bars, R export methods, and R
parallel-cluster wrappers are not part of the Fortran library. The numerical
algorithms behind the package are included. The current Monte Carlo driver is
serial and portable; callers can parallelize independent simulation batches.

## Requirements

- a Fortran 2018 compiler;
- Fortran Package Manager (`fpm`) for the normal build workflow.

No external Fortran dependencies are required.

## Build and test

```text
fpm build
fpm test
fpm run --example demo_gcpm
```

The translation was also compiled and checked directly with GNU Fortran 14.2
using runtime bounds and consistency checks.

## Basic workflow

```fortran
use gcpm

 type(credit_portfolio) :: portfolio
 type(gcpm_model) :: model
 type(loss_distribution) :: distribution
 integer :: status
 character(len=:), allocatable :: message

 call allocate_portfolio(portfolio, n_counterparties=100, n_sectors=3)
 ! Fill EAD, LGD, PD, default_kind, and weight.

 model%model_kind = model_analytical_crp
 model%link_kind = link_crp
 model%loss_unit = 1000.0_dp
 model%alpha_max = 0.9999_dp
 allocate(model%sector_variance(3))
 model%sector_variance = [0.8_dp, 1.0_dp, 1.2_dp]

 call calculate_portfolio_statistics(portfolio, model, status, message)
 call analytical_creditrisk_plus(portfolio, model, distribution, status, message)
```

See `example/demo_gcpm.f90` for complete analytical and Monte Carlo examples.

## Main modules

- `gcpm_types`: portfolio, model, distribution, and risk-result types.
- `gcpm_portfolio`: validation, discretization, sector statistics, and SD
  contributions.
- `gcpm_analytical`: analytical CreditRisk+ recursion.
- `gcpm_simulation`: CRP and CreditMetrics Monte Carlo engines.
- `gcpm_risk`: VaR, EC, ES, threshold search, and risk contributions.
- `gcpm_csv`: original-format portfolio reader.
- `gcpm_math`: distribution functions and random-number utilities.
- `gcpm`: umbrella module re-exporting the public API.

## Input conventions

`read_gcpm_portfolio` expects the original column order:

```text
Number;Name;Business;Country;EAD;LGD;PD;Default;Sector1;Sector2;...
```

The default separator is `;`, and decimal commas are accepted. Business and
country fields are read but not retained because they do not enter the
numerical model. The original compressed example portfolios are included in
`data/`, together with a small uncompressed test file.

For simulation, factor rows are scenarios and factor columns are sectors. If
`model%n_simulations` is zero, all supplied rows are used once. A positive
value recycles rows as necessary, matching the behavior of the R package.

## Numerical notes

- Potential losses are rounded to positive integer multiples of `loss_unit`.
- The transformed default probability preserves expected loss after
  discretization.
- Analytical CreditRisk+ uses Poisson default counts, as in the original
  analytical implementation.
- The analytical recursion stops when `alpha_max` is reached or the optional
  exposure-band limit is exhausted.
- Simulation risk contributions require a finite `loss_threshold` and enough
  `max_stored_scenarios` capacity to retain the required tail scenarios.
- CreditMetrics simulation requires each counterparty's systematic variance
  `w' Sigma w` to be less than one.

## Original references

- Jakob, K. and Fischer, M., "GCPM: A flexible package to explore credit
  portfolio risk," Austrian Journal of Statistics 45.1 (2016), 25-44.
- J.P. Morgan, *CreditMetrics Technical Document*, 1997.
- Credit Suisse First Boston, *CreditRisk+*, 1997.
- Gundlach and Lehrbass, *CreditRisk+ in the Banking Industry*, 2003.
