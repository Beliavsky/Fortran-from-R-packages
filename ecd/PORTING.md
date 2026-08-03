# Porting notes

## Scope

The R package declares 279 `export(...)` entries, but many are aliases, S4/S3
methods, vector/MPFR dispatch wrappers, database operations, data loaders,
plots, or report helpers. This port translates the reusable numerical content
rather than creating hundreds of one-line Fortran aliases.

A common Fortran procedure is used when several R exports select the same
formula, integration route, precision mode, or plotting-data calculation.

## Rmpfr replacement

The original package can route selected functions through Rmpfr. This port is
strictly double precision. Numerical stability is obtained through:

- transformed finite/infinite adaptive quadrature
- incomplete-gamma algorithms and root inversion
- scaled `erfcx` and quartic `erfq` evaluations
- closed-form quartic MGF, IMGF, and OGF expressions
- safeguarded bracketed roots
- explicit support and convergence checks

Results can differ for parameter combinations for which the R package requires
hundreds of bits of precision.

## Numerical substitutions

- `optimx` is replaced by a bounded Nelder-Mead optimizer.
- `stabledist` random generation is replaced by a portable
  Chambers-Mallows-Stuck implementation.
- GSL Dawson and related functions are replaced by self-contained
  approximations and quadrature.
- R `integrate`, `uniroot`, `pgamma`, and `qgamma` calls are replaced by native
  adaptive integration, root solving, and incomplete-gamma algorithms.
- R global RNG state is replaced by an explicit `rng_state` argument.

Random streams therefore do not reproduce R bit for bit.

## Distribution conventions

The Fortran port retains the package parameterizations for:

- elliptic/cusp and lambda-only ECD models
- ECLD and SGED
- standardized Lihn-Laplace
- stable-count distributions
- SLD/QSLD
- Levy lambda/skewed mixtures
- LAMP random walks

The source package often exposes the same distribution through several
specialized wrapper names. The Fortran API uses typed models and optional
arguments instead.

## Quartic option calculations

A naive infinite-tail integration of `exp(x) f(x)` over a lambda-4 distribution
can overflow even when the package's finite-precision option construction is
well defined. The port therefore translates the original closed-form quartic
MGF, IMGF, and OGF formulas based on `erfq`. The quartic V x OGF composite and
Q/Qp smile calculations use these formulas directly.

## Fitting

The port supplies bounded maximum-likelihood fitting for ECD, ECLD, and SLD
models. The optimizer and stopping path differ from `optimx`, so fitted values
can vary slightly, especially for flat or multimodal likelihoods.

The routine named `fit_ecld_moments` uses moments only to initialize a full
likelihood fit. `fit_ecld_mle` is provided as a clearer compatibility alias.

## Omitted R infrastructure

Not translated:

- S4/S3 classes, method dispatch, and string formatting
- `ecdb` SQLite persistence and transaction helpers
- YAML/CSV market-symbol configuration and package data objects
- xts/zoo index metadata
- ggplot2/base plotting and xtable output
- multicore wrapper selection
- Rmpfr class conversion helpers

The complete original tree is included under `original/ecd-master` so omitted
infrastructure and provenance remain available.

## Source corrections and hardening

The Fortran implementation includes several numerical hardening changes:

- Infinite integration endpoints are represented portably rather than by
  arithmetic with `huge()` sentinels.
- Recursive root and quadrature paths are explicitly declared recursive.
- Quartic nested root searches are restructured into standard-conforming
  module procedures.
- Exponential-tail products use analytic quartic formulas where the direct
  product can overflow.
- Invalid optimization probes return a finite penalty instead of triggering a
  floating-point exception.
