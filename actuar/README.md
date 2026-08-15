# actuar-fortran

Modern Fortran/FPM translation of the computational core of the R package
`actuar` 3.3-7.

The upstream package is by Vincent Goulet and contributors and is licensed
under GPL-2.0-or-later. This translation preserves that license and retains
selected upstream R/C sources under `upstream/` for algorithm provenance.

## Scope of v0.3.0

v0.3 retains the v0.1-v0.2 distribution, aggregate-risk, ruin, grouped-data
and credibility functionality and closes the principal computational gaps
identified after v0.2.

### Minimum-distance estimation

The matrix-first MDE API translates the numerical objectives in upstream
`mde()`:

- `mde_cvm` for individual Cramer-von Mises fitting;
- `mde_grouped_cvm` for grouped-data Cramer-von Mises fitting;
- `mde_grouped_chisq` for the modified chi-square criterion;
- `mde_grouped_las` for layer-average-severity fitting.

The model CDF or limited-expected-value function is supplied as a Fortran
procedure. Scalar/vector weights and parameter bounds are supported. A
self-contained bounded Nelder-Mead optimizer is used in this release; the
upstream R objective definitions themselves are preserved.

### Coverage transformations

`coverage_spec_t`, `coverage_cdf` and `coverage_pdf` implement the numerical
behavior of upstream `coverage()` for combinations of:

- ordinary or franchise deductibles;
- policy limits;
- coinsurance;
- inflation;
- per-loss or per-payment variables.

`coverage_pdf` deliberately returns probability mass at the zero-payment and
policy-limit endpoints when the transformed law is mixed, matching the R
helper's density/mass convention.

### Hachemeister barycenter model

`hachemeister_barycenter_fit` implements the upstream barycenter variant:

- average period weights across contracts;
- weighted QR/orthogonalization of the regression design;
- contract-specific weighted least-squares fits;
- diagonal between-contract variance estimation in the orthogonal basis;
- unbiased or iterative Buhlmann-Straub-style variance estimation;
- diagonal credibility matrices;
- transition of collective, individual and adjusted coefficients back to the
  original design basis.

The original-basis and orthogonal-basis coefficients are both retained in the
result type.

### Exact iterative hierarchical credibility

`hierarc_exact_fit` is a direct Fortran translation of the numerical recursion
split between upstream `R/hierarc.R` and compiled `src/hierarc.c`.

It supports:

- Buhlmann-Gisler structure-parameter estimators;
- Ohlsson estimators;
- the compiled iterative estimator;
- arbitrary hierarchy depth from an integer classification matrix;
- exact parent mappings, level weights, weighted means and credibility
  factors;
- recursively calculated credibility premiums.

In the iterative method the Fortran port preserves the upstream compiled rule
that aggregation uses a child's credibility factor when nonzero and otherwise
falls back to its current/natural weight.

### Hierarchical portfolio simulation

`rcomphierarc_simulate` is the matrix/callback-driven Fortran counterpart of
upstream `rcomphierarc()`.

Instead of evaluating R expressions, users provide:

- per-level node counts;
- optional frequency mixing callbacks;
- a terminal frequency callback;
- optional severity mixing callbacks;
- a terminal severity callback.

The returned `hierarchical_portfolio_t` exposes terminal paths, frequencies,
flattened claim severities, claim-to-node assignments and aggregate claims per
terminal node. The callback interface propagates all previously generated
mixing parameters down the appropriate hierarchy.

## Existing functionality retained from v0.1-v0.2

The package continues to include:

- heavy-tailed loss families including Feller-Pareto, transformed-beta/gamma,
  Burr, Pareto variants and inverse families;
- raw and limited moments;
- zero-truncated/zero-modified count laws and Poisson-inverse-Gaussian;
- phase-type distributions;
- Panjer recursion and exact aggregate convolution;
- normal and normal-power aggregate approximations;
- compound/mixture simulation;
- VaR and CTE;
- phase-type Cramer-Lundberg/Sparre-Andersen ruin algorithms;
- Buhlmann-Straub, Bayesian, hierarchical and Hachemeister-origin credibility;
- grouped moments, ogive, empirical LEV and grouped quantiles.

The supplied `expint-fortran` translation remains a vendored FPM dependency
for extended incomplete-gamma calculations.

## Umbrella API

```fortran
use actuar
```

## Building

With FPM:

```text
fpm build
fpm test
fpm run --example v03_remaining
```

FPM was unavailable in the translation environment. The complete package and
vendored dependency were therefore compiled directly with GNU Fortran using:

```text
-std=f2018 -Werror=implicit-interface -Werror=trampolines -fcheck=all -O0
```

No unlimited-free-form-line-length compiler option is required.

## Deliberate API differences

R S3 classes, formula evaluation, plotting, printing and expression evaluation
are not replicated. `mde`, `coverage` and `rcomphierarc` therefore use explicit
Fortran callbacks and arrays instead of R closures/calls. This keeps the
numerical algorithms portable while preserving their statistical meaning.

The large computational targets listed after v0.2 are now represented. Any
remaining gaps are primarily R convenience/dispatch layers and specialized
aliases rather than distinct major numerical engines.
