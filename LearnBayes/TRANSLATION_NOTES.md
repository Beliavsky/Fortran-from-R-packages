# Translation notes

## Upstream

- Package: `LearnBayes`
- Upstream version: 2.15.2
- Upstream date: 2026-01-29
- Author/maintainer: Jim Albert
- Upstream license: GPL (>= 2)
- Upstream compiled code: none

The full source used for the translation is retained in
`upstream/LearnBayes-2.15.2/`.

## Translation strategy

The package consists primarily of short R implementations of Bayesian formulas,
simulation algorithms, and teaching utilities. Rather than translating every
function as an isolated procedure, the Fortran version first provides shared
probability, linear-algebra, RNG, callback, and Monte Carlo infrastructure and
then expresses the original LearnBayes routines on top of it.

All numerical code is native Fortran; no R runtime is required.

## Functional correspondence

### Discrete/posterior summaries

Translated directly or through equivalent native result types:

- `beta.select` -> `beta_select`
- `normal.select` -> `normal_select`
- `pbetap` -> `pbetap`
- `pbetat` -> `pbetat`
- `pdisc` -> `pdisc`
- `pdiscp` -> `pdiscp`
- `ctable` -> `ctable`
- `discint` -> `discint`
- `summary.bayes` -> `summarize_discrete`
- `normal.normal.mix` -> `normal_normal_mix`
- `binomial.beta.mix` -> `binomial_beta_mix`
- `poisson.gamma.mix` -> `poisson_gamma_mix`
- `prior.two.parameters` -> `prior_two_parameters`
- `histprior` -> `histprior`
- `discrete.bayes` -> `discrete_bayes`
- `discrete.bayes.2` -> `discrete_bayes_2`

### Densities and posterior kernels

Translated:

- `dmnorm`, `dmt`
- `betabinexch`, `betabinexch0`, `bfexch`
- `bradley.terry.post`
- `cauchyerrorpost`
- `groupeddatapost`
- `howardprior`
- `lbinorm`
- `logctablepost`
- `logisticpost`
- `logpoissgamma`, `logpoissnormal`
- `normchi2post`, `normnormexch`
- `poissgamexch`
- `reg.gprior.post`
- `transplantpost`
- `weibullregpost`
- `mnormt.onesided`, `mnormt.twosided`

### Generic simulation and approximation

Translated:

- `gibbs`
- `rwmetrop`
- `indepmetrop`
- `laplace`
- `impsampling`
- `rejectsampling`
- `sir`
- `simcontour`
- `rdirichlet`, `rigamma`, `rmnorm`, `rmt`

The R routines accepting arbitrary function objects are represented by typed
Fortran callback objects. Callback data/constants may be carried in allocatable
numeric components.

`laplace()` uses a native Nelder-Mead maximizer, corresponding to the default
method used by R `optim()`, followed by a central finite-difference Hessian and
the usual Gaussian/Laplace normalization. The optimizer has its own deterministic
termination details, so last-bit equality with a particular R version of
`optim()` is not claimed on difficult objectives.

### Regression

Translated:

- `blinreg`
- `blinregexpected`
- `blinregpred`
- `bayesresiduals`
- `bayes.probit`
- `bprobit.probs`
- `bayes.model.selection`

The probit implementation uses Albert-Chib latent-normal Gibbs sampling. Model
selection exhaustively enumerates subsets and evaluates the package's Laplace
approximation natively.

### Hierarchical/specialized simulation

Translated:

- `normpostsim`
- `normpostpred`
- `robustt`
- `hiergibbs`
- `ordergibbs`
- `bfindep`
- `bayes.influence`

The specialized routines retain the dimensions and priors used by their
upstream teaching examples where those are intrinsic to the original function.

### Data/presentation calculations

Translated computationally:

- `regroup`
- `careertraj.setup` -> `careertraj_setup`
- `mycontour` -> `contour_grid`
- `triplot` -> `triplot_data`

`careertraj_setup` uses integer player identifiers instead of R character/factor
labels.

`predplot()` is a plotting wrapper around `pbetap()` and therefore needs no
separate numerical kernel. `plot.bayes`, `plot.bayes2`, and `print.bayes` are
presentation-only and are omitted.

## Numerical infrastructure

The translation supplies native implementations of:

- normal CDF and quantile;
- regularized incomplete beta;
- beta log function;
- type-7 sample quantiles;
- normal/gamma/beta/Student-t/Poisson/binomial log densities or masses;
- dense Cholesky factorization, SPD determinant, linear solve/inverse, and OLS;
- deterministic uniform/normal/gamma/chi-square/discrete RNGs;
- multivariate normal/t simulation and density evaluation.

These avoid dependencies on BLAS/LAPACK or external statistics libraries and
make the FPM package self-contained. For very large dense linear-algebra
problems, a future optional BLAS/LAPACK backend would be a performance
improvement rather than a missing LearnBayes feature.

## Random-number parity

The Fortran port deliberately uses its own deterministic RNG state rather than
R's RNG. Distributional algorithms follow the same statistical models, but a
seed value does not produce R's exact sequence. Tests therefore use analytic
fixtures where possible and fixed native seeds for reproducible Monte Carlo
smoke checks.

## `rtruncated`

R's `rtruncated()` takes arbitrary CDF and quantile functions through dynamic
function objects and `...`. A direct equivalent would be awkward and unsafe in
a statically typed API. The translated library provides `rtruncated_normal`,
which covers the package's internal order-restricted Gibbs use. Typed custom
interfaces can be added if another concrete distribution is required.

## Validation

The strict suite includes:

- analytic probability/distribution references;
- the retained LearnBayes posterior/Bayes-factor formula fixtures;
- Laplace approximation of a normalized standard Gaussian;
- random-walk, coordinate Gibbs, independence Metropolis, importance,
  rejection, SIR, and contour simulation checks;
- regression, prediction, Bayesian residual, probit, and model-selection tests;
- hierarchical, robust-t, order-restricted, independence-BF, and influence
  tests;
- mixture, discrete Bayes, regroup, career-trajectory, contour, and triplot
  calculations.

The maintained code is also built with bounds checking and traps for invalid,
divide-by-zero, and overflow floating-point exceptions.
