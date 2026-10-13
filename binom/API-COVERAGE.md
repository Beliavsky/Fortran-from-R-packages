# Computational API coverage

Package status: **substantial**. **18 of 18 (100%)** distinct exported
computational R functions are mapped. This fraction measures mapped functions,
not complete R compatibility. The denominator excludes the presentation-only
exports `binom.bayes.densityplot`, `tkbinom.power` and `binom.plot`.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
|---|---|---|---|
| `binom.asymp` | `binom_intervals` | `binom_asymp` | complete |
| `binom.agresti.coull` | `binom_intervals` | `binom_agresti_coull` | complete |
| `binom.exact` | `binom_intervals` | `binom_exact` | complete |
| `binom.wilson` | `binom_intervals` | `binom_wilson` | complete |
| `binom.prop.test` | `binom_intervals` | `binom_prop_test` | substantial |
| `binom.bayes` | `binom_intervals` | `binom_bayes` | substantial |
| `binom.cloglog` | `binom_intervals` | `binom_cloglog` | complete |
| `binom.probit` | `binom_intervals` | `binom_probit` | complete |
| `binom.profile` | `binom_profile_module` | `binom_profile` | substantial |
| `binom.confint` | `binom_intervals` | `binom_confint`, `binom_confint_all` | substantial |
| `binom.coverage` | `binom_design` | `binom_coverage` | substantial |
| `binom.lrt` | `binom_intervals` | `binom_lrt` | substantial |
| `binom.logit` | `binom_intervals` | `binom_logit` | complete |
| `cloglog.sample.size` | `binom_design` | `cloglog_sample_size` | substantial |
| `binom.power` | `binom_design` | `binom_power` | substantial |
| `binom.sim` | `binom_design` | `binom_sim`, `binom_sim_counts` | substantial |
| `binom.length` | `binom_design` | `binom_length` | substantial |
| `binom.optim` | `binom_optimization` | `binom_optim` | substantial |

## Mapping notes

- **binom.asymp**: Wald interval; the unbounded limits are deliberately not clipped.
- **binom.agresti.coull**: Agresti-Coull adjusted-count interval, without clipping.
- **binom.exact**: Clopper-Pearson interval with explicit boundary values and stable beta inversion.
- **binom.wilson**: Wilson score interval, with exact 0/1 observed-boundary limits to correct floating-point endpoint coverage artifacts.
- **binom.prop.test**: One-sample two-sided stats::prop.test confidence interval, including null-dependent continuity correction; R warnings omitted.
- **binom.bayes**: Central and highest-density beta intervals, prior shapes, tolerance and iteration cap. Pure equal-density tail-allocation solver replaces C density-height roots; R data frames and warnings omitted.
- **binom.cloglog**: Complementary-log-log interval and exact boundary correction.
- **binom.probit**: Probit interval and exact boundary correction.
- **binom.profile**: Sampled signed deviance, FMM spline, three-times-knot grid and inverse linear interpolation; maxsteps, del and Bayesian edge adjustment. Plotting and R objects omitted; explicit resource and bracket status.
- **binom.confint**: All eleven numerical methods. Scalar/elemental selected-method API plus all-method array; no partial-name matching, recycling, warnings, data frames or ellipsis forwarding.
- **binom.coverage**: Exact binomial-mass enumeration, inclusive endpoints. One named method per call; R callbacks, implicit probability grid and table construction omitted.
- **binom.lrt**: Likelihood-ratio roots, Bayesian pseudocounts, conf.adj and explicit tol/maxit. Upstream Brent defaults and boundary clipping retained; stable likelihood differences avoid common-coefficient cancellation. Plotting and R object handling omitted.
- **binom.logit**: Logit interval and exact boundary correction.
- **cloglog.sample.size**: Sample size, attained power and detectable-alternative branches, alternatives, dispersion, real-valued rounding and exact.n. One parameter combination per elemental call; R expand.grid table omitted.
- **binom.power**: All six power methods and three alternatives; retains upstream beta-tail formula for exact/lrt. Scalar/elemental evaluation avoids R vector-recycling bug; nonlinear links require interior probabilities, while asymptotic/exact/lrt include endpoints.
- **binom.sim**: Full capture-count confidence interval and mean-width calculation from caller-supplied uniforms or binomial counts. Caller controls RNG; no R RNG stream or implicit M/default draws.
- **binom.length**: Exact binomial-weighted interval length for every named method; R callbacks and expand.grid output omitted.
- **binom.optim**: Symmetric confidence optimization, beta-integrated squared coverage loss and shared L-BFGS-B with R-style fixed finite differences. Supports k, p0, transform, start and bounded solver controls. Raw closed-form trial algebra, partial-NaN removal and empty-count integration failures follow R; all-invalid zero-loss artifacts are rejected and invalid final intervals return failure. Public intervals retain (0,1). No R callbacks, plotting or general optim scaling arguments.
