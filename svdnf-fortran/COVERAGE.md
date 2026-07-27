# Computational Coverage

## Exported R API

| Upstream routine | Fortran implementation | Status |
|---|---|---|
| `dynamicsSVM` | `dynamics_svm`, `set_custom_dynamics` | Implemented |
| `DNF.dynamicsSVM` | `dnf_filter`, alias `dnf` | Implemented |
| `DNFOptim.dynamicsSVM` | `dnf_optimize`, alias `dnf_optim` | Implemented for built-ins |
| custom `DNFOptim` | `dnf_optimize_custom` plus parameter-setter callback | Implemented |
| `modelSim.dynamicsSVM` | `model_simulate`, alias `model_sim` | Implemented |
| `extractVolPerc.SVDNF` | `extract_vol_percentile`, alias `extract_vol_perc` | Implemented |
| `predict.SVDNF` | `predict_filter`, alias `predict_svdnf` | Implemented |
| `predict.DNFOptim` | call `predict_filter` on `optimization_result%filter` | Implemented |
| `pars.dynamicsSVM` | `model_parameter_names`, `parameter_vector` | Implemented |
| `gridMaker.*` | `grid_maker` | Implemented for every built-in model |
| `probCalculator` | `probability_components` | Implemented |
| `initGuess.dynamicsSVM` | `initial_guess` | Implemented with native heuristics |
| `logLik.*` | `filter_result%log_likelihood` | Implemented as typed field |
| `summary.DNFOptim` | Hessian and `standard_errors` in `optimization_result` | Numerical content implemented |

## Built-in models

- Duffie-Pan-Singleton stochastic variance, return jumps, and variance jumps.
- Bates stochastic variance and return jumps.
- Heston stochastic variance.
- Pitt-Malik-Doucet log-volatility with Bernoulli jumps.
- Taylor log-volatility.
- Taylor log-volatility with leverage.
- CAPM-SV factor-adjusted log-volatility.

## Internal Rcpp helpers

The operations performed by `Cpp_rowSums_modN`, `Cpp_prodfun`, and
`dnorm_cpp_prod` are incorporated directly into `dnf_filter` and
`probability_components`. They are not exposed as separate public functions
because the Fortran implementation performs the same contraction without
constructing R-style flattened `expand.grid` arrays.

## Numerical support

The project includes self-contained implementations of:

- normal density, CDF, and inverse CDF;
- gamma density, CDF, and random generation;
- Poisson and binomial probabilities;
- Poisson, Bernoulli, normal, and gamma random generation;
- discrete weighted sampling;
- Nelder-Mead optimization;
- central-difference Hessians and matrix inversion.

## Excluded presentation/integration code

- S3 print and plot methods;
- `xts` and `zoo` date alignment;
- graphics and confidence-band plotting;
- R formula/list introspection;
- vignette build artifacts.

Their underlying numerical quantities are available from typed results.
