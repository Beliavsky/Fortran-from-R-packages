# Porting notes

## API translation

R list results are represented by Fortran derived types, including
`lmcd_result`, `auto_vr_result`, `bootstrap_result`, `wright_result`,
`dl_result`, and `chen_deo_result`. Variable-length outputs are allocatable.

The Fortran routine names are descriptive snake-case names rather than names
containing R punctuation. `use vrtest` imports the complete public API.

## Numerical dependencies

The project is self-contained. It includes:

- inverse normal CDF using the Acklam rational approximation and refinement;
- regularized incomplete gamma and chi-square CDF/quantile calculations;
- R type-7 sample quantiles;
- average ranks for tied observations;
- partial-pivoting linear solves and matrix inversion;
- normal, Mammen, and Rademacher random weights;
- Fisher-Yates random permutations;
- Simpson integration and stable log-mean-exp evaluation.

## Deliberate corrections and safeguards

### AR1 indexing

The upstream `AR1.R` sets `T <- length(x) - 1`, then uses `x[2:T]` and
`x[1:(T-1)]`. This omits the final adjacent observation pair. `ar1_fit`
uses all adjacent pairs by default. Passing `legacy_indexing=.true.` reproduces
the upstream indexing, and `automatic_variance_ratio` exposes the corresponding
`legacy_ar1` option.

### Ordinary bootstrap indices

The upstream ordinary bootstrap creates indices with
`as.integer(runif(n, min=1, max=n))`, which never selects the final observation
because truncation is applied to a half-open uniform draw. The Fortran routine
samples all indices from 1 through n.

### Wald covariance matrix

The upstream nested loop overwrites symmetric entries repeatedly. The Fortran
implementation evaluates the final intended symmetric formula directly using
the smaller and larger holding periods.

### Generalized spectral bootstrap

As in the R routine, bootstrap statistics reuse the weight matrix constructed
from the original observations.

### Numerical stability

The average exponential test uses log-mean-exp rather than directly averaging
large exponentials. Matrix solves return an integer `solve_info` field instead
of continuing with a singular system. Small denominators and degenerate input
series are guarded to prevent division by zero.

## Behavioral differences

- Bootstrap streams are reproducible within a compiler/runtime after
  `seed_random`, but they are not bit-identical to R's RNG.
- Exact p-values from bootstraps depend on the selected seed and iteration count.
- `generalized_spectral_test` returns the statistic and bootstrap critical
  values in addition to the p-value; the R function computes these values but
  returns only the p-value.
- `variance_ratio_curve` returns values and confidence bands without drawing a
  graph.
