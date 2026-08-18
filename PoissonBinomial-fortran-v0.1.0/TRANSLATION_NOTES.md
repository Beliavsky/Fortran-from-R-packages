# Translation notes

## Scope

Source translated: R package `PoissonBinomial` 1.2.8 (2026-03-06).

The upstream exported computational API consists of:

- `dpbinom`, `ppbinom`, `qpbinom`, `rpbinom`
- `dgpbinom`, `pgpbinom`, `qgpbinom`, `rgpbinom`

All of those computational capabilities are represented in the Fortran API,
including vector queries, integer multiplicity weights, lower/upper tails,
log probabilities, complete support tables, and both RNG strategies.

R documentation, Rcpp registration, package unload hooks, and R object
presentation code are not translated.

## Ordinary Poisson-binomial algorithms

### Convolve

`dpb_convolve` is a direct translation of the Bernoulli polynomial/dynamic
convolution.

### DivideFFT

Upstream uses FFTW and switches to a grouped divide-and-conquer FFT only for
large inputs. The Fortran version removes the external FFTW dependency and
implements an in-package radix-2 Cooley-Tukey FFT. It recursively divides the
probability vector and merges probability polynomials with FFT convolution.
Small subproblems use direct convolution.

This computes the same probability polynomial; the grouping threshold and
performance profile are not intended to reproduce FFTW exactly.

### Characteristic

The characteristic-function formula is retained. The upstream FFTW forward
transform is replaced by direct DFT inversion, which avoids an external
library and is straightforward to validate. This has O(n^2) transform cost
instead of FFTW's O(n log n) cost.

### Recursive

The upstream two-column recursive dynamic program is translated directly,
including its reduced-memory structure.

### Approximations

Arithmetic-mean binomial, geometric-mean binomial, geometric counter-mean
binomial, Poisson, normal and refined-normal formulas follow the upstream
parameterizations and continuity corrections. The package is self-contained:
regularized beta/gamma functions used by binomial/Poisson calculations are
implemented in Fortran.

## Generalized distribution

For each trial, the implementation first converts the two possible integer
outcomes to a lower value plus a nonnegative jump. If `val_p` is the lower
outcome, its probability is complemented, matching the upstream
transformation. Deterministic trials are removed and folded into the support
shift.

The exact generalized methods preserve the upstream GCD compression of jump
sizes. `Convolve` uses direct sparse-step convolution; `DivideFFT` recursively
merges sparse Bernoulli polynomials using the native FFT; `Characteristic`
uses characteristic-function DFT inversion. Normal and refined-normal methods
use the weighted mean, variance, skewness correction and continuity correction
from the upstream implementation.

## Multiplicity weights

Weights are expanded internally exactly as the R wrappers intend. The Fortran
`qpbinom` determines its support from the expanded weighted distribution. This
also avoids an apparent upstream wrapper inconsistency in which the ordinary R
`qpbinom` computes its final support bounds from unexpanded `probs` when
`wts` is supplied.

## Fortran API choices

Fortran does not have R's `NULL` argument convention. Therefore complete
ordinary distributions are returned by `dpbinom`/`ppbinom`, while scalar and
vector observations use explicit `_at` and `_values` routines.

Generalized complete distributions use:

```fortran
type(gpb_table)
  integer :: lower
  integer :: upper
  real(dp), allocatable :: values(:)
end type
```

This makes negative support values safe and explicit. Outcome `x` maps to
`values(x-lower+1)`.

## Validation

The test suite includes:

1. fixed independent PMF references for the ordinary distribution;
2. agreement of all four exact ordinary methods;
3. fixed independent references for all six ordinary approximations;
4. a 96-trial FFT/direct-convolution stress comparison;
5. multiplicity weights, deterministic probabilities and quantiles;
6. both ordinary RNG strategies;
7. fixed enumerated references for a generalized law with negative outcomes;
8. agreement of all three exact generalized methods;
9. generalized GCD compression and weighted multiplicities;
10. independent normal/refined-normal generalized reference values;
11. generalized-to-ordinary reduction; and
12. both generalized RNG strategies.

The release build uses GNU Fortran 14.2 with Fortran 2018, optimization,
warnings promoted to errors, and runtime/bounds checking.
