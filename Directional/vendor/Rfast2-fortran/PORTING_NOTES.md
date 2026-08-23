# Porting notes

## API design

Rfast2 exposes R functions, S3 methods, Rcpp objects, and C++ kernels. The
Fortran API is array based. Dotted R names are represented with underscores,
and structured R lists are represented by derived result types.

The top-level module is:

```fortran
use rfast2
```

The supplied Rfast translation is an FPM path dependency rather than copied
into the Rfast2 source modules.

## Deliberate corrections

Several upstream defects were not reproduced:

1. **Cauchy RNG**: upstream `Random.h` applies `tan(pi*z)` to a normal draw.
   The Fortran code uses the inverse Cauchy transform of a uniform draw.
2. **Integer sampling with replacement**: the C++ modulus is based on the
   upper bound rather than the number of admissible values and can exclude an
   endpoint. The Fortran sampler uses all `1:n` values.
3. **Sampling without replacement**: the upstream implementation can allocate
   only `size` output slots but initialize `n` entries. The Fortran version
   uses a bounded partial Fisher-Yates shuffle.
4. **Gamma shape < 1**: the boosting transformation uses an independent
   uniform draw in the Fortran implementation.
5. **`col.waldpoisrat`**: the R source contains `Rfast::colmeans/n2` and also
   divides already computed column means by sample sizes. The scalar/vector
   Fortran implementation follows the intended Poisson-rate Wald formula.
6. **`gammapois.mle`**: an upstream R branch refers to an undefined `ea`.
   The Fortran routine optimizes the intended gamma-Poisson likelihood in
   positive log-parameters.
7. **Kaplan-Meier port safety**: tied-time loops are written with explicit
   bound checks because Fortran does not guarantee short-circuit evaluation.
8. **Zero-truncated Poisson regression**: `log(exp(mu)-1)` is evaluated with a
   stable large-`mu` branch.

## RNG compatibility

The PCG32 recurrence is implemented in portable Fortran integer/limb
arithmetic. Normal generation in this Rfast2 layer uses Box-Muller rather than
the external R `zigg` interface, so normal streams are distribution-compatible
but not bit-for-bit compatible with upstream ziggurat streams. The vendored
Rfast dependency retains its own translated zigg code.

## R-specific features

Formula parsing, data frames, environments/hash objects, S3 printing, R
parallel-cluster orchestration, and package benchmarking/attachment messages
are intentionally not reproduced. Computational kernels are exposed directly.

## Vendored Rfast warning cleanup

The active vendored Rfast source has a small set of exact real equality/inequality checks rewritten with ordering/absolute-value equivalents so the complete dependency graph can be built with GNU Fortran's `-Wcompare-reals -Werror`. These are semantics-preserving source cleanups only. The exact user-supplied Rfast archive is retained unchanged at
`upstream/dependencies/Rfast-fortran-v0.3.0.zip`.
