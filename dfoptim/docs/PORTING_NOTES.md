# Porting notes

## Modified Nelder-Mead

The package intentionally uses the cross-product surrogate

```text
sgrad = V_difference * function_difference
```

rather than solving the simplex linear system for a conventional simplex
gradient. This behavior was restored deliberately in the R package in 2020
and is preserved here.

The original nonstandard shrink/restart formulas, regular-simplex
initialization, convergence criteria, and restart test are also retained.

## Bound transformation

`nmkb` uses the same transformations as the R implementation:

- two finite bounds: hyperbolic tangent transformation
- finite lower bound: exponential transformation
- finite upper bound: reflected exponential transformation
- no bounds: identity

The R documentation says transformed starting values must be strictly inside
the bounds, although the source initially performs an inclusive check and can
then evaluate `atanh(1)` or `log(0)`. The Fortran port enforces the documented
strict requirement and returns `dfoptim_invalid_input` instead of generating
an infinity.

## Randomized Hooke-Jeeves ordering

The R implementation calls `sample.int` for every exploratory move. The
Fortran port preserves randomized coordinate order but uses a package-local
xorshift generator. Consequently, a numeric seed does not reproduce R's RNG
sequence, but it is portable and reproducible within this package.

## Bounded MADS correction

The supplied R code transforms a bounded problem to `[-1,1]^n`, but then tests
poll and line-search points against the original physical bounds. It also uses
the physical `scale` directly in normalized coordinates. For general bounds,
this can accept physical points outside the requested box or reject valid
points.

The Fortran port checks normalized candidates against `[-1,1]^n` and converts
physical scales with

```text
normalized_scale = 2 * physical_scale / (upper - lower)
```

before constructing polls. This preserves the intended box constraints and
the documented interpretation of `scale`.

## Non-finite objective values

MADS treats NaN or infinity as a strict-barrier rejection. Other methods also
replace non-finite trial values by a very large objective value. A non-finite
starting objective is reported as an error. Every actual callback invocation
is counted as a function evaluation.

## Omitted code

The source package has no plotting or compiled native code. R-specific list,
data-frame, warning, and dynamic function-dispatch behavior is not reproduced.
