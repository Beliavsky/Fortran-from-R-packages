# Porting notes

## Why this is a backend interface

The attached `highs` package contains R and Rcpp wrappers around HiGHS. Rewriting only a small optimizer and calling it a translation would remove the defining LP/MIP/QP algorithms and produce incompatible results. The Fortran project therefore wraps the same bundled HiGHS 1.14.0 engine.

## Runtime loading

FPM cannot build a large CMake C++ dependency automatically and then infer its library directory. A normal external `link = ["highs"]` manifest would reproduce the common `cannot find -lhighs` failure. The C loader avoids that problem:

- `fpm build` links only ordinary Fortran and C objects;
- the HiGHS bridge is loaded on first solver use;
- the backend can be built once, moved, or selected with an environment variable;
- missing backend errors are reported at runtime with the operating-system loader message.

## Bundled R-specific HiGHS modifications

The HiGHS snapshot shipped by the R package contains a few `Rcpp.h`, `Rcpp::Rcout`, `Rcpp::stop`, and `Rprintf` references. `bridge/compat/Rcpp.h` supplies only the tiny stream, exception, and formatted-output compatibility surface needed to compile that snapshot outside R. It does not emulate R or Rcpp.

## Sparse matrices

HiGHS consumes 0-based sparse arrays. Public dense/triplet constructors accept normal 1-based Fortran indices, then store canonical 0-based CSC/CSR arrays. Duplicate triplets are summed.

The Hessian constructor stores the upper triangle and uses HiGHS triangular format. HiGHS evaluates the quadratic objective as:

```text
0.5 * transpose(x) * Q * x + transpose(c) * x + offset
```

## Bounds and infinity

Use `highs_default_infinity` before the backend is loaded. `highs_infinity()` returns the backend value when available. Values at or beyond the HiGHS infinity threshold are interpreted as unbounded.

## Statuses

`highs_solution%call_status` is the API call status (`-1`, `0`, or `1` for error, OK, or warning). `model_status` is the optimization result, such as `highs_model_optimal`, `highs_model_infeasible`, or `highs_model_time_limit`.

## Threading

`highs_control%threads` is passed to the HiGHS `threads` option. `parallel` is enabled automatically when more than one thread is requested. Multiple solver objects remain subject to HiGHS' own scheduler and thread-safety rules.

## Licensing

The Fortran frontend and bridge preserve the R wrapper's GPL-2.0-or-later terms. HiGHS itself remains under its MIT license. The original source snapshots and notices are retained.
