# Translation coverage

## Translated

- Stark-Parker BVLS active-set method.
- Cold start from lower bounds.
- Warm start using `key` and `istate`.
- Bound/active working-set state and output.
- Kuhn-Tucker bound release test.
- Short-cycle protection using the last-bound-variable guard.
- Householder QR least-squares solve.
- Rejection of linearly dependent candidate active columns.
- Feasible interpolation when an unconstrained active-set solution crosses a
  bound.
- Roundoff repair of variables at/outside their bounds.
- Original `3*n` major-loop limit.
- R-wrapper-style fitted values, residuals, and deviance.

## Modernization

- Fixed-form Fortran 77 is replaced by free-form Fortran 2018 modules.
- All interfaces are explicit.
- Work arrays are allocated from array shapes rather than exposed to callers.
- Exact-zero tests from the old source are written without `==`/`.eq.` so
  strict `-Werror` builds do not trigger real-comparison diagnostics.
- `bvls_fit` returns a typed result rather than an R list.

## Small deliberate robustness difference

The original routine returns immediately when every variable has identical
lower and upper bounds, before assigning the solution. The R wrapper happened
to initialize `x` to zero, which is wrong when the fixed bounds are nonzero.
The Fortran translation returns `x=lower=upper` and status
`bvls_no_free_variables` for this case.

R printing/S3 assessor methods are not translated; their numerical outputs are
fields of `bvls_result`.
