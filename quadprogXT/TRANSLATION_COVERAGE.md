# Translation coverage

## Translated

All computational R code in `quadprogXT` 0.0.6 is represented:

| Upstream R function | Fortran |
|---|---|
| `solveQPXT` | `solve_qp_xt` |
| `buildQP` | `build_qp_xt` |
| `convertToCompact` | `convert_to_compact` |
| `normalizeConstraints` | `normalize_constraints` |

The exact auxiliary-variable construction for absolute values is preserved,
including:

- optional `n` auxiliary variables for `abs(b)`;
- optional `n` auxiliary variables for `abs(b-b0)`;
- the source `MAP` transformation from positive/negative coordinates;
- slack constraints;
- linear penalties in positive/negative coordinates;
- the small `tol` diagonal regularization on auxiliary variables;
- factorized expanded matrices;
- normalized constraints;
- compact constraints.

## R-only infrastructure omitted

There is no R formula/model object layer, `do.call`, list construction, roxygen
runtime behavior, or R error-condition object.  Fortran uses derived types and
explicit status/message fields instead.

## Small defensive differences

The R source relies on downstream `quadprog` errors for several malformed
argument combinations.  The Fortran translation validates dimensions and
missing companion arguments earlier and returns a `qpxt_problem` build error.

A completely unconstrained problem is accepted and solved directly even when
`compact=.true.`.  This is useful in Fortran and avoids the undefined
`max(integer(0))` behavior that the R `convertToCompact` implementation would
encounter on a zero-column constraint matrix.

## Supplied quadprog dependency

The user-supplied `quadprog-fortran` translation is vendored unchanged except
for exact-real comparison syntax in `quadprog_core.f90`.  `.EQ. 0/1` tests were
rewritten as exact `abs(x-c) <= 0` tests solely to compile under
`-Werror=compare-reals`.  No tolerance was introduced.
