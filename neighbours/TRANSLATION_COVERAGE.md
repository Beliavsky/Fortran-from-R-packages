# Translation coverage

Source package: `neighbours` 0.1-5 (2025-12-17).

## Directly translated

### `R/neighbourfun.R`

- numeric neighbourhood with scalar-budget/zero-sum behavior
- numeric neighbourhood with `sum=FALSE`
- numeric neighbourhood with two-element budget interval
- random/fixed step size
- scalar/vector bounds
- active-coordinate selection
- `update="Ax"`
- permutation neighbourhood
- logical unconstrained toggling
- logical fixed-cardinality exchange
- logical bounded-cardinality branch
- `type="5/10/40"`
- computational output of `compare_vectors`
- internal `random_vector` for logical/numeric solutions

### `R/next_subset.R`

- NEXKSB next-combination algorithm

## Representation differences

R `neighbourfun()` returns closures whose bodies are constructed using
`substitute`. Modern Fortran has no equivalent need for runtime expression
rewriting, so configuration derived types hold the captured values and explicit
procedures apply the moves.

R stores updated `A*x` as an attribute on vector `x`. Fortran arrays do not have
arbitrary R-style attributes, so `numeric_neighbour` accepts an optional `ax`
argument and updates it in place.

R `next_subset` returns `NULL` after the final subset. The Fortran routine keeps
the final subset in `a` and returns `has_next=.false.`.

`compare_vectors` is primarily a console-formatting helper. The Fortran routine
returns consecutive Hamming distances and an optional difference mask; R's
printing/message formatting is intentionally not reproduced.

R permutation neighbourhoods can operate on arbitrary lists. Fortran generic
procedures cover the practical atomic cases: real, integer, logical and
character arrays. Heterogeneous R lists have no direct Fortran analogue.

## Source behaviors intentionally retained

The upstream `kmin < kmax` logical branch ignores `active`; the translation does
the same. The upstream internal numeric `random_vector` sets selected entries to
zero even if zero is outside the requested numeric interval; that behavior is
also retained.

For active sets in other branches, the Fortran API uses actual 1-based indices.
This represents the documented meaning of `active` directly and avoids R's
subsetting-expression implementation detail for non-contiguous numeric active
vectors.

## NMOF and quadprog

Both are only `Suggests` dependencies in the R package; neither is called by the
package computational source. NMOF is used in the vignette to demonstrate local
search and portfolio routines, while quadprog appears only in a vignette
comparison. Therefore neither is a mandatory core dependency.

The supplied NMOF Fortran translation is vendored for provenance/integration,
and `integration/nmof-demo` demonstrates direct callback interoperability.
quadprog is not needed to build or use this translation.
