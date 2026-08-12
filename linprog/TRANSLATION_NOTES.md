# Translation notes

## Source version

- R package: `linprog` 0.9-6
- Package date: 2026-01-19
- License: GPL (>= 2)
- Dependency: user-supplied `lpSolve-fortran` 0.1.0

## Representation changes

R lists/data frames and row/column names are represented by explicit Fortran
derived types and numeric arrays.  `linprog_result` exposes the numerical
components of the original `solveLP` object directly:

- `solution`
- `opt`
- `iter1`, `iter2`
- `basis`, `basic_values`
- `allvar_opt`, `allvar_cvec`, `allvar_min_c`, `allvar_max_c`
- `allvar_marg`, `allvar_marg_reg`
- `con_actual`, `con_bvec`, `con_free`, `con_dual`, `con_dual_reg`
- final `tableau`

The lpSolve-backed path intentionally remains less detailed, matching the R
package, which only returns the solution/objective and constraint activity for
that branch.

## Numerical fidelity

The internal simplex follows `R/linprog.R`, including:

- slack-variable construction;
- two-phase initialization when the zero solution is infeasible;
- artificial objective row;
- the package's objective-decrease pivot-column heuristic;
- coefficient zeroing in Phase II;
- detailed sensitivity formulas;
- explicit dual construction;
- historical equality-constraint behavior;
- the literal `99` / `77` minimization sensitivity placeholders.

Fortran does not guarantee short-circuit evaluation of `.and.` or `.or.`.
Index/ratio guards are therefore implemented as explicit nested tests rather
than literal translations of R logical expressions.

## MPS

`readMps` follows the deliberately limited R parser.  Equality rows and
`LO`/`FX`/`FR` bounds are rejected.  `G` rows are multiplied by -1 and returned
as ordinary `<=` constraints.  `UP` bounds are converted to extra rows.

`writeMps` sanitizes labels to unique 8-character MPS names and writes an MPS
file that round-trips through the translated reader.
