# Porting notes

## Scope

The upstream computational implementation is the single R file
`R/spIndexTrack.R`. It contains the main continuation/MM solver, four MM update
rules, and a KKT projection routine. All of those numerical components are
translated.

The bundled `INDEX_2010` R data object and presentation-oriented `xts` examples
are preserved in the upstream snapshot but are not converted into compiled
Fortran data.

## Data representation

R matrices become rank-2 `real(dp)` arrays with observations in rows and assets
in columns. Benchmark returns and weights are rank-1 arrays. Allocatable output
arrays replace R vectors returned from a function.

## Projection

The upstream sorted bisection routine solves

```text
minimize ||w-z||^2
subject to sum(w)=1 and 0 <= w_i <= u.
```

The translated `bisection` solves the same scalar KKT equation directly by a
monotone bisection on the multiplier. `project_capped_simplex` exposes the more
natural `z`-based interface.

## Largest eigenvalues

Every matrix passed to the upstream symmetric `eigen()` call is a weighted
cross-product and is positive semidefinite. The port therefore uses a
self-contained deterministic power iteration. This changes only roundoff and
convergence details and avoids an external LAPACK requirement.

## Accelerated MM loop

The two successive MM updates, squared iterative acceleration parameter,
projection, objective backtracking, continuation schedule, and relative-change
tolerances follow the R source. The iteration bound is enforced safely across
all continuation stages.

## HDR source defect

R's `&&` operator is scalar, but the upstream HDR objective uses it with vector
conditions. The intended elementwise loss is used by default. The optional
`source_compatible_hdr_objective` flag reproduces the actual R objective for
parity investigations. The MM update itself already uses elementwise indexing
in the upstream source and is translated directly.

## Input validation

The Fortran library makes numerical requirements explicit:

- at least two assets,
- matching observation counts,
- finite values,
- strictly positive `lambda`,
- feasible `n * upper_bound >= 1`,
- positive Huber parameter for HETE/HDR, and
- a feasible supplied initial portfolio.

Errors are returned as status codes and messages rather than R exceptions.
