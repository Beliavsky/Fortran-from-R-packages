# Porting notes

## Array and missing-value conventions

Fortran arrays are one-based, matching R indexing. Missing results use IEEE
quiet NaNs. `na_rm=.false.` propagates missing values within a window;
`na_rm=.true.` ignores them when a valid statistic can still be calculated.

## Window conventions

A scalar width uses the trailing window ending at each observation. A width
vector with the same length as the input supplies an adaptive trailing width.
A shorter vector is treated as endpoints; consecutive endpoints define
inclusive intervals, matching the examples in the upstream documentation.

## Numerical implementation

The upstream vectorized rolling implementation is replaced with direct calls
to the scalar estimators on each requested window. This favors clarity and
cross-compiler consistency. It produces the same formulas but is not intended
to match the highly optimized `data.table` throughput for enormous R tables.

The simulator uses an internal xorshift generator and Box-Muller normals, so a
Fortran seed is reproducible across supported compilers but does not reproduce
R's random-number stream.

## Portability correction

The upstream algorithm itself is unchanged. The Fortran simulator avoids
logical expressions that assume short-circuit evaluation; this prevents an
out-of-bounds access when a period contains no observed trades.
