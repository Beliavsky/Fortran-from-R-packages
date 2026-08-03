# Porting notes

## R objects

The R functions return unnamed lists and print values as a side effect. The
Fortran procedures return a typed `jdmbs_result` and perform no printing.

## Plotting

The `ggplot2`, `png`, graphics, and `igraph` display code is omitted. Complete
paths can be requested with `control%store_paths = .true.` and plotted by the
calling application.

## Random numbers

R's session RNG is replaced with a deterministic Park-Miller generator and a
Box-Muller normal generator. Supplying the same 64-bit seed reproduces the same
Fortran result on supported compilers; it does not reproduce R's random stream.

## Time scaling

`legacy_mode = .true.` preserves the upstream formulas, including their missing
`sqrt(dt)` Brownian scaling and off-by-one time handling. Corrected mode uses a
standard recursive daily GBM discretization and applies all jumps in an
interval.

## Jump interpretation

The R code rescales `lambda` by the contract length. Consequently the supplied
value is the expected number of events over the full horizon, rather than a
per-year intensity. The Fortran translation makes this interpretation explicit.

The shared-company model's `correlation_matrix` controls jump propagation only.
It does not correlate Brownian shocks and is therefore described in the API as
a jump-transmission matrix.

## Numerical safety

Very large or very negative exponential arguments are clipped before calling
`exp`, avoiding floating-point overflow/underflow traps. A numerical-warning
status and count are returned when clipping occurs.
