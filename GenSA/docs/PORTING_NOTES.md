# Porting notes

## Directly preserved

- Generalized visiting-distribution formula
- Temperature schedule
- Generalized acceptance law
- Full-vector followed by coordinatewise proposals
- Periodic wrapping at finite bounds
- Temperature restart
- `ran2` uniform generator and Marsaglia polar normal generator
- Core controls and stopping conditions
- Hybrid global/local search design

## Adaptations

The upstream C++ engine invokes an embedded L-BFGS-B translation for smooth
functions and R's `constrOptim` for nonsmooth functions. This port replaces
those runtime-specific paths with native Fortran projected BFGS and bounded
coordinate-pattern search. Therefore, even with the same random seed, complete
optimization trajectories need not be bit-for-bit identical.

The upstream trace records selected improvement events. The Fortran trace
records the initial state and every completed outer iteration. Its columns are
the same: step, temperature, current objective, and current best objective.

The R wrapper currently fixes its constraint callback to `NULL`, although the
C++ engine contains constraint support. The Fortran API exposes that capability
as an optional logical procedure.

The hidden R control `high.dim` is always true upstream. In Fortran it is
represented more clearly by `control%local_search`.

## Array and callback conventions

All vectors are rank-one `real(dp)` arrays. Bounds and the optional initial
point must have identical size. Objective and constraint procedures receive an
`intent(in)` vector and may use host or module association for extra data.

## Reproducibility

The `ran2` sequence is tested against fixed reference values. Local-search
adaptations, compiler floating-point transformations, and user objective code
can still produce platform-level differences in final paths.
