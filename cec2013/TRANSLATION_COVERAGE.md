# Translation coverage

## Translated

The complete computational surface of the upstream R package is translated:

- R `cec2013(i, x)` validation semantics at the numerical-array level;
- C `test_func` dispatcher for functions 1-28;
- Sphere;
- Ellipsoidal;
- Bent Cigar;
- Discus;
- Different Powers;
- Rosenbrock;
- Schaffer F7;
- Ackley;
- Weierstrass;
- Griewank;
- Rastrigin and noncontinuous Rastrigin;
- Schwefel;
- Katsuura;
- Lunacek bi-Rastrigin;
- expanded Griewank-Rosenbrock;
- expanded Scaffer F6;
- all eight composition functions;
- shift, rotation, asymmetry, oscillation, and composition-weight helpers;
- vector and matrix-point evaluation;
- all dimension-specific rotation data and shared shift data.

## R-only functionality omitted

The following are not computational benchmark algorithms and are omitted:

- `.C` registration and R package loading;
- `R_CheckUserInterrupt`;
- R vector/matrix type inspection;
- R error-object construction.

They are replaced by typed Fortran arrays and integer status codes.

## Deliberately preserved source behaviors

The goal is to translate the supplied package, not replace it with a different
CEC2013 implementation.  Therefore several unusual source behaviors are kept:

1. `dif_powers_func` uses the C integer expression `2 + 4*i/(nx-1)`.
2. `asyfunc` only writes output entries for positive inputs; nonpositive entries
   retain their previous destination values.
3. `grie_rosen_func` computes a rotation but then executes the upstream
   equivalent of `z = y + 1`, discarding the rotation.
4. Data loading reads the first `10*n` whitespace-separated shift values exactly
   like the C implementation, rather than interpreting file line boundaries as
   fixed 100-dimensional records.

These details are covered by cross-language regression tests.

## Representation differences

- Upstream uses global C work arrays; Fortran uses local work arrays.
- Upstream caches global data based on dimension/problem; Fortran stores the
  loaded data explicitly in `cec2013_context`.
- R matrices store evaluation points by row.  The Fortran batch API stores one
  point per column to follow column-major numerical conventions.
- Runtime data files are kept external rather than converted into giant Fortran
  constant arrays.
