# Translation coverage

## Translated

- All 25 CEC-2005 benchmark objective functions from the supplied C code.
- Basic Sphere, Schwefel 1.2, Rosenbrock, Griewank, Ackley, Rastrigin, and
  Weierstrass kernels.
- Shift/rotation transformations.
- Hybrid-composition weighting and normalization.
- Non-continuous Rastrigin/Scaffer transformations.
- F5 and F12 specialized data matrices.
- Dimension-specific rotations for dimensions 2, 10, 30, and 50.
- Stochastic noise in F4, F17, F24, and F25.
- Scalar and matrix/batch evaluation.
- Official benchmark data files and upstream test/reference data.

## R-only pieces omitted

- `.C` registration/marshalling.
- `R_CheckUserInterrupt`.
- R vector/matrix type checks.
- R namespace/package installation machinery.

These do not change the numerical benchmark definitions.

## Deliberate Fortran differences

### Random-number stream

The original C wrapper obtains normal deviates from R (`norm_rand`).  The
Fortran version uses the intrinsic RNG plus a Box-Muller transform.  Thus the
noise law is the same but seeded stochastic trajectories are not expected to
match R bit-for-bit.  With noise disabled, deterministic values match the
supplied upstream reference data to roundoff.

### Cached context

The R/C wrapper allocates and reloads data for every call.  The Fortran API
adds `cec2005_context` so optimization programs can load a function/dimension
once and evaluate it repeatedly.  The convenience `cec2005_eval` routine
retains the one-call behavior.

### Data files

Official data are kept as run-time files rather than compiled into source.
This makes the translation auditable and avoids embedding several megabytes of
numeric constants in Fortran modules.  Applications should deploy the `data/`
directory and pass its path to `ctx%init` when the working directory differs
from the package root.

## Source fidelity details

The translation preserves the supplied C implementation's orientation of
rotation matrices, hybrid weights, integer truncation in non-continuous
rounding, F5/F8/F20 optimum-on-bound modifications, and the identical numeric
formula used by functions 24 and 25.  F24 and F25 differ in the official search
bounds, not in this package's evaluation routine.
