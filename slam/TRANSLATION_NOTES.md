# Translation notes

## Upstream

- R package: `slam`
- Version translated: 0.1-56
- Upstream authors: Kurt Hornik, David Meyer, Christian Buchta
- Upstream license declaration: GPL-2

The original package DESCRIPTION is retained as `ORIGINAL_DESCRIPTION`.

## Design choices

The original code mixes high-level R methods and C kernels. The Fortran port
folds those layers into ordinary modules and derived types. The public umbrella
module is `slam`; applications normally need only `use slam`.

`dp` is defined as `kind(1.0d0)`. Sparse indices are ordinary Fortran integers;
64-bit integers are used internally for linearized positions and products of
dimensions where overflow risk is greater.

The R package permits explicit stored zero values in some malformed or internal
objects. Constructors in the Fortran API validate coordinates and uniqueness
but do not automatically reject explicit zero values, because several original
algorithms and tests rely on being able to sanitize such representations.
Conversion from dense form stores nonzero values and NaNs only.

## Validation performed

The project was compiled using GNU Fortran 14.2.0 with Fortran 2018 mode,
`-Wall -Wextra -Wconversion -fcheck=all -fbacktrace -O0`. The included tests and
example both pass under that build.

FPM itself was not installed in the execution environment. The FPM manifest is
included and has no external dependencies.
