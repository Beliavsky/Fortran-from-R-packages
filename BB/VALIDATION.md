# Validation

The package was compiled with GNU Fortran 14.2.0 using strict interface and
runtime checking:

```text
-std=f2018 -Wall -Wextra -Wimplicit-interface -Werror=implicit-interface -fcheck=all
```

The following tests pass:

1. `test_spg` - Rosenbrock minimization with analytic gradient and an exact-start maximization case.
2. `test_sane` - SANE on a two-equation linear system.
3. `test_dfsane` - DF-SANE on the same system.
4. `test_projection` - ports the BB `projectLinear` examples for inequalities, an equality, and box-equivalent linear constraints.
5. `test_drivers` - `BBoptim`, `BBsolve`, and multistart optimization.

The two examples also build and run successfully.

FPM was not installed in the validation container.  The tree follows standard
FPM layout and the explicit test/example targets are declared in `fpm.toml`;
validation was performed by compiling the same source set directly with
`gfortran`.

Warnings from `quadprog_core.f90` about exact floating-point comparisons are
inherited unchanged from the supplied quadprog translation.  The BB algorithms
also intentionally retain several exact-zero denominator guards from the R
source.  No implicit-interface diagnostics occur.
