# Validation

The v0.2.0 source tree was compiled with GNU Fortran 14.2 using strict runtime
checking, bounds checking, implicit-interface enforcement, and floating-point
traps for invalid operations, division by zero, and overflow.

The trawl regression programs pass:

- `test_deoptim_integration`
- `test_distributions`
- `test_fit_intersection`
- `test_functions`
- `test_simulation`

`test_deoptim_integration` compares the trawl optimizer wrapper with a direct
call to the supplied `deoptim_solve` using identical bounds, control values,
and a fixed seed. The returned best member and objective value are required to
match exactly. It also checks reproducibility through `set_trawl_seed()`.

The three tests shipped with the supplied DEoptim Fortran translation were also
run against the embedded `src/deoptim*.f90` sources and pass:

- `test_initial_and_stop`
- `test_mapping_storage`
- `test_strategies`

The example `demo_trawl` also runs successfully.

FPM was not installed in the validation environment, so the FPM source/test
layout was compiled directly with gfortran rather than through the `fpm`
executable.
