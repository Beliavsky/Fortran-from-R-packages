# Build report

## Environment

- GNU Fortran 14.2.0
- Python 3.13 used only for archive/manifest checks
- FPM executable: not installed in the validation container

## Checked build

Command:

```sh
make clean
make check
```

Flags:

```text
-std=f2018 -Wall -Wextra -Werror -pedantic -Wno-compare-reals
-O0 -g -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow
```

Result: all five tests passed; demonstration passed.

## Optimized build

Command:

```sh
make optimized
```

Flags:

```text
-std=f2018 -Wall -Wextra -Werror -pedantic -Wno-compare-reals -O3
```

Result: all five tests passed; demonstration passed.

`-Wno-compare-reals` is required only because the preserved Goldfarb-Idnani
core uses exact comparisons to the algorithmic sentinel values `0` and `1`.
The translated NlcOptim module itself builds with the remaining warnings
promoted to errors.

## Demonstration result

The first upstream documentation example converged to a feasible point with an
objective near `5.395e-2` and maximum nonlinear equality violation below
`1e-6` in both build modes.

## FPM target-name validation

Version 0.1.1 fixes the duplicate module/program units previously present in
`app/demo_nlcoptim.f90` and `example/basic_nonlinear.f90`. The example now uses
`basic_nonlinear_problem` and `basic_nonlinear`, while the application retains
`demo_problem` and `demo_nlcoptim`.

The validation container does not provide an FPM executable. The `fpm.toml`
manifest was parsed successfully, all module names under FPM's standard source
roots (`src`, `app`, `example`, and `test`) were checked for uniqueness, and both
the application and example targets were compiled and run independently against
the library sources.
