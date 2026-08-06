# Build report

Validation environment:

- GNU Fortran 14.2.0
- GNU Make
- Fortran standard mode: `-std=f2018 -pedantic`

Validated configurations:

- Checked: `-O0 -g -Wall -Wextra -Werror -fcheck=all -fbacktrace`
- Optimized: `-O3 -Wall -Wextra -Werror`

GNU Fortran 14.2 emits optimizer-only false positives for allocatable descriptors
inside the vendored robust-covariance module. The optimized target disables only
`-Wuninitialized` and `-Wmaybe-uninitialized`; the checked target retains all
warnings and runtime checks. All other diagnostics remain errors.

Both configurations passed six test executables. The demonstration executable
also ran successfully. FPM was not installed in the validation environment;
`fpm.toml` was parsed as TOML and the same ordered source graph was built through
GNU Make.
