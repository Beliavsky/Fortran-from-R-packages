# Build report

Validated on GNU Fortran 14.2.0 with BLAS/LAPACK.

- Checked flags: `-std=f2018 -O0 -g -Wall -Wextra -Wimplicit-interface -Werror -fcheck=all -fbacktrace`
- Optimized flags: `-std=f2018 -O2 -Wall -Wextra -Wimplicit-interface -Werror -fbacktrace`
- Five test programs passed in both configurations.
- The demonstration program passed in both configurations.
- The FPM manifest was parsed as TOML; FPM itself was unavailable.

GNU `ld` emits an executable-stack warning for internal-procedure optimizer
callbacks. See `README.md` for the portability note.
