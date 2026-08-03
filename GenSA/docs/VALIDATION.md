# Validation

The release archive was validated with GNU Fortran 14.2.0 in two modes.

Checked mode:

```text
-std=f2018 -Wall -Wextra -Wpedantic -Wconversion-extra
-O0 -g -fcheck=all -fbacktrace
```

Optimized mode:

```text
-std=f2018 -Wall -Wextra -Wpedantic -Wconversion-extra -O3
```

Both modes build and run:

- six regression programs;
- five examples;
- the `demo_gensa` application.

The FPM executable was not installed in the validation environment. The
`fpm.toml` file was parsed successfully with a TOML 1.0 parser and uses a
numeric version string.

GNU `ld` can print an executable-stack notice when an internal Fortran
procedure is passed as a callback. The examples intentionally use internal
objective procedures for concise demonstrations. All resulting executables
ran successfully in both validation modes.
