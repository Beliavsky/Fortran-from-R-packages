# Build report

Validated compiler:

```text
GNU Fortran (Debian 14.2.0-19) 14.2.0
```

Validated commands:

```sh
make clean
make check
make release
make app example
./build/check/app/waveslim_demo
./build/check/example/image_packet_example
```

Both checked and optimized suites pass all seven test programs with no compiler
warnings. All project Fortran source lines are at most 132 characters, so no
nonstandard free-line-length option is required.

An FPM executable was not installed in the validation environment. The
`fpm.toml` file was parsed as TOML, the same source graph was built through the
included Makefile, and source roots were checked for duplicate module and
program names.
