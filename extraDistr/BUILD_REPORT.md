# Build report

## Environment

- GNU Fortran 14.2.0
- GNU Make
- Linux x86-64
- Fortran standard mode: Fortran 2018
- FPM executable: not installed in the validation environment

## Checked configuration

```text
-O0 -g -std=f2018 -ffree-line-length-none -fimplicit-none
-Wall -Wextra -Werror -Wno-compare-reals -Wno-integer-division
-fcheck=all -fbacktrace
```

Result: all five test programs passed.

## Optimized configuration

```text
-O3 -std=f2018 -ffree-line-length-none -fimplicit-none
-Wall -Wextra -Werror -Wno-compare-reals -Wno-integer-division
```

Result: all five test programs passed without warnings.

## Additional checks

- Demonstration program compiled and ran successfully.
- `fpm.toml` parsed successfully as TOML.
- All 194 upstream NAMESPACE exports were found in the translated Fortran
  sources.
- All translated source, build, test, example, and documentation files are
  ASCII-only; the upstream snapshot is preserved unchanged.
- No external numerical library is required.
