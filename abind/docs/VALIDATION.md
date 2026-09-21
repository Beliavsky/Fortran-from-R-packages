# Validation

Validation used GNU Fortran 14.2.0 on Linux.

## Successful direct compiler validation

A strict build used:

```text
gfortran -std=f2018 -Wall -Wextra -Werror -fcheck=all -fbacktrace
```

It compiled the library, both deterministic test executables, and the example.
Both tests passed:

```text
All abind unit tests passed.
All abind metadata tests passed.
```

The example produced a `2 2 2` result containing values 1 through 8 in
column-major order.

A second release-style build used:

```text
gfortran -std=f2018 -O2 -Wall -Wextra -Werror
```

Both tests and the example passed again.

The tests cover vector and mixed-rank binding, inserted dimensions from
fractional/out-of-range `along`, `rev_along`, zero-length arrays, arbitrary
subsetting and drop behavior, singleton-dimension removal, leading/trailing
`acorn` selections, name-driven multidimensional `afill` replication,
hierarchical join labels, selected labels, and synthetic `acorn` labels.

## FPM command attempts

The required literal FPM commands were attempted from the package root:

```text
fpm build
fpm test
fpm run --example basic
fpm clean --all
```

Each returned exit status 127 because this sandbox does not contain an `fpm`
executable (`fpm: command not found`). An attempt to obtain the official Linux
FPM release binary was also blocked by the sandbox download restrictions.
Accordingly, this document does not claim an FPM execution pass. The manifest
was parsed as TOML and the same source/test/example targets were compiled and
run directly with GNU Fortran as described above.

After the failed `fpm clean --all` attempt, generated compiler products from
manual validation were removed explicitly before packaging.
