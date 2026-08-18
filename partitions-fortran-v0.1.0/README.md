# partitions-fortran

Modern free-format Fortran 2018 translation of the computational code in the
R package `partitions` 1.10-9.

The library enumerates integer partitions, distinct and restricted
partitions, bounded/block partitions, compositions, set partitions,
permutations, multiset arrangements, and riffle shuffles.  It also provides
exact 64-bit partition/count functions, conjugates, Durfee sizes, and binary
composition conversions.

## Build with FPM

```text
fpm build
fpm test
fpm run --example demo
```

All source is free-format `.f90`; implicit typing and implicit external
interfaces are disabled in `fpm.toml`.

## Small example

```fortran
program example
    use partitions
    implicit none
    integer, allocatable :: a(:,:)

    print *, p(100)              ! 190569292
    print *, q(100)              ! 444793
    print *, r(5,12)             ! 13

    a = parts(5)
    print *, size(a,2)           ! 7
end program example
```

The package also exposes explicit names such as `partition_count`,
`distinct_parts`, `restricted_parts`, `set_partitions`, and
`multiset_permutations`.  See `TRANSLATION_NOTES.md` for the full mapping and
for representation differences from R.

## License and attribution

The upstream package declares `License: GPL`; this translation is distributed
under the same license terms.  Original metadata, citation, and computational
sources are retained in `upstream/`.
