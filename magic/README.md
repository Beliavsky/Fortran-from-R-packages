# magic-fortran

A modern Fortran translation of the computational code in the R package
`magic` 1.6-1-1 by Robin K. S. Hankin.

The library creates and investigates magic squares, panmagic squares, magic
hypercubes, Latin squares, incidence tensors, Hadamard matrices, sparse
antimagic squares, and multiplicative magic squares. It also provides an
arbitrary-rank integer tensor type for the package's array-manipulation
operations.

Plotting, R object dispatch, replacement functions, dimnames, and packaged R
data files are not translated. The complete upstream source is retained under
`upstream/magic-master` for provenance.

## Build with FPM

```text
fpm build
fpm test
fpm run
fpm run --example create_magic_square
```

The project has no external library dependencies.

## Minimal example

```fortran
program example
   use magic, only : ik, magic_square_of_order, is_magic
   implicit none
   integer(ik), allocatable :: square(:, :)

   square = magic_square_of_order(7)
   print *, is_magic(square)
end program example
```

## Main modules

- `magic_tensor`: arbitrary-rank integer arrays and transformations
- `magic_square`: square constructors, sums, tests, and products
- `magic_hypercube`: magic cubes and general hypercube tests
- `magic_combinatorics`: Latin/incidence, Hadamard, SAM, and related methods
- `magic`: umbrella module re-exporting the complete public API

See `API.md`, `API_MAP.md`, and `PORTING_NOTES.md` for details.

## License

The upstream package declares `GPL-2`. This translation is distributed under
`GPL-2.0-only`; see `LICENSE`. Upstream authorship and provenance are retained
in `NOTICE.md` and the `upstream` directory.
