# matlab-fortran

A modern Fortran 2018 translation of the non-graphical computational parts of
the R package `matlab` 1.0.4.1.

The upstream package emulates commonly used MATLAB utility functions in R. This
port provides the corresponding numerical, matrix, array, number-theory, path,
and timing routines through an FPM-ready Fortran library.

## Implemented functionality

- Numeric helpers: `ceil`, `fix`, `matlab_mod`, `rem`, `nextpow2`, `pow2`,
  `linspace`, `logspace`, standard deviation, and column sums.
- Array helpers: `zeros`, `ones`, `eye`, `find`, `fliplr`, `flipud`, `rot90`,
  `repmat`, `reshape2d`, two- and three-dimensional `meshgrid`, and `padarray`.
- Shape helpers: `shape_of`, `size_dim`, `numel`, `ndims`, and `isempty`.
- Special matrices: Hilbert, Vandermonde, magic, Pascal, and Rosser matrices.
- Number theory: prime generation, primality testing, and integer factorization.
- Path and utility helpers: `fileparts`, `fullfile`, `filesep`, `pathsep`,
  `strcmp`, `tic`, and `toc`.

The umbrella module is named `matlab`.

```fortran
program example
    use matlab, only : dp, linspace, magic
    implicit none

    real(dp), allocatable :: x(:), m(:, :)

    x = linspace(0.0_dp, 1.0_dp, 6)
    m = magic(5)
end program example
```

## Build

With FPM:

```text
fpm test
fpm run --example demo_matlab
```

With GNU Make:

```text
make MODE=debug test
make clean
make MODE=release test
make MODE=release example
```

The checked build enables bounds, allocation, floating-point, and uninitialized
value diagnostics. The release build uses `-O3` with warnings treated as errors.

## Porting scope

Graphics (`imagesc`, `colorbar`, color-map helpers) are omitted. R's
heterogeneous `cell` arrays and S4 dispatch machinery do not have direct
Fortran equivalents and are also omitted. Rank-generic inquiries are provided,
while allocation-producing routines focus on the rank-1/rank-2 cases most
useful in Fortran; `meshgrid3` supplies the package's rank-3 grid operation.

See `docs/API_MAP.md` and `docs/PORTING_NOTES.md` for exact mappings and known
source-compatibility details.

## License

The upstream package is licensed under the Artistic License 2.0. This
translation is distributed under the same license. Original copyright and
provenance information are retained in `COPYRIGHTS`, `LICENSE`, and `upstream/`.
