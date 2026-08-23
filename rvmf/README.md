# rvMF-fortran

Modern Fortran translation of the computational core of the R package **rvMF 0.1.2**.

The upstream package implements the non-rejection von Mises--Fisher sampler of Kang and Oh (2024). This port translates the core algorithm, including the discrete mixture of symmetric beta distributions and the condensed base-64 lookup-table sampler used by upstream `rvMF64.cpp`.

## API

- `rvmf_sample(x, mu, kappa)` -- vMF random vectors; each row of `x` is a sample.
- `rvmf_angle_sample(w, p, kappa)` -- samples the inner product with the mean direction.
- `dvmf_angle(r, p, kappa)` -- density of that inner product.
- `log_chf(kappa, d1)` -- log confluent hypergeometric normalization `log M(d1/2,d1,2*kappa)`.

The implementation is self-contained. The attached Rfast Fortran translation was reviewed because upstream `rvMF()` uses `Rfast::matrnorm`; this port uses a local standard-normal generator instead, avoiding a large runtime dependency for one elementary operation.

## Build

```sh
fpm test
```

or, with gfortran directly:

```sh
gfortran -std=f2018 -Wall -Wextra -fcheck=all src/rvmf.f90 test/test_rvmf.f90 -o test_rvmf
./test_rvmf
```

## License

GPL-3.0-or-later, following the upstream rvMF package. See `LICENSE` and `NOTICE.md`.
