# spantest modern Fortran

A self-contained modern Fortran translation of the computational code in the
R package `spantest` 1.4-0. The library implements classical, Monte Carlo, and
high-dimensional mean-variance spanning tests and the package's simulation
engine. R objects, Rcpp bindings, documentation-generation hooks, and other
R-runtime infrastructure are intentionally omitted.

## Implemented procedures

- `span_bj`: Britten-Jones tangency-portfolio spanning test
- `span_f1`: Kan-Zhou F1 alpha-spanning test
- `span_f2`: Kan-Zhou F2 variance-spanning test
- `span_grs`: Gibbons-Ross-Shanken alpha test
- `span_hk`: Huberman-Kandel joint mean-variance test
- `span_km`: Kempf-Memmel GMVP spanning test
- `span_py`: Pesaran-Yamagata alpha test
- `span_gl_a`: Gungor-Luger alpha-only sign-flip Monte Carlo test
- `span_gl_ad`: Gungor-Luger joint sign-flip Monte Carlo test
- `span_as`: Ardia-Sessinou subseries Cauchy-combination test
- `span_simulate`: twelve normal, Student-t, and skew-t IID/AR/GARCH DGPs

All public procedures are re-exported by `use spantest`. Results use typed
Fortran structures with a numeric status and diagnostic message.

## Build

With FPM:

```text
fpm test
fpm run --example spantest_demo
```

With GNU Fortran and Make:

```text
make checked
make optimized
```

On Windows, run `tools\run_tests.bat checked` or
`tools\run_tests.bat optimized` from a command prompt with `gfortran` on PATH.

## Minimal example

```fortran
program demo
  use spantest
  implicit none
  type(simulation_result) :: sim
  type(span_result) :: ans

  sim = span_simulate(250, 3, 8, ncp=0.2_dp, dgp=1, seed=123)
  ans = span_grs(sim%r1, sim%r2)
  write(*,'(a,es14.6)') 'GRS p-value: ', ans%pval
end program demo
```

## Numerical notes

The formulas and finite-sample degrees of freedom follow the upstream package.
Linear systems, probability functions, random-number generation, and
Student-t/skew-t simulation are implemented locally to keep the package
self-contained. Seeded runs are reproducible within this Fortran library, but
Monte Carlo draws are not bit-for-bit identical to R's Mersenne-Twister stream.
See `docs/PORTING_NOTES.md` for details.

## License

GPL-3.0-only, matching the upstream package. The complete upstream source and
original archive are retained under `upstream/` for provenance.
