# longmemo-modern-fortran

A standalone modern Fortran translation of the computational parts of the R package `longmemo` 1.1-4, which implements methods from Jan Beran's *Statistics for Long-Memory Processes* and related work.

Plotting code was intentionally omitted. The project does not require R, BLAS, LAPACK, FFTW, or another external numerical library.

## Implemented functionality

- Fractional Gaussian noise autocovariances and simulation
- Fractional ARIMA(0,d,0) autocovariances and simulation
- Paxson's FFT-based fGn simulation
- Fractional Gaussian noise spectral density
- Fractional ARIMA(p,d,q) spectral density
- Paxson approximation and direct summation for `B(lambda,H)`
- Periodograms and Fourier frequencies
- Whittle estimation for fGn and fARIMA(p,d,q)
- Beran goodness-of-fit statistic
- Numerical covariance matrices for Whittle estimates
- Polynomial FEXP estimation using a Gamma GLM with log link
- The five package datasets converted to CSV

## Source layout

```text
src/longmemo_kinds.f90       floating-point kinds and constants
src/longmemo_fft.f90         radix-2 and Bluestein FFT implementation
src/longmemo_linalg.f90      Cholesky solves and SPD matrix inversion
src/longmemo_stats.f90       distributions and random-number helpers
src/longmemo_optimize.f90    scalar minimization and Nelder-Mead
src/longmemo_io.f90          reader for bundled index/value CSV files
src/longmemo.f90             public long-memory API
examples/demo_longmemo.f90   simulation and estimation example
examples/fit_csv.f90         estimate a model from a CSV series
tests/test_longmemo.f90      numerical regression tests
data/*.csv                   datasets from the R package
```

All real calculations use

```fortran
integer, parameter :: dp = kind(1.0d0)
```

and all modules use `implicit none`.

## Build with Make

```sh
make
make test
make example
```

To fit a bundled series:

```sh
make build/fit_csv
./build/fit_csv data/NileMin.csv
```

The Nile example gives approximately:

```text
n = 663
Whittle H = 0.837424
Whittle SE = 0.026138
FEXP H = 0.871798
FEXP order = 1
```

The original R package reports FEXP `H = 0.87178` and selected order 1 for this series.

## Build with fpm

```sh
fpm test
fpm run demo_longmemo
fpm run fit_csv -- data/NileMin.csv
```

## Small API example

```fortran
program example
    use longmemo_kinds, only : dp
    use longmemo_stats, only : set_random_seed
    use longmemo
    implicit none

    real(dp), allocatable :: x(:)
    type(whittle_result) :: fit

    call set_random_seed(12345)
    call sim_fgn_fft(1024, 0.75_dp, x)
    call whittle_estimate(x, "fGn", fit, start_eta=[0.6_dp])

    print *, fit%eta(1)
    print *, fit%std_error(1)
end program example
```

## Compatibility notes

- R names are translated to lower-case snake case. See `API_MAP.md`.
- The slow `B.specFGN(..., k.approx=NA)` path is selected with `k_approx <= 0`.
- Fortran random streams are not identical to R's random streams, so simulated paths differ even with the same integer seed.
- FEXP uses an orthonormal polynomial basis spanning the same polynomial space as R's `poly()`. Fitted values, H estimates, and model order closely reproduce the R implementation, but the scaling and signs of individual nuisance-polynomial coefficients can differ.
- The FFT implementation uses radix-2 Cooley-Tukey for power-of-two lengths and Bluestein's algorithm otherwise.
- Invalid arguments currently terminate with `error stop`; callers that need recoverable error handling can wrap or adapt the validation layer.
- Whittle fARIMA optimization follows the R package's unconstrained approach and uses Nelder-Mead. It does not impose AR/MA stationarity or invertibility constraints.

## Verification

`make test` checks FFT round trips, autocovariances, spectra, periodograms, simulation, Whittle estimation, FEXP estimation, and saved R regression values. In particular:

```text
CetaFGN(H=0.7, m=256)       0.4862165912
CetaARIMA(0,0, m=256)       0.7163356699
CetaARIMA(1,1)[1,1], m=256 14.48898898
```

These agree with the saved tests in the attached R package.

## License and attribution

The original package is licensed under GPL version 2 or later. This translation is a derivative work and is distributed under the same terms. Copies of GPL-2 and GPL-3 are included.

Original authors and contributors credited by the package:

- Jan Beran
- Martin Maechler
- Brandon Whitcher

Dataset-specific sources and literature references are preserved in
`DATA_SOURCES.md`.
