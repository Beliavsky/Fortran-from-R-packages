# waveslim-fortran

**Official CRAN title:** Basic Wavelet Routines for One-, Two-, and Three-Dimensional Signal Processing

Modern Fortran 2018 computational port of the R package `waveslim` 1.8.5.
The library provides one-, two-, and three-dimensional wavelet transforms,
wavelet packets, dual-tree and Hilbert transforms, denoising, wavelet
statistics, long-memory models, and simulation utilities in an FPM-ready
package.

## Highlights

- DWT, MODWT, inverse transforms, MRA, phase correction, and boundary handling
- DWPT and MODWPT with packet tests, basis selection, bootstrap, and 2-D packets
- 2-D and 3-D DWT/MODWT transforms and reconstruction
- Dual-tree 1-D/2-D transforms and Hilbert wavelet pairs
- Wavelet variance, covariance, correlation, spectra, tapers, and change tests
- Thresholding and 2-D wavelet denoising
- FDP/SPP/SPP2 spectra, estimation, adaptive bases, and Hosking simulation
- Twenty-one standard wavelet filters plus six Hilbert-pair filter sets
- No required BLAS, LAPACK, FFT, or R runtime dependency

## Build with FPM

```sh
fpm build
fpm test
fpm run
fpm run --example image_packet_example
```

## Build with GNU Make

```sh
make check
make release
make app example
./build/check/app/waveslim_demo
```

On Windows with FPM:

```bat
build_and_test.bat
```

## Minimal example

```fortran
program example_dwt
  use waveslim
  implicit none
  real(dp) :: x(128)
  real(dp), allocatable :: reconstructed(:)
  type(wavelet_transform) :: wt
  integer :: i

  do i = 1, size(x)
    x(i) = sin(0.1_dp*real(i,dp))
  end do

  wt = dwt(x, 'la8', 4)
  if (.not. wt%status%ok()) error stop trim(wt%status%message)
  reconstructed = idwt(wt)
  print *, maxval(abs(reconstructed-x))
end program example_dwt
```

## Design

R lists and S3 classes are represented by typed Fortran derived types. Routine
names use underscores instead of dots, for example `wave.variance` becomes
`wave_variance`. Errors are reported through `status_type` fields where a
result object is returned.

Common finite-filtered descriptive statistics, type-7 quantiles, normal
distribution functions, and chi-square probabilities delegate to the
shared MIT-licensed `rfortran-core` dependency. Waveslim's lag-dependent
autocovariance and cross-correlation conventions remain package-local.

See `API_MAP.md`, `PORTING_NOTES.md`, and `VALIDATION.md` for coverage and
numerical details.

## License

BSD 3-Clause, preserving the upstream license and copyright of Brandon
Whitcher. The complete upstream source and original archive are retained under
`upstream/` for provenance.
