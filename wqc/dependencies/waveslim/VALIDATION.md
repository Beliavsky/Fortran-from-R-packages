# Validation

Seven deterministic test programs are included.

1. `test_transforms_1d`
   - DWT and MODWT reconstruction
   - non-dyadic transforms
   - MRA identities
   - Hilbert-pair reconstruction and phase utilities
2. `test_packets`
   - DWPT inverse identity
   - MODWPT dimensions
   - entropy, cumulative-spectrum, cumulative-periodogram, and portmanteau tests
   - packet-basis selection
3. `test_multidimensional`
   - 2-D and 3-D DWT/MODWT reconstruction
   - MRA reconstruction
   - convolution, shifting, dual-tree, and Hilbert transforms
4. `test_statistics_denoise`
   - wavelet variance/covariance/correlation
   - tapers and periodograms
   - thresholding and 2-D denoising
   - cascade and squared-gain calculations
5. `test_long_memory`
   - FDP/SPP/SPP2 spectra and bandpass integrals
   - hypergeometric evaluation
   - Hosking and wavelet-domain simulation
   - long-memory estimators and adaptive bases
6. `test_extended`
   - 2-D packet reconstruction
   - seeded packet bootstrap
   - variance-change detection
7. `test_errors`
   - invalid filters, levels, dimensions, and arguments

The checked build uses bounds checking, backtraces, warnings as errors, and
floating-point traps. The optimized build uses `-O3` with warnings as errors.
All tests pass under GNU Fortran 14.2.0.

Representative errors from the examples:

- 1-D LA8 DWT reconstruction: approximately `2e-12`
- 2-D Haar packet reconstruction: approximately `2e-15`
- 1-D and 2-D dual-tree reconstruction tests: below `1e-8`

The denoising demonstration reduces deterministic test-signal RMSE from about
`0.387` to `0.103` using SURE soft thresholding.
