# API coverage

Coverage basis: **exported computational R functions**. Package status: **substantial**.

Mapped: **161 of 162 (99.4%)**. The percentage counts functions with meaningful mappings; it is not an assertion of full R/S3, boundary, RNG, optimizer, or option compatibility.

## Major translated families

- Real periodic DWT, adaptive preconditioned wavelets on an interval, and upstream-style recursive stationary `wst` transforms and inverse transforms in 1-D, 2-D, and 3-D, including one-level `conbar` synthesis, recursive `av.basis`, minimum-entropy stationary/packet basis selection, and selected-basis reconstruction.
- Decimated and stationary packet-tree construction, typed 1-D and 2-D packet access/mutation, packet rotations, raw stationary packet addressing, packet feature matrices, regression selection, discrimination feature selection, pooled-covariance LDA fitting, and classification of new stationary packet features.
- Periodic Geronimo/Donovan3 multiple-wavelet decomposition/reconstruction, all upstream pre/postfilter branches, typed multiple-coefficient access/mutation, and multiple-wavelet thresholding.
- Thresholding and cross-validation, including manual/universal/SURE paths, variance-adjusted irregular-grid shrinkage, Ogden-Parzen rules, scalar and level-specific stationary CV, and the periodic empirical-Bayes `BAYES.THR` posterior-median rule.
- Filter selection, support, transform matrices, scaling-function refinement, wavelet/scaling moments, grid interpolation, and ordinary/multiple-wavelet first-last bookkeeping.
- LSW/autocorrelation-wavelet construction, raw local spectra, unsmoothed `ewspec` correction, seeded `LSWsim`/`checkmyews`, Fourier helpers, standard signals, simulations, claw-distribution utilities, and high-resolution density projection with covariance bands.

## Untranslated computational exports

These functions remain in the denominator and are not credited as translated:

- `putC.wp`

## Main remaining work families

- **Irregular transforms:** `irregwd`, `accessc`, and the manual/universal paths of `threshold.irregwd` are translated for the typed interpolation map and real periodic DWT. Symmetric boundaries, complex filters, and additional threshold policies remain limitations rather than separate unmapped exports.
- **Density estimation:** the Daubechies-Lagarias projection kernel, zero-boundary `denwd`, sampled-basis `CWavDE`, and covariance propagation used by `dencvwd` are translated. The covariance implementation favors transparent dense matrix products over the upstream specialized band-storage acceleration.
- **Complex thresholding:** `cthresh` provides the default periodic Lina-Mayrand 3.1 decimated transform, robust noise covariance, multiwavelet-style threshold, and levelwise empirical-Bayes mixture paths. `find.parameters` uses the upstream likelihood and bounds with a deterministic coordinate optimizer. Translation does not yet expose other complex filters, TI transforms, filter averaging, plotting, or R/NAG optimizer objects.
- **API guard:** upstream `putC.wp` is an error-only guard because packet objects do not have separate father-coefficient levels; it remains intentionally uncredited rather than inflating computational coverage.
