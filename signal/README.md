# signal

Modern free-form Fortran translation of the computational portions of R package
**signal 1.8-1**. The package is intended to live as the top-level `signal/`
directory in `Beliavsky/Fortran-from-R-packages` and builds with the Fortran
Package Manager (FPM).

The translation focuses on numerical signal-processing operations: FIR/IIR
filter design and filtering, windows, frequency and impulse responses,
interpolation/resampling, Savitzky-Golay filters, polynomial utilities,
chirps, and spectrogram calculations. Plotting, printing, interactive display,
R S3 dispatch machinery, attributes, and other R-only presentation behavior are
not reproduced.

## Build and test

```text
fpm build
fpm test
fpm run --example filter_design
fpm run --example resampling
```

No system BLAS, LAPACK, ARPACK, or external FFT installation is required. This
translation is deliberately standalone: the upstream package's single MASS
numerical use (`MASS::ginv` in `sgolay`) is replaced by a package-local pivoted
linear solve, so adding MASS or a linear-algebra package solely for that call
would not improve reuse. FFT-backed upstream operations use a self-contained
DFT/direct-convolution implementation; this avoids a new dependency but is
slower for large transforms.

The source uses one real kind, `dp`, defined once from `iso_fortran_env::real64`
in `signal_kinds` and re-exported by the public `signal` module.

## Translation coverage

**Package status: substantial. Coverage: 58 of 58 (100%).**

The fraction measures mapped **exported computational R functions**, not complete
R compatibility. `freqs_plot`, `freqz_plot`, and `zplane` are excluded from the
denominator because they are plotting/presentation-only exports. Internal R
helpers are also excluded. A mapped function may still be marked `partial` or
`substantial` where R S3 behavior, matrix/recycling semantics, exact backend
algorithms, warnings, or acceleration paths differ. See
[`docs/API_COVERAGE.md`](docs/API_COVERAGE.md) for material compatibility notes.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
|---|---|---|---|
| `bartlett` | `signal_windows` | `bartlett_window` | complete |
| `bilinear` | `signal_iir_design` | `bilinear_transform` | substantial |
| `blackman` | `signal_windows` | `blackman_window` | complete |
| `boxcar` | `signal_windows` | `boxcar_window` | complete |
| `butter` | `signal_iir_design` | `butter_filter` | substantial |
| `buttord` | `signal_iir_design` | `butter_order` | substantial |
| `chebwin` | `signal_windows` | `chebyshev_window` | substantial |
| `cheb1ord` | `signal_iir_design` | `cheby1_order` | substantial |
| `cheby1` | `signal_iir_design` | `cheby1_filter` | substantial |
| `cheby2` | `signal_iir_design` | `cheby2_filter` | substantial |
| `chirp` | `signal_analysis` | `chirp_signal` | substantial |
| `decimate` | `signal_interpolation` | `decimate_signal` | substantial |
| `ellip` | `signal_iir_design` | `ellip_filter` | substantial |
| `ellipord` | `signal_iir_design` | `ellip_order` | substantial |
| `FftFilter` | `signal_filters` | `make_fft_filter` | partial |
| `filter` | `signal_filters` | `filter_signal` | substantial |
| `fftfilt` | `signal_filters` | `fft_filter_signal` | substantial |
| `MedianFilter` | `signal_filters` | `make_median_filter` | partial |
| `medfilt1` | `signal_filters` | `median_filter_signal` | partial |
| `spencerFilter` | `signal_filters` | `make_spencer_filter` | substantial |
| `spencer` | `signal_filters` | `spencer_smooth` | substantial |
| `FilterOfOrder` | `signal_filters` | `make_filter_order` | partial |
| `an` | `signal_filters` | `unit_phasor_degrees` | complete |
| `roots` | `signal_utils` | `polynomial_roots` | partial |
| `Arma` | `signal_filters` | `make_arma` | partial |
| `as.Arma` | `signal_filters` | `arma_from_zpg`, `make_ma` | partial |
| `Ma` | `signal_filters` | `make_ma` | partial |
| `Zpg` | `signal_filters` | `make_zpg` | partial |
| `as.Zpg` | `signal_filters` | `zpg_from_arma` | partial |
| `polyval` | `signal_utils` | `polyval_real_complex` | substantial |
| `conv` | `signal_utils` | `convolve` | complete |
| `ifft` | `signal_utils` | `inverse_dft` | substantial |
| `filtfilt` | `signal_filters` | `zero_phase_filter` | substantial |
| `fir1` | `signal_fir_design` | `fir1_filter` | substantial |
| `fir2` | `signal_fir_design` | `fir2_filter` | partial |
| `flattopwin` | `signal_windows` | `flattop_window` | substantial |
| `freqs` | `signal_analysis` | `analog_frequency_response` | substantial |
| `freqz` | `signal_analysis` | `digital_frequency_response`, `digital_frequency_response_at` | substantial |
| `gausswin` | `signal_windows` | `gaussian_window` | complete |
| `grpdelay` | `signal_analysis` | `group_delay` | substantial |
| `hamming` | `signal_windows` | `hamming_window` | complete |
| `hanning` | `signal_windows` | `hanning_window` | complete |
| `impz` | `signal_analysis` | `impulse_response` | substantial |
| `interp1` | `signal_interpolation` | `interp1_signal` | partial |
| `interp` | `signal_interpolation` | `interpolate_signal` | substantial |
| `kaiser` | `signal_windows` | `kaiser_window` | substantial |
| `kaiserord` | `signal_fir_design` | `kaiser_order` | substantial |
| `levinson` | `signal_filters` | `levinson_durbin` | substantial |
| `pchip` | `signal_interpolation` | `pchip_interpolate` | substantial |
| `poly` | `signal_utils` | `polynomial_from_roots` | partial |
| `remez` | `signal_fir_design` | `remez_filter` | partial |
| `resample` | `signal_interpolation` | `resample_signal` | substantial |
| `sftrans` | `signal_iir_design` | `splane_frequency_transform` | substantial |
| `sgolay` | `signal_sgolay` | `savitzky_golay` | substantial |
| `sgolayfilt` | `signal_sgolay` | `savitzky_golay_filter` | substantial |
| `specgram` | `signal_analysis` | `spectrogram` | substantial |
| `triang` | `signal_windows` | `triangular_window` | complete |
| `unwrap` | `signal_filters` | `unwrap_phase` | partial |

## Important compatibility notes

- Derived types (`arma_filter`, `zpg_filter`, `filter_order`, and related result
  types) replace R lists/S3 objects.
- `remez_filter` currently uses a deterministic weighted least-squares dense-grid
  design. It is useful computationally but is **not** an exact translation of
  the upstream Parks-McClellan C backend; the mapping is therefore partial.
- `fftfilt`, `ifft`, `freqz`, `fir2`, `chebwin`, and `specgram` use direct DFT or
  convolution kernels rather than an external FFT package. Numerical intent is
  preserved, but asymptotic performance differs.
- `interp1_signal(..., method="spline")` uses PCHIP as a deterministic fallback;
  linear, nearest, PCHIP, and local cubic modes are implemented directly.
- `roots` uses Durand-Kerner iteration and does not reproduce R's optional
  eigenvalue-based method exactly.
- Plotting, print methods, R formula/list coercion, names/dimnames, and warning
  text are intentionally outside scope.

## License and provenance

The upstream `signal` package declares `GPL-2`. This translation is distributed
under GPL-2.0-only. Upstream copyright, authorship, citation, Octave-derived
provenance, and Parks-McClellan notices are retained in `NOTICE.md` and the
`upstream/` directory. See those files before redistributing derived work.
