# API coverage

Coverage basis: distinct exported computational R functions from upstream
`NAMESPACE`. Presentation-only exports (`freqs_plot`, `freqz_plot`, `zplane`),
registered print/plot methods, datasets, and internal helpers are excluded.

Package status: **substantial**. Mapped functions: **58 of 58 (100%)**.
The 100% fraction means every function in this coverage basis has a meaningful
Fortran mapping. It does **not** mean exact R API, S3, backend, or bit-for-bit
compatibility.

| R function | Status | Material compatibility note |
|---|---|---|
| `bartlett` | complete | Direct numerical translation of the Bartlett window for vector output. |
| `bilinear` | substantial | Bilinear zero/pole/gain transform is mapped; R S3 dispatch and object construction are represented by typed Fortran values. |
| `blackman` | complete | Direct numerical translation of the Blackman window. |
| `boxcar` | complete | Direct numerical translation of the rectangular window. |
| `butter` | substantial | Low/high/bandpass/bandstop Butterworth design is mapped for analog and digital filters; R S3 dispatch is omitted. |
| `buttord` | substantial | Filter-order computation is mapped; R warnings and list/class construction are omitted. |
| `chebwin` | substantial | Dolph-Chebyshev window computation is mapped using the package-local DFT; numerical roundoff and runtime differ from R fft. |
| `cheb1ord` | substantial | Chebyshev-I order computation is mapped; band-reject handling remains limited as in the upstream implementation. |
| `cheby1` | substantial | Chebyshev-I prototype and frequency transforms are mapped; R S3 dispatch is omitted. |
| `cheby2` | substantial | Chebyshev-II prototype and frequency transforms are mapped; R S3 dispatch is omitted. |
| `chirp` | substantial | Linear, quadratic, and logarithmic chirp calculations are mapped; R argument matching and warning text are omitted. |
| `decimate` | substantial | FIR/IIR antialias filtering and downsampling are mapped; endpoint behavior follows the translated filters rather than R object methods. |
| `ellip` | substantial | Elliptic prototype/design and transforms are mapped; R S3 dispatch and warning behavior are omitted. |
| `ellipord` | substantial | Elliptic order computation is mapped with package-local elliptic-integral evaluation; R warnings/list construction are omitted. |
| `FftFilter` | partial | The computational filter object is represented by a Fortran derived type; requested FFT length is retained as metadata, without an R S3 object. |
| `filter` | substantial | Direct-form FIR/IIR filtering and initial-state inputs are mapped; R S3 dispatch, recycling, and exact stats::filter state conventions are not reproduced. |
| `fftfilt` | substantial | FIR convolution result is mapped, but this standalone translation uses direct convolution instead of overlap-add FFT acceleration. |
| `MedianFilter` | partial | The running-median filter configuration is represented by a Fortran derived type rather than an R S3 object. |
| `medfilt1` | partial | Running median computation is mapped; boundary handling is local-window median and is not an exact clone of stats::runmed endrule behavior. |
| `spencerFilter` | substantial | The 15-term Spencer FIR coefficients are mapped into the Fortran ARMA filter type; R class construction is omitted. |
| `spencer` | substantial | Spencer smoothing and NaN endpoints are mapped; R names/attributes are omitted. |
| `FilterOfOrder` | partial | Filter-order metadata is represented by a Fortran derived type; arbitrary R list fields in ... are not supported. |
| `an` | complete | Degree-to-unit-phasor computation is mapped as a pure elemental function. |
| `roots` | partial | Polynomial roots are computed with Durand-Kerner iteration; the alternate R eigenvalue method and exact polyroot behavior are not reproduced. |
| `Arma` | partial | ARMA coefficient storage is represented by a typed Fortran value rather than an R S3 list. |
| `as.Arma` | partial | Zero/pole/gain and moving-average coefficient conversions are supported numerically; R S3 identity/dispatch behavior is omitted. |
| `Ma` | partial | Moving-average coefficients are represented by the shared Fortran ARMA filter type with denominator one. |
| `Zpg` | partial | Zero-pole-gain storage is represented by a typed Fortran value rather than an R S3 list. |
| `as.Zpg` | partial | Polynomial-to-zero/pole/gain conversion is mapped; R S3 identity/dispatch behavior is omitted. |
| `polyval` | substantial | Polynomial evaluation at complex points is mapped; R vector recycling/coercion rules are omitted. |
| `conv` | complete | Linear convolution of real vectors is mapped directly. |
| `ifft` | substantial | Inverse complex transform is mapped with a package-local O(N^2) DFT instead of R fft, so performance and roundoff differ. |
| `filtfilt` | substantial | Forward/reverse zero-phase filtering with upstream zero padding is mapped; R S3 dispatch is omitted. |
| `fir1` | substantial | Windowed FIR low/high/bandpass/bandstop design and scaling are mapped; R argument matching and attributes are omitted. |
| `fir2` | partial | Frequency-sampling FIR design is mapped with a local DFT; duplicate-knot ramp handling is numerically approximated rather than cloning all R interpolation details. |
| `flattopwin` | substantial | Symmetric and periodic flat-top windows are mapped, including the upstream periodic n=1 edge case. |
| `freqs` | substantial | Analog complex frequency response is mapped; R S3 class, printing, and plotting behavior are omitted. |
| `freqz` | substantial | Digital complex frequency response on generated or supplied grids is mapped; R S3 class, printing, plotting, and FFT acceleration are omitted. |
| `gausswin` | complete | Direct numerical translation of the Gaussian window. |
| `grpdelay` | substantial | Group-delay computation and frequency grid are mapped; R S3 class, printing, and plotting are omitted. |
| `hamming` | complete | Direct numerical translation of the Hamming window. |
| `hanning` | complete | Direct numerical translation of the Hann window exposed upstream as hanning. |
| `impz` | substantial | Impulse response and time grid are mapped; R S3 printing/plotting and exact default-length heuristics are simplified. |
| `interp1` | partial | Linear, nearest, PCHIP, and local cubic interpolation are mapped; spline currently falls back to PCHIP, and R matrix/recycling behavior is omitted. |
| `interp` | substantial | Integer-rate zero insertion, low-pass FIR filtering, gain, and delay trimming are mapped; FFT acceleration is not used. |
| `kaiser` | substantial | Kaiser window computation is mapped with a package-local modified-Bessel approximation. |
| `kaiserord` | substantial | Kaiser order/beta/cutoff estimation is mapped; R list/class representation and warnings are omitted. |
| `levinson` | substantial | Levinson-Durbin recursion, prediction variance, and reflection coefficients are mapped; R return-list structure is replaced by output arguments. |
| `pchip` | substantial | Shape-preserving piecewise cubic Hermite interpolation is mapped numerically; upstream piecewise-polynomial R object construction is omitted. |
| `poly` | partial | Polynomial construction from roots is mapped; the R matrix-characteristic-polynomial input mode is not translated. |
| `remez` | partial | API-level FIR design is provided using weighted least squares on a dense grid; it is not an exact Parks-McClellan/Remez equiripple translation. |
| `resample` | substantial | Windowed-sinc resampling, including non-integer p/q ratios and antialias filtering, is mapped; FFT acceleration and some R endpoint details differ. |
| `sftrans` | substantial | S-plane low/high/bandpass/bandstop zero/pole/gain transforms are mapped; R S3 dispatch is omitted. |
| `sgolay` | substantial | Savitzky-Golay coefficient construction is mapped with a package-local pivoted linear solver rather than MASS::ginv; R class attributes are omitted. |
| `sgolayfilt` | substantial | Savitzky-Golay filtering and derivative scaling are mapped; R S3 filter dispatch is omitted. |
| `specgram` | substantial | Windowing, overlap, per-frame DFT, and numeric spectrogram outputs are mapped; R plotting/printing and FFT acceleration are omitted. |
| `triang` | complete | Direct numerical translation of the triangular window. |
| `unwrap` | partial | One-dimensional phase unwrapping is mapped; R multi-dimensional dim selection and array attributes are omitted. |

## Additional translated internal computational helper

`signal:::fractdiff` is translated as `signal_filters::fractional_difference`.
It is intentionally not included in the exported-R-function coverage denominator.
