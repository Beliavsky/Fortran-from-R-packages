# Porting notes

## Representation

The R package stores coefficients in several packed S3 structures with first/last boundary databases. The Fortran translation uses typed allocatable arrays. Ordinary and stationary 1-D details use R-style level numbering (0 = coarsest detail), while the internal packet tree is root-first; explicit R-style packet adapters are provided.

## Boundaries and filters

The principal translated transform path is periodic and real-valued. Symmetric and zero-boundary first/last coefficient bookkeeping is translated, including `first.last.dh`, but the corresponding transform kernels are not yet reproduced. Real Daubechies/Coiflet/Yates/Littlewood-Paley filters are represented locally. The default Lina-Mayrand 3.1 periodic decimated complex transform and inverse are translated for complex thresholding; alternate complex filters and translation-invariant complex transforms are not yet exposed. Geronimo and Donovan3 multiple-wavelet filter tables, `mfirst.last` bookkeeping, all upstream pre/postfilter branches, periodic `mwd`/`mwr` transforms, and the principal `threshold.mwd` hard/soft paths are translated. The Fortran representation uses typed per-level matrices rather than R's packed first/last coefficient vectors.

## Shared transform backend

Low-level ordinary decimated 1-D steps and multidimensional transform kernels reuse the sibling `waveslim` package. The 1-D `wst` path follows upstream wavethresh's recursive unshifted/one-sample-shifted `wavepackst` construction and its exact `convolveC`/`convolveD` indexing and signs. Stationary reconstruction uses the matching `conbar` convention through `wst_conbar`; this is deliberately separate from the phase convention used by the shared ordinary DWT backend. `wr_wst` then follows recursive average-basis synthesis. The autocorrelation-wavelet spectrum correction reuses sibling `rfortran-linalg::inverse_matrix` rather than embedding a matrix solver.

## Thresholding

Manual hard/soft, universal, LSuniversal, and SURE-style paths are represented. The Ogden-Parzen rule-1 and rule-2 sequential threshold procedures (`TOgetthrda*`, `TOonebyone*`, and typed `TOthreshda*`) are also translated, including the upstream rule-2 empirical calibration table. `threshold.mwd` covers componentwise and joint multiple-wavelet hard/soft shrinkage with manual/universal thresholds and robust or ordinary covariance. `WaveletCV`, `FullWaveletCV`, `GetRSSWST`, scalar and level-specific `wvcvlrss`, scalar-threshold `wstCV`, and vector-threshold `wstCVl` are translated without plotting or verbose R output. The `wstCVl` search uses bounded coordinate-wise golden-section minimization instead of R's `nlminb`, so optimizer trajectories are not expected to be identical. `BAYES.THR` implements the periodic real-valued empirical-Bayes shrinkage path, including the upstream likelihood grid and posterior-median calculation. Complex `cthresh` supports the default Lina-Mayrand 3.1 periodic transform, robust covariance construction, mws hard/soft thresholding, and empirical-Bayes hard/soft/posterior-mean rules. Its `find_parameters` path replaces R's `optim`/`nlminb` and historical external NAG executable with a deterministic bounded coordinate search, so fitted parameters need not follow identical optimizer trajectories.

## Additional numerical utilities

`compgrot`, `guyrot`, `rotateback`, `cns`, `first.last`, `first.last.dh`, `mfirst.last`, `makegrid`, `getarrvec`, `numtonv`, and scaling-function refinement are translated. Stationary packet access includes `getpacket.wpst`, `wpst2m`, `wpst2discr`, and the numerical extraction core of `wpstREGR`; typed `makewpstRO` regression and `makewpstDO` discrimination objects reuse those packet features. `BMdiscr` supplies a pooled-covariance LDA core and `wpstCLASS` applies it to selected packet features; this intentionally does not reproduce every rank/SVD edge case of `MASS::lda`. `Best1DCols` and `bestm` expose their correlation-selection cores without recreating R packet/S3 result classes. Decimated `AutoBasis`/`InvBasis.wp`, minimum-entropy `MaNoVe.wp` and `MaNoVe.wst`, stationary `InvBasis.wst`, recursive `av.basis`, seeded `LSWsim`, deterministic `checkmyews`, autocorrelation-wavelet construction (`PsiJ`/`PsiJmat`/`ipndacw`), raw local periodograms, and the unsmoothed stationary `ewspec` correction are also available. Specialized 2-D stationary-packet access/mutation remains untranslated. `wvmoments` computes moments by direct interpolation and trapezoidal quadrature over the translated refinement samples rather than routing through R's `draw.default`, `approx`, and `integrate`.

## Density projection

`Chires5`, `Chires6`, and `denproj` use a direct free-form Fortran translation of the upstream Daubechies-Lagarias binary transition-matrix product. The typed `density_projection_t` stores the coefficient translation range, filter and resolution metadata, and the same upper covariance-band layout produced by `Chires6`. `evaluate_density` supplies the non-plotting numerical evaluation performed by R's `denplot`. `cwavde` translates the separate sampled-basis estimator, including support bounds, coefficient hard thresholding, the upstream interpolation convention, and density-grid evaluation. Its scaling and wavelet samples are explicit typed inputs instead of being generated through R's plotting-oriented `draw.default`. `denwd` translates the arbitrary-index zero-boundary decomposition into `density_wavelet_t`, retaining every scaling/detail level and its integer bounds. `dencvwd` expands the compact input bands to a symmetric covariance matrix and propagates it with the zero-boundary scaling/detail operators. This dense formulation matches the upstream variance diagonals and is easier to audit, at the cost of more memory and arithmetic than the specialized C band algorithm.

`wd_int` and `wr_int` directly translate the upstream adaptive interval-wavelet C implementation. The generated `wavethresh_interval_coefficients` module preserves the GPL-covered interior, left/right boundary, and preconditioning tables for all eight supported filter orders; `scripts/generate_interval_coefficients.py` regenerates it deterministically from the retained upstream C source. The typed `interval_wavelet_t` records the filter order used at each scale so reconstruction can reverse both adaptive order changes and optional boundary preconditioning.

The stationary 2-D representation now retains the smooth matrix at every level. `getpacket_wst2d` and `putpacket_wst2d` translate the upstream base-four path-to-coordinate rule and select the requested S, H, V, or D block from the corresponding typed coefficient matrix. This avoids reproducing the upstream monolithic three-dimensional packed array while preserving packet addressing and mutation semantics.

`irregwd` uses the interpolation weights retained by `makegrid` as an explicit linear map. It transforms each observation's influence vector and sums squared detail responses to obtain the coefficient variance factors returned by `accessc`. This clear reference implementation replaces the upstream specialized `computec` storage kernel and supports variance-adjusted manual and universal `threshold_irregwd` paths. Its setup cost grows with the number of original observations; a future matrix or sparse batched implementation could preserve the same semantics more efficiently.

## Dense 2-D storage

The typed `imwd_t` representation is already dense and contains no separate compressed boundary payload. Consequently `compress_imwd`/`uncompress_imwd` are identity copies, and `imwr.imwdc`/`threshold.imwdc` map partially to the same reconstruction and threshold kernels used for dense objects. The R compressed-list representation and memory-size behavior are not reproduced.

## Fourier calculations

`rfft`, `rfftinv`, and `rfftwt` use direct dependency-free Fourier sums instead of calling an FFT library. Their coefficient normalization follows upstream `rfft`/`rfftinv`, including the unpaired Nyquist cosine term for even-length transforms. Complexity is O(n^2).

## RNG

Simulation routines use the Fortran intrinsic RNG with optional deterministic seed expansion. They do not reproduce R's RNG stream bit-for-bit.

## 3-D filters

The 3-D adapter uses the filter names currently supported by sibling `waveslim`; not every wavethresh filter family/number has a 3-D backend mapping.
