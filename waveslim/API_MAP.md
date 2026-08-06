# API map

R names use dots; corresponding Fortran names use underscores. Internal R
helpers are generally incorporated into the public Fortran routine that uses
them rather than exported separately.

## Filters and basic transforms

| R API | Fortran API | Status |
|---|---|---|
| `wave.filter`, selector helpers | `wave_filter` | Implemented for all standard package filters |
| `qmf` | `qmf` | Implemented |
| `hilbert.filter` | `hilbert_filter` | Implemented for K3L3, K3L5, K4L2, K4L4, K5L7, K6L6 |
| `dwt`, `idwt` | `dwt`, `idwt` | Implemented |
| `dwt.nondyadic` | `dwt_nondyadic` | Implemented |
| `modwt`, `imodwt` | `modwt`, `imodwt` | Implemented |
| `mra` | `mra` | Implemented |
| `brick.wall` | `brick_wall` | Implemented |
| `phase.shift` | `phase_shift` | Implemented |
| `up.sample` | `up_sample` | Implemented |

## Wavelet packets

| R API | Fortran API | Status |
|---|---|---|
| `dwpt`, `idwpt` | `dwpt`, `idwpt` | Implemented |
| `modwpt` | `modwpt` | Implemented |
| `basis`, `ortho.basis` | `packet_basis`, `ortho_basis` | Implemented |
| `phase.shift.packet` | `phase_shift_packet` | Implemented |
| `dwpt.brick.wall` | `dwpt_brick_wall` | Implemented |
| `css.test` | `css_test` | Implemented |
| `entropy.test` | `entropy_test` | Implemented |
| `cpgram.test` | `cpgram_test` | Implemented |
| `portmanteau.test` | `portmanteau_test` | Implemented |
| `dwpt.boot` | `dwpt_boot` | Implemented with seeded Fortran RNG |
| `dwpt.2d`, `idwpt.2d` | `dwpt_2d`, `idwpt_2d` | Implemented |

## Multidimensional transforms

| R API | Fortran API | Status |
|---|---|---|
| `dwt.2d`, `idwt.2d` | `dwt_2d`, `idwt_2d` | Implemented |
| `modwt.2d`, `imodwt.2d` | `modwt_2d`, `imodwt_2d` | Implemented |
| `mra.2d` | `mra_2d` | Implemented |
| `dwt.3d`, `idwt.3d` | `dwt_3d`, `idwt_3d` | Implemented |
| `modwt.3d`, `imodwt.3d` | `modwt_3d`, `imodwt_3d` | Implemented |
| `mra.3d` | `mra_3d` | Implemented |
| `shift.2d` | `shift_2d` | Implemented |
| `convolve2D` | `convolve_2d` | Implemented |
| `brick.wall.2d` | `brick_wall_2d` | Implemented |

## Dual-tree and Hilbert transforms

| R API | Fortran API | Status |
|---|---|---|
| `farras`, `FSfarras`, `dualfilt1` | `farras_filters`, `fsfarras_filters`, `dualfilt1_filters` | Implemented |
| `afb`, `sfb` | `afb`, `sfb` | Implemented |
| `dualtree`, `idualtree` | `dualtree`, `idualtree` | Implemented |
| `afb2D`, `sfb2D` | `afb2d`, `sfb2d` | Implemented |
| `dualtree2D`, `idualtree2D` | `dualtree_2d`, `idualtree_2d` | Implemented |
| `dwt.hilbert`, `idwt.hilbert` | `dwt_hilbert`, `idwt_hilbert` | Implemented |
| `modwt.hilbert`, `imodwt.hilbert` | `modwt_hilbert`, `imodwt_hilbert` | Implemented |
| `modhwt.coh`, `modhwt.phase` | `modhwt_coherence`, `modhwt_phase` | Implemented |
| seasonal coherence/phase | `modhwt_coherence_seasonal`, `modhwt_phase_seasonal` | Implemented |
| `cplxdual2D`, `icplxdual2D` | `dualtree_2d`, `idualtree_2d` | Stable two-tree equivalent; upstream routines call `pm` with the wrong arity |
| `modwpt.hilbert` | none | Not separately exposed; DWT/MODWT Hilbert and packet transforms are available |

## Denoising

| R API | Fortran API | Status |
|---|---|---|
| `soft` | `soft` | Implemented |
| `bishrink` | `bishrink` | Implemented |
| `manual.thresh` | `manual_thresh` | Implemented |
| `universal.thresh` | `universal_thresh` | Implemented |
| `universal.thresh.modwt` | `universal_thresh_modwt` | Implemented |
| `sure.thresh` | `sure_thresh` | Implemented |
| `hybrid.thresh` | `hybrid_thresh` | Implemented |
| `da.thresh` | `da_thresh` | Implemented |
| `denoise.dwt.2d` | `denoise_dwt_2d` | Implemented |
| `denoise.modwt.2d` | `denoise_modwt_2d` | Implemented |

## Statistics and spectral routines

| R API | Fortran API | Status |
|---|---|---|
| `wave.variance` | `wave_variance` | Implemented |
| `wave.covariance` | `wave_covariance` | Implemented |
| `wave.correlation` | `wave_correlation` | Implemented |
| `spin.covariance` | `spin_covariance` | Implemented |
| `spin.correlation` | `spin_correlation` | Implemented |
| `my.acf`, `my.ccf` | `my_acf`, `my_ccf` | Implemented |
| `per` | `periodogram` | Implemented |
| `sine.taper` | `sine_taper` | Implemented |
| `dpss.taper` | `dpss_taper` | Implemented with self-contained eigensolver |
| `wave.variance.2d` | `wave_variance_2d` | Implemented |
| `rotcumvar` | `rotcumvar` | Implemented |
| `wavelet.filter`, `cascade` | `wavelet_filter` | Implemented |
| `squared.gain` | `squared_gain` | Implemented |
| `testing.hov`, `mult.loc` | `testing_hov`, `mult_loc` | Implemented |

## Long-memory and simulation routines

| R API | Fortran API | Status |
|---|---|---|
| `fdp.sdf`, `bandpass.fdp` | `fdp_sdf`, `bandpass_fdp` | Implemented |
| `spp.sdf`, `bandpass.spp` | `spp_sdf`, `bandpass_spp` | Implemented |
| `spp2.sdf`, `bandpass.spp2` | `spp2_sdf`, `bandpass_spp2` | Implemented |
| `sfd.sdf` | `sfd_sdf` | Implemented |
| `spp.var` | `spp_variance` | Implemented |
| `Hypergeometric` | `hypergeometric_2f1` | Implemented |
| `hosking.sim` | `hosking_sim` | Implemented |
| `fdp.mle` | `fdp_mle` | Implemented numerical equivalent |
| `spp.mle` | `spp_mle` | Implemented numerical equivalent |
| `spp2.mle` | `spp2_mle` | Implemented numerical equivalent |
| `find.adaptive.basis` | `find_adaptive_basis` | Implemented |
| `bandpass.var.spp` | `bandpass_var_spp` | Implemented |
| `dwpt.sim` | `dwpt_sim` | Implemented |

## Omitted R infrastructure

`stackPlot`, `addmain`, `plot.dwt.2d`, datasets, S3 print/plot methods, and
R namespace/runtime helpers are omitted. Internal confidence-interval and
thresholding helpers are incorporated into the corresponding public numerical
routines rather than exposed as standalone procedures.
