# wavethresh for modern Fortran

Modern free-form Fortran/FPM translation of a substantial computational subset of the R package **wavethresh 4.7.3** (Guy Nason and contributors). The port emphasizes portable numerical wavelet transforms and shrinkage without R object, plotting, or interactive infrastructure.

## Build

Place this directory beside the shared dependencies in `Fortran-from-R-packages` and run:

```text
fpm build
fpm test
fpm run wavethresh_demo
fpm run --example basic_dwt
fpm run --example denoising
```

The manifest uses sibling path dependencies on `../rfortran-core`, `../rfortran-linalg`, and `../waveslim`. `rfortran-linalg` supplies the shared matrix inverse used by spectrum correction and, in the repository layout, obtains LAPACK through its pinned FPM dependency rather than a system `-llapack`/`-lblas` link. No system BLAS, LAPACK, ARPACK, or FFTW link is declared by this package.

## Implemented numerical scope

- Daubechies extremal-phase and least-asymmetric filters, Coiflets, Yates, Littlewood-Paley, selected complex descriptors, and complete Geronimo/Donovan3 multiple-wavelet filter tables.
- Periodic real 1-D decimated and stationary transforms and inverse transforms.
- Periodic Geronimo and Donovan3 multiple-wavelet transforms, inverse transforms, and their upstream pre/postfilter families.
- Decimated and stationary wavelet-packet trees with typed coefficient access/update, packet rotations, and stationary feature-matrix/regression extraction.
- Periodic real 2-D decimated/stationary transforms and periodic real 3-D transforms.
- Manual, universal, levelwise, LSuniversal, SURE-style, and periodic complex empirical-Bayes thresholding for the translated object families.
- Wavelet cross-validation (`rsswav`, `Crsswav`, `WaveletCV`, `FullWaveletCV`, stationary `wstCV`) and legacy reconstruction (`denwr`).
- Real Fourier-series helpers, SURE, robust MAD/covariance utilities, entropy/norms, support calculations, transform matrices, scaling-function refinement, adaptive interval wavelets, high-resolution density projection, grid interpolation, and first/last bookkeeping.
- Donoho-Johnstone/test-data signals, chirps, Haar simulations, LSW simulation/checking, autocorrelation-wavelet construction and unsmoothed evolutionary-spectrum correction, and claw-distribution density/CDF/simulation.

## API conventions

R lists/S3 objects are replaced by typed derived types such as `wd_t`, `wp_t`, `mwd_t`, `imwd_t`, and `wd3d_t`. Detail levels follow wavethresh conventions where level 0 is the coarsest detail for ordinary/stationary transforms. Packet wrappers `getpacket_r`, `putpacket_r`, and `access_d_wp` preserve wavethresh packet-level addressing even though the internal tree is stored root-first.

Periodic real transforms remain the principal translated path. Symmetric/zero boundary bookkeeping is translated, while symmetric/zero transform kernels remain gaps. Adaptive interval transforms, irregular-grid decomposition/thresholding, periodic Lina-Mayrand 3.1 complex transforms and thresholding, Geronimo/Donovan3 multiple-wavelet transforms, stationary packet discrimination, and the core unsmoothed LSW/evolutionary-spectrum path are translated with the limitations documented in `API_COVERAGE.md` and `PORTING_NOTES.md`.

## Translation coverage

Package status: **substantial**.

**161 of 162 (99.4%)** exported computational R functions have a meaningful complete, substantial, or partial mapping. This fraction measures mapped computational functions, not complete R compatibility or exact parity of every option. Plotting, printing, formatting, interactive, migration-only, and dispatch-only exports are excluded from the denominator.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
|---|---|---|---|
| `test.dataCT` | `wavethresh_signals` | `test_data_ct` | substantial |
| `filter.select` | `wavethresh_filters` | `filter_select` | substantial |
| `accessC.wd` | `wavethresh_access` | `access_c_wd` | substantial |
| `accessC.wst` | `wavethresh_access` | `access_c_wd` | substantial |
| `accessD.wd` | `wavethresh_access` | `access_d_wd` | substantial |
| `accessD.wd3D` | `wavethresh_access` | `access_d_wd3d` | substantial |
| `accessD.wp` | `wavethresh_access` | `access_d_wp` | substantial |
| `accessD.wst` | `wavethresh_access` | `access_d_wd` | substantial |
| `AvBasis.wst` | `wavethresh_transform_1d` | `wr_wst` | substantial |
| `AvBasis.wst2D` | `wavethresh_transform_nd` | `iwst2d` | substantial |
| `c2to4` | `wavethresh_stats` | `c2to4` | complete |
| `dclaw` | `wavethresh_density` | `dclaw` | complete |
| `denwr` | `wavethresh_cv` | `denwr` | substantial |
| `pclaw` | `wavethresh_density` | `pclaw` | complete |
| `rclaw` | `wavethresh_density` | `rclaw` | substantial |
| `compare.filters` | `wavethresh_filters` | `compare_filters` | substantial |
| `compress.imwd` | `wavethresh_access` | `compress_imwd` | partial |
| `convert.wd` | `wavethresh_access` | `convert_wd` | partial |
| `convert.wst` | `wavethresh_access` | `convert_wd` | partial |
| `Crsswav` | `wavethresh_cv` | `crsswav` | substantial |
| `Cthreshold` | `wavethresh_threshold` | `threshold_wd` | substantial |
| `rsswav` | `wavethresh_cv` | `rsswav` | substantial |
| `WaveletCV` | `wavethresh_cv` | `wavelet_cv` | substantial |
| `DJ.EX` | `wavethresh_signals` | `dj_ex` | substantial |
| `doppler` | `wavethresh_signals` | `doppler` | complete |
| `simchirp` | `wavethresh_signals` | `simchirp` | complete |
| `dof` | `wavethresh_stats` | `dof` | substantial |
| `support` | `wavethresh_basis` | `wavelet_support` | complete |
| `firstdot` | `wavethresh_stats` | `firstdot` | complete |
| `GenW` | `wavethresh_basis` | `gen_w` | substantial |
| `getpacket.wp` | `wavethresh_access` | `getpacket_r` | substantial |
| `HaarConcat` | `wavethresh_signals` | `haar_concat` | substantial |
| `HaarMA` | `wavethresh_signals` | `haar_ma` | substantial |
| `imwd` | `wavethresh_transform_nd` | `imwd` | substantial |
| `imwr.imwd` | `wavethresh_transform_nd` | `imwr` | substantial |
| `IsPowerOfTwo` | `wavethresh_transform_1d` | `nlevels_from_length` | complete |
| `l2norm` | `wavethresh_stats` | `l2norm` | complete |
| `linfnorm` | `wavethresh_stats` | `linfnorm` | complete |
| `levarr` | `wavethresh_stats` | `levarr` | complete |
| `logabs` | `wavethresh_stats` | `logabs` | complete |
| `madmad` | `wavethresh_stats` | `madmad` | complete |
| `ssq` | `wavethresh_stats` | `ssq` | complete |
| `rcov` | `wavethresh_stats` | `robust_covariance` | substantial |
| `newsure` | `wavethresh_stats` | `newsure` | complete |
| `sure` | `wavethresh_stats` | `sure` | complete |
| `nullevels.imwd` | `wavethresh_access` | `nullevels_imwd` | substantial |
| `nullevels.wd` | `wavethresh_access` | `nullevels_wd` | substantial |
| `nullevels.wst` | `wavethresh_access` | `nullevels_wd` | substantial |
| `putC.wd` | `wavethresh_access` | `put_c_wd` | substantial |
| `putC.wst` | `wavethresh_access` | `put_c_wd` | substantial |
| `putD.wd` | `wavethresh_access` | `put_d_wd` | substantial |
| `putD.wd3D` | `wavethresh_access` | `put_d_wd3d` | substantial |
| `putD.wst` | `wavethresh_access` | `put_d_wd` | substantial |
| `putpacket.wp` | `wavethresh_access` | `putpacket_r` | substantial |
| `rfft` | `wavethresh_fourier` | `rfft` | complete |
| `rfftinv` | `wavethresh_fourier` | `rfftinv` | complete |
| `rfftwt` | `wavethresh_fourier` | `rfftwt` | complete |
| `Shannon.entropy` | `wavethresh_stats` | `shannon_entropy` | complete |
| `threshold.imwd` | `wavethresh_threshold` | `threshold_imwd` | substantial |
| `threshold.wd` | `wavethresh_threshold` | `threshold_wd` | substantial |
| `threshold.wd3D` | `wavethresh_threshold` | `threshold_wd3d` | substantial |
| `threshold.wp` | `wavethresh_threshold` | `threshold_wp` | partial |
| `threshold.wst` | `wavethresh_threshold` | `threshold_wd` | substantial |
| `TOkolsmi.chi2` | `wavethresh_stats` | `kolsmi_chi2` | substantial |
| `TOshrinkit` | `wavethresh_stats` | `shrink_soft` | complete |
| `tpwd` | `wavethresh_transform_nd` | `tpwd` | substantial |
| `tpwr` | `wavethresh_transform_nd` | `tpwr` | substantial |
| `wd` | `wavethresh_transform_1d` | `wd` | substantial |
| `wd3D` | `wavethresh_transform_nd` | `wd3d` | substantial |
| `wp` | `wavethresh_transform_1d` | `wp` | substantial |
| `wpst` | `wavethresh_transform_1d` | `wpst` | partial |
| `wr.wd` | `wavethresh_transform_1d` | `wr_wd` | substantial |
| `wr3D` | `wavethresh_transform_nd` | `wr3d` | substantial |
| `wst` | `wavethresh_transform_1d` | `wst` | substantial |
| `wst2D` | `wavethresh_transform_nd` | `wst2d` | substantial |
| `Best1DCols` | `wavethresh_utilities` | `best_1d_cols` | partial |
| `bestm` | `wavethresh_utilities` | `bestm` | partial |
| `cns` | `wavethresh_utilities` | `cns` | substantial |
| `compgrot` | `wavethresh_utilities` | `compgrot` | complete |
| `guyrot` | `wavethresh_utilities` | `guyrot` | complete |
| `rotateback` | `wavethresh_utilities` | `rotateback` | complete |
| `ScalingFunction` | `wavethresh_utilities` | `scaling_function` | substantial |
| `first.last` | `wavethresh_utilities` | `first_last` | substantial |
| `makegrid` | `wavethresh_utilities` | `makegrid` | substantial |
| `TOgetthrda1` | `wavethresh_threshold_extra` | `to_getthrda1` | complete |
| `TOgetthrda2` | `wavethresh_threshold_extra` | `to_getthrda2` | complete |
| `TOonebyone1` | `wavethresh_threshold_extra` | `to_one_by_one1` | complete |
| `TOonebyone2` | `wavethresh_threshold_extra` | `to_one_by_one2` | complete |
| `TOthreshda1` | `wavethresh_threshold_extra` | `to_threshda1, to_threshda1_thresholds` | substantial |
| `TOthreshda2` | `wavethresh_threshold_extra` | `to_threshda2, to_threshda2_thresholds` | substantial |
| `wvmoments` | `wavethresh_utilities` | `wvmoments` | partial |
| `accessD.wpst` | `wavethresh_access` | `access_d_wpst` | substantial |
| `getarrvec` | `wavethresh_utilities` | `getarrvec` | complete |
| `make.dwwt` | `wavethresh_utilities` | `make_dwwt` | partial |
| `numtonv` | `wavethresh_utilities` | `numtonv` | substantial |
| `putD.wp` | `wavethresh_access` | `put_d_wp_level` | substantial |
| `wpst2m` | `wavethresh_utilities` | `wpst2m` | substantial |
| `wpstREGR` | `wavethresh_utilities` | `wpst_regr` | partial |
| `conbar` | `wavethresh_transform_1d` | `conbar` | substantial |
| `first.last.dh` | `wavethresh_utilities` | `first_last_dh` | substantial |
| `getpacket.wst` | `wavethresh_access` | `getpacket_wst` | substantial |
| `putpacket.wst` | `wavethresh_access` | `putpacket_wst` | substantial |
| `rm.det` | `wavethresh_access` | `rm_det` | partial |
| `imwr.imwdc` | `wavethresh_transform_nd` | `imwr` | partial |
| `threshold.imwdc` | `wavethresh_threshold` | `threshold_imwd` | partial |
| `uncompress.imwdc` | `wavethresh_access` | `uncompress_imwd` | partial |
| `wd.dh` | `wavethresh_transform_1d` | `wd` | partial |
| `mfilter.select` | `wavethresh_filters` | `mfilter_select` | complete |
| `mfirst.last` | `wavethresh_utilities` | `mfirst_last` | substantial |
| `accessC.mwd` | `wavethresh_access` | `access_c_mwd` | substantial |
| `accessD.mwd` | `wavethresh_access` | `access_d_mwd` | substantial |
| `mpostfilter` | `wavethresh_multiwavelet` | `mpostfilter` | substantial |
| `mprefilter` | `wavethresh_multiwavelet` | `mprefilter` | substantial |
| `mwd` | `wavethresh_multiwavelet` | `mwd` | substantial |
| `mwr` | `wavethresh_multiwavelet` | `mwr` | substantial |
| `putC.mwd` | `wavethresh_access` | `put_c_mwd` | substantial |
| `putD.mwd` | `wavethresh_access` | `put_d_mwd` | substantial |
| `wr.mwd` | `wavethresh_multiwavelet` | `mwr` | substantial |
| `AutoBasis` | `wavethresh_basis` | `auto_basis` | substantial |
| `av.basis` | `wavethresh_transform_1d` | `av_basis` | substantial |
| `LSWsim` | `wavethresh_signals` | `lsw_sim` | substantial |
| `FullWaveletCV` | `wavethresh_cv` | `full_wavelet_cv` | substantial |
| `GetRSSWST` | `wavethresh_cv` | `get_rss_wst` | substantial |
| `wstCV` | `wavethresh_cv` | `wst_cv` | substantial |
| `wvcvlrss` | `wavethresh_cv` | `wvcvlrss, wvcvlrss_levels` | substantial |
| `InvBasis.wp` | `wavethresh_basis` | `inv_basis_wp` | substantial |
| `threshold.mwd` | `wavethresh_threshold` | `threshold_mwd` | substantial |
| `wpst2discr` | `wavethresh_utilities` | `wpst2discr` | substantial |
| `PsiJ` | `wavethresh_spectrum` | `psi_j` | substantial |
| `PsiJmat` | `wavethresh_spectrum` | `psi_j_mat` | substantial |
| `ipndacw` | `wavethresh_spectrum` | `ipndacw` | substantial |
| `LocalSpec.wd` | `wavethresh_spectrum` | `local_spec_wd` | partial |
| `LocalSpec.wst` | `wavethresh_spectrum` | `local_spec_wst` | partial |
| `ewspec` | `wavethresh_spectrum` | `ewspec` | partial |
| `CWCV` | `wavethresh_cv` | `cwcv` | partial |
| `checkmyews` | `wavethresh_signals` | `check_my_ews` | substantial |
| `getpacket.wpst` | `wavethresh_access` | `getpacket_wpst` | substantial |
| `InvBasis.wst` | `wavethresh_basis` | `inv_basis_wst` | substantial |
| `makewpstRO` | `wavethresh_utilities` | `makewpst_ro` | substantial |
| `MaNoVe.wp` | `wavethresh_basis` | `manove_wp` | substantial |
| `MaNoVe.wst` | `wavethresh_basis` | `manove_wst` | substantial |
| `BMdiscr` | `wavethresh_utilities` | `bm_discr` | substantial |
| `makewpstDO` | `wavethresh_utilities` | `makewpst_do` | substantial |
| `wpstCLASS` | `wavethresh_utilities` | `wpst_class` | substantial |
| `wstCVl` | `wavethresh_cv` | `wst_cvl` | partial |
| `BAYES.THR` | `wavethresh_bayes` | `bayes_thr` | substantial |
| `Chires5` | `wavethresh_density` | `chires5` | complete |
| `Chires6` | `wavethresh_density` | `chires6` | complete |
| `denproj` | `wavethresh_density` | `denproj` | complete |
| `CWavDE` | `wavethresh_density` | `cwavde` | substantial |
| `denwd` | `wavethresh_density` | `denwd` | complete |
| `dencvwd` | `wavethresh_density` | `dencvwd` | substantial |
| `wd.int` | `wavethresh_interval` | `wd_int` | complete |
| `wr.int` | `wavethresh_interval` | `wr_int` | complete |
| `getpacket.wst2D` | `wavethresh_access` | `getpacket_wst2d` | substantial |
| `putpacket.wst2D` | `wavethresh_access` | `putpacket_wst2d` | substantial |
| `irregwd` | `wavethresh_irregular` | `irregwd` | substantial |
| `accessc` | `wavethresh_irregular` | `accessc` | substantial |
| `threshold.irregwd` | `wavethresh_irregular` | `threshold_irregwd` | substantial |
| `cthresh` | `wavethresh_complex_threshold` | `cthresh` | substantial |
| `find.parameters` | `wavethresh_complex_threshold` | `find_parameters` | substantial |

## License and provenance

The upstream package declares `GPL (>= 2)`; this translation is distributed as **GPL-2.0-or-later**. The GPL version 2 text is in `LICENSE`. The supplied upstream source tree is retained under `upstream/` except for its generated `build/partial.rdb` artifact, which is excluded by the release hygiene rule and recorded by path/hash in `provenance/OMITTED_UPSTREAM_BUILD_ARTIFACTS.txt`. Source and archive checksums are under `provenance/`.

Upstream attribution includes Guy Nason and the contributors named in `DESCRIPTION`; source-specific notices for Tim Downie, Arne Kovac, Piotr Fryzlewicz, Markus Monnerjahn, and others are preserved verbatim in the retained upstream source. See `NOTICE.md`.
