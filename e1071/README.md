# e1071-fortran

Modern Fortran/FPM translation of the computational functionality of the R package `e1071` 1.7-17.

The library provides native Fortran implementations of support vector machines, fuzzy c-means/c-shell clustering, shortest paths, latent-class analysis, naive Bayes, generalized k-nearest neighbours, independent component analysis, bagged clustering, class/control matching, sparse CSR I/O, short-time Fourier transforms, discrete distributions, probability-plot calculations, tuning helpers, and the package's numerical utilities.

## Basic SVM example

```fortran
use e1071, only: dp, svm_model, svm_options, svm_c_classification, svm_radial, &
                 svm_fit_classification, svm_predict_classification

real(dp) :: x(6, 2)
integer :: y(6)
integer, allocatable :: prediction(:)
type(svm_options) :: options
type(svm_model) :: model

x = reshape([0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, &
             3.0_dp, 3.0_dp, 3.0_dp, 4.0_dp, 4.0_dp, 3.0_dp], &
            [6, 2], order=[2, 1])
y = [-1, -1, -1, 1, 1, 1]

options = svm_options(svm_type=svm_c_classification, kernel=svm_radial, &
                      gamma=0.5_dp, scale=.false.)
call svm_fit_classification(x, y, model, options)
call svm_predict_classification(model, x, prediction)
```

## Major translated areas

### Support vector machines

The native dense solver covers:

- C-SVC and nu-SVC;
- one-class SVM;
- epsilon-SVR and nu-SVR;
- linear, polynomial, radial-basis and sigmoid kernels;
- multiclass one-versus-one voting and decision values;
- predictor/response scaling;
- class weights;
- binary probability calibration and multiclass probability coupling;
- SVR probability-scale estimation;
- cross-validation;
- linear coefficients for compatible models.

The solver follows the LIBSVM 3.23 algorithms retained in the upstream `e1071` source. It uses dense `real(dp)` kernel matrices rather than LIBSVM's single-precision `Qfloat` cache/shrinking machinery, so optimization paths and support-vector coefficients are not claimed bit-for-bit identical even when predictions and objectives agree within solver tolerance.

CSR SVM fitting now constructs the kernel directly from sparse rows and follows upstream `e1071` by disabling predictor and response scaling for sparse inputs. The full dense training matrix is not materialized. Selected support vectors are still stored densely in `svm_model` for compatibility with the common prediction/inspection API, and probability calibration may materialize temporary dense pair/fold subsets. CSR prediction streams one sparse row at a time instead of expanding the full prediction matrix.

Fitted SVMs can be written with `svm_write_libsvm` and loaded with `svm_read_libsvm`. The writer emits standard LIBSVM model syntax plus optional predictor/response scaling sidecars. Multiclass pair models are expanded into a valid LIBSVM coefficient layout; duplicate support-vector rows are permitted, so the file need not be byte-identical to `e1071::write.svm` while preserving the decision functions. Files produced by the retained LIBSVM 3.23 C++ implementation are imported by the Fortran reader, and files produced by the Fortran writer have been loaded successfully by the retained C++ `svm_load_model`.

### Fuzzy clustering

`cmeans_fit`/`cmeans_fit_k` support Euclidean and Manhattan fuzzy c-means, observation weights, batch updates and the upstream on-line UFCL path. The Manhattan branch retains the upstream weighted-median prototype update.

`cshell_fit`/`cshell_fit_k` provide fuzzy c-shell memberships, radii and native nonlinear shell refinement.

The nine `fclustIndex` statistics are available through `fclust_indices`.

### Other statistical and numerical APIs

The public `e1071` module also exposes native equivalents of the computational parts of:

- discrete probability mass/CDF/quantile/random generation;
- moments, skewness and kurtosis;
- Hanning, Hamming and rectangular windows;
- STFT with radix-2 FFT and general DFT fallback;
- `rwiener` and `rbridge` simulation;
- all-pairs shortest paths and path extraction;
- mixed numeric/categorical naive Bayes;
- generalized k-NN using the translated `proxy` package;
- latent-class EM, prediction, goodness-of-fit and bootstrap;
- ICA;
- `classAgreement`, `matchClasses`, `compareMatchedClasses` and control matching;
- bagged clustering and agglomerative summaries;
- dense/CSR conversion plus LIBSVM-style CSR text I/O;
- normal and callback-defined probability-plot calculations;
- numerical imputation, scaling and interpolation;
- `bincombinations`, `permutations`, Hamming distance and column-major `element()` indexing;
- native SVM and g-KNN tuning grids.

The translated `proxy-fortran` dependency is vendored under `dependencies/proxy-fortran`, making the release self-contained.

## Deliberate interface adaptations

The Fortran library does not reproduce R formula parsing, S3 methods, data frames, plotting or namespace machinery. `hsv_palette` is graphical and is omitted.

The generic R `tune()` dispatcher can call arbitrary R functions and external packages. The native library instead supplies typed tuning APIs for translated SVM and g-KNN models. R convenience wrappers that tune `nnet`, `randomForest`, `rpart`, or `knn` are external-package orchestration and are not duplicated.

`probplot_normal` provides the normal-quantile path and `probplot_custom` accepts a typed `probplot_distribution` quantile callback. This covers the computational role of R's function-valued `qdist` without dynamic-language argument dispatch.

`svm_write_libsvm` and `svm_read_libsvm` provide LIBSVM-compatible model interchange. The multiclass writer expands independent native pair models into a valid global coefficient layout and therefore may emit duplicate support-vector rows rather than reproducing `e1071::write.svm` byte-for-byte. Optional sidecars preserve native predictor and SVR response scaling.

## Build

With FPM:

```text
fpm test
fpm run --example svm_example
```

The release is also validated directly with GNU Fortran using strict standard, warning, interface, bounds and floating-point checks.

All maintained reals use the single public `dp = real64` kind from `e1071_kinds`. Every maintained dummy argument has explicit `INTENT`/`VALUE`, is declared on its own declaration line and has a meaningful trailing FORD `!!` comment. Maintained source is free-form, stays within 132 columns and is formatted to be compatible with `fprettify`.
