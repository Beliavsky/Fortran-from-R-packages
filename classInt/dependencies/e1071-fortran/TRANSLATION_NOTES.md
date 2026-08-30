# Translation notes

## Scope

The translation targets the non-R computational behavior of `e1071` 1.7-17. It includes the compiled fuzzy-clustering, c-shell, Floyd-Warshall and LIBSVM kernels plus substantial numerical/statistical functionality implemented in the package's R sources.

Translated areas include:

- weighted fuzzy c-means, UFCL and fuzzy c-shell clustering;
- the `fclustIndex` family;
- C-SVC, nu-SVC, one-class SVM, epsilon-SVR and nu-SVR;
- linear, polynomial, RBF and sigmoid SVM kernels;
- SVM probabilities, decision values, scaling, class weights and cross-validation;
- latent-class analysis and bootstrap;
- naive Bayes and generalized k-NN;
- ICA;
- shortest paths;
- bagged clustering;
- matching/agreement functions;
- discrete distributions;
- moments, windows, STFT, Wiener/bridge simulation and interpolation utilities;
- sparse CSR conversion/read/write plus direct sparse-kernel SVM fitting and streamed prediction;
- normal and typed callback-defined probability-plot calculations;
- typed SVM and g-KNN tuning grids.

R formulas, S3 printing/summary/plotting, data-frame/list dispatch and graphical utilities are language-interface concerns and are not translated literally.

## LIBSVM translation

The upstream package contains LIBSVM 3.23 C++ code. The Fortran solver follows the same dual optimization structure and supports the five upstream SVM modes and four kernels.

The native implementation deliberately uses `real(dp)` throughout to satisfy the package-wide real-kind rule. Upstream LIBSVM stores cached kernel rows in single-precision `Qfloat`; the Fortran solver instead uses a dense double-precision kernel matrix and does not reproduce LIBSVM's cache/shrinking execution path. As a result, support-vector coefficients and rho can differ slightly while remaining within the solver tolerance and producing the same classifications on parity fixtures.

For a retained deterministic RBF C-SVC fixture, direct compilation of the upstream LIBSVM source gives decision magnitudes essentially equal to one. The Fortran solver reproduces the labels and decision magnitudes within 0.02. The sign convention can differ because the Fortran multiclass layer sorts integer class labels before constructing pair models.

Probability fitting uses the LIBSVM five-fold calibration structure and native Platt sigmoid fitting/multiclass coupling. Fold shuffling uses the deterministic native `rng_state` (`probability_seed`) rather than R's current RNG state, so probability coefficients are not expected to be bit-for-bit R-identical for a given R session.

SVR probability scale uses five-fold residuals and the upstream outlier rule based on `5*sqrt(2)*MAE`.

## Exact compiled-kernel fixture

The translated Euclidean fuzzy c-means kernel was compared directly against the retained upstream `src/cmeans.c` on a deterministic eight-row fixture. It reproduces the upstream iteration count, final centers, memberships and objective to the printed double-precision digits. The retained objective fixture is:

```text
0.0343428119711938
```

The Manhattan/weighted-median path was also compared against the same upstream C kernel. It reproduces the one-iteration stop, unchanged weighted-median centers and objective:

```text
0.21844433535388025
```

These values are enforced by `test/test_fuzzy.f90`.

## Native API adaptations

### Fuzzy clustering

The R APIs accept a number of clusters and randomly choose distinct initial centers. `cmeans_fit_k` and `cshell_fit_k` provide that behavior using the package's native deterministic RNG. Lower-level APIs also accept explicit initial centers for deterministic parity work.

### g-KNN and matching

The translated `proxy-fortran` library supplies the proximity calculations used by generalized k-NN and numeric control matching. Mixed-variable control matching uses native Gower similarity with explicit per-column type codes rather than an R data frame plus `cluster::daisy` dispatch.

### Sparse matrices

`read_matrix_csr`, `write_matrix_csr` and dense/CSR conversion are native. CSR SVM fitting now builds kernel matrices directly from sparse rows and does not materialize the full dense training matrix. This also fixes an important parity point: upstream `e1071::svm` forcibly disables scaling for `matrix.csr`, so the Fortran CSR entry points now ignore the compatibility `scale_mask` values and use identity predictor/response scaling.

The common `svm_model` still stores the selected support vectors densely, keeping dense and sparse fits interoperable with the same prediction, inspection, and serialization routines. Probability calibration is optional and may materialize temporary dense pair/fold subsets. The dual solver still uses dense kernel/Hessian working matrices rather than LIBSVM's kernel cache, so very-large-`n` sparse-memory complexity is not claimed. CSR prediction expands only one query row at a time.

### Tuning

R's generic `tune()` can dynamically invoke arbitrary methods, including functions from other R packages. Static Fortran cannot reproduce that unrestricted dynamic interface directly. Native typed tuning grids are provided for SVM classification/regression and g-KNN classification/regression. User-supplied fold IDs are supported; default folds are deterministic rather than R-RNG-driven.

### Probability plots

`probplot_normal` implements the original normal path. `probplot_custom` accepts a typed `probplot_distribution` containing a procedure pointer to a scalar quantile function, so arbitrary reference distributions can be used from Fortran without reproducing R's dynamic `...` argument machinery.

### Serialization

`svm_write_libsvm` emits standard LIBSVM model syntax and optional two-column predictor/response scaling sidecars. The native one-versus-one representation is converted into LIBSVM's global coefficient rows by emitting each pair's support-vector contribution in the appropriate class block. This may duplicate the same geometric support vector across pair models, but the serialized decision functions are identical.

`svm_read_libsvm` imports standard C-SVC, nu-SVC, one-class, epsilon-SVR and nu-SVR files for the four kernels supported by `e1071`. A retained LIBSVM 3.23 C++ multiclass model is included as a test fixture; the Fortran reader reproduces its training predictions. Conversely, a Fortran-emitted multiclass model was loaded by the retained C++ `svm_load_model` and reproduced the expected class predictions. Byte-for-byte equivalence with `e1071::write.svm` is not claimed because support-vector deduplication/layout and text formatting can differ.

## Validation

The strict suite covers:

- discrete distributions, moments and utility combinatorics;
- STFT and shortest paths;
- CSR round trips, sparse SVM dense-parity checks, and normal/custom probability-plot calculations;
- Euclidean and Manhattan fuzzy c-means parity fixtures;
- c-shell and fuzzy-cluster indices;
- naive Bayes, g-KNN, LCA, ICA, matching and bagged clustering;
- mixed Gower control matching;
- all five SVM modes;
- all four SVM kernels;
- SVC/SVR probability paths;
- dense and direct sparse-kernel CSR SVM entry points;
- LIBSVM model import/export round trips, including an upstream-C++ fixture;
- SVM cross-validation and tuning.

No R executable is installed in the translation environment. Direct reference comparisons for compiled kernels were therefore performed by compiling retained upstream C/C++ sources with minimal R-runtime stubs where practical.

## Source conventions

- `dp = real64` is defined exactly once in `e1071_kinds` and re-exported by `e1071`.
- Every maintained dummy argument has explicit `INTENT` or `VALUE`.
- Every dummy argument is declared on its own declaration line.
- Every dummy declaration has a meaningful trailing FORD `!!` documentation comment.
- Maintained source is free-form, stays within the normal 132-column limit and is formatted to be compatible with `fprettify`.
- `fprettify` is not installed in the release-validation environment, so compatibility is enforced by source style/audit rather than by running that executable.
