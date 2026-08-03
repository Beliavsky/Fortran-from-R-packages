# Testing

The test suite contains five independent programs:

- `test_statistics`: robust scales, chi-square inversion, matrix square root,
  correlation conversion, and ILR transformation;
- `test_covariance_estimators`: classical, MCD, MVE, OGK, M, S, MM, SDE, and
  MRCD estimators under gross contamination;
- `test_pca`: all translated PCA families and distance outputs;
- `test_discriminant`: classical and robust LDA/QDA fitting and prediction;
- `test_multivariate_tests`: one-/two-sample Hotelling tests and Wilks MANOVA.

Run with `fpm test` or `scripts/build_checked.sh`.

The checked script uses GNU Fortran options including `-std=f2018`,
`-Wall`, `-Wextra`, `-Werror`, `-fcheck=all`, and `-fimplicit-none`.
