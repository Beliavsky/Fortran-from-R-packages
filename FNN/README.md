# FNN-fortran

Modern Fortran/FPM translation of the computational code in the R package
`FNN` 1.1.4.1 (Fast Nearest Neighbor Search Algorithms and Applications).

The library is self-contained and uses Fortran 2018.  It provides exact
nearest-neighbor search, information-theoretic estimators, k-NN
classification/regression, and optimal weighted nearest-neighbor classifiers.
R object dispatch and printing infrastructure are not reproduced.

## Build

```text
fpm build
fpm test
fpm run --example basic
```

The source is also directly compilable with a Fortran 2018 compiler.  The
validation release was checked with gfortran 14.2 using `-Wall -Wextra -Werror
-Wimplicit-interface -fcheck=all`.

## Main interface

```fortran
use fnn
```

Arrays use the natural Fortran/statistical layout: observations are rows and
variables are columns.  Neighbor indices are 1-based, matching R.

### Neighbor search

```fortran
type(knn_result) :: z

z = get_knn(x, 10, "kd_tree")
z = get_knnx(reference, query, 10, "cover_tree")
```

Available algorithms are:

- `"kd_tree"`: exact Euclidean kd-tree search;
- `"cover_tree"`: native scaled cover hierarchy with exact subtree-radius
  branch-and-bound search;
- `"brute"`: exact linear reference search;
- `"CR"`: the upstream correlation-style distance `1 - dot_product(x,y)`.

Convenience functions `knn_index`, `knn_dist`, `knnx_index`, and `knnx_dist`
mirror the R package.

### Information measures

```fortran
h  = entropy(x, k, "kd_tree")
hc = crossentropy(x, y, k, "kd_tree")
d  = kl_divergence(x, y, k, "kd_tree")
ds = kl_dist(x, y, k, "kd_tree")
mi = mutinfo(x, y, k)
```

`entropy`, `crossentropy`, and the KL routines return one estimate for each
neighbor order `1:k`, matching FNN.  `mutinfo` implements the direct KSG
maximum-norm estimator.  `mutual_information_entropy` provides the upstream
`direct = FALSE` entropy-combination route.

### Classification and regression

```fortran
type(classification_result) :: cls
type(regression_result) :: reg

cls = knn_classify(train, test, class_id, 5, "kd_tree")
cls = knn_cv(train, class_id, 5, "cover_tree")
reg = knn_reg(train, response, 5, test=test, algorithm="kd_tree")
reg = knn_reg(train, response, 5, algorithm="kd_tree")  ! leave-one-out CV
```

Class labels are integer IDs.  The classification result always includes the
winning-class probability and neighbor index/distance matrices.

### OWNN / BNN

```fortran
type(ownn_result) :: ow

ow = ownn(train, test, class_id, k=15, algorithm="kd_tree")
ow = ownn(train, test, class_id, algorithm="kd_tree", seed=12345)
```

If `k` is omitted, five-fold CV chooses it before the OWNN/BNN weighting is
applied.  An optional `testcl=` argument computes the three test accuracies.

## `mvtnorm`

The upstream `FNN` package lists `mvtnorm` only under `Suggests`; it is used in
examples to generate multivariate-normal samples and is not called by FNN's
runtime computational code.  Consequently this port has no `mvtnorm` build
dependency.  The supplied `mvtnorm-fortran` translation can be used separately
to reproduce those simulation examples.  See `OPTIONAL_MVTNORM.md`.

## Deliberate port corrections

A few behaviors were made safer while preserving the intended mathematics:

- self-neighbor exclusion is by observation index, so duplicated observations
  do not cause the wrong zero-distance point to be discarded;
- cross-entropy/KL neighbor counts are validated against the actual reference
  sample size;
- automatic OWNN cross-validation withholds observations explicitly rather
  than masking neighbor-rank columns;
- OWNN/BNN weight dimensions are bounded to valid Fortran array ranges.

Details are in `PORTING_NOTES.md`.

## License

`FNN` declares GPL version 2 or later.  This translation is distributed under
`GPL-2.0-or-later`; the complete GPL v2 text is in `LICENSE`.

The upstream ANN 1.1.2 code embedded in FNN is copyright the University of
Maryland, Sunil Arya, and David Mount and is LGPL 2.1 or later.  Its copyright
notice is preserved in `upstream/FNN-1.1.4.1/inst/COPYRIGHTS`, and the LGPL 2.1
text is included as `ANN-LICENSE`.

The original FNN R/C/C++ sources used for this translation are retained under
`upstream/FNN-1.1.4.1/` for provenance and license auditing.
