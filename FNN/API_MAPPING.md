# API mapping

| R FNN interface | Fortran interface | Notes |
|---|---|---|
| `get.knn` | `get_knn` | Returns `knn_result` |
| `knn.index` | `knn_index` | 1-based indices |
| `knn.dist` | `knn_dist` | Euclidean distances |
| `get.knnx` | `get_knnx` | Reference/query search |
| `knnx.index` | `knnx_index` | 1-based indices |
| `knnx.dist` | `knnx_dist` | Distances |
| `entropy` | `entropy` | Kozachenko-Leonenko-style vector for orders `1:k` |
| `crossentropy` | `crossentropy` | Vector for orders `1:k` |
| `KL.divergence` | `kl_divergence` | Vector for orders `1:k` |
| `KL.dist` | `kl_dist` | Symmetric KL vector |
| `KLx.divergence` | `klx_divergence` | Same exact estimator, shared native kernel |
| `KLx.dist` | `klx_dist` | Same exact symmetric estimator, shared native kernel |
| `mutinfo(..., direct=TRUE)` | `mutinfo` | Direct KSG maximum-norm estimator |
| `mutinfo(..., direct=FALSE)` | `mutual_information_entropy` | Entropy-combination vector |
| `knn` | `knn_classify` | Integer class IDs replace R factors |
| `knn.cv` | `knn_cv` | Leave-one-out classification |
| `knn.reg` | `knn_reg` | Optional `test=`; absent means LOOCV |
| `ownn` | `ownn` | KNN, OWNN and BNN predictions; optional automatic k |

## Native registered routines

The R package's registered C/C++ entry points are covered as follows:

- `get_KNN_kd`, `get_KNNX_kd`: `get_knn/get_knnx(...,"kd_tree")`;
- `get_KNN_cover`, `get_KNNX_cover`: `...("cover_tree")`;
- `get_KNN_brute`, `get_KNNX_brute`: `...("brute")`;
- `get_KNN_CR`, `get_KNNX_CR`: `...("CR")`;
- `KNN_MLD_kd`, `KNN_MLD_brute`: `mean_log_knn_distance`;
- `KL_divergence`, `KL_dist`: `klx_divergence`, `klx_dist`;
- `mutinfo`, `mdmutinfo`: the dimension-general `mutinfo` implementation.

The old R `.C` ABI, R factor/class metadata, S3 print methods, and R memory
management are intentionally not replicated.
