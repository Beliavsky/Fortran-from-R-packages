# Translation coverage

| R export | Fortran procedure | Status |
|---|---|---|
| `HRP_Portfolio` | `HRP_Portfolio` | Translated |
| `HCAA_Portfolio` | `HCAA_Portfolio` | Translated |
| `HERC_Portfolio` | `HERC_Portfolio` | Translated |
| `DHRP_Portfolio` | `DHRP_Portfolio` | Translated |

Supporting computations included in the library:

- covariance-to-correlation distance conversion;
- Euclidean distances between correlation-distance rows;
- single, complete, average, and Ward.D2 agglomerative clustering;
- DIANA divisive clustering;
- hierarchy ordering and `cutree` labels;
- deterministic gap-statistic reference simulation;
- inverse-variance cluster portfolios;
- box-simplex projection for constrained weights.

Not translated into compiled Fortran:

- dendrogram plotting;
- R data-frame/list/S3 presentation;
- direct loading of bundled `.RData` datasets.
