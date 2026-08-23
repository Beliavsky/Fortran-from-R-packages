# Translation coverage

| R export | Fortran coverage |
|---|---|
| `dsgt` | `dsgt`, `sgt_pdf`, `sgt_logpdf` |
| `psgt` | `psgt`, `sgt_cdf` |
| `qsgt` | `qsgt`, `sgt_quantile` |
| `rsgt` | `rsgt` |
| `sgt.mle` | `sgt_mle_constant`, `sgt_mle_model` |

The `print`, `summary`, and S3 print methods are presentation machinery rather
than numerical exports. Their numerical summary quantities (standard errors,
z scores, and normal-approximation p-values) are fields of `sgt_mle_result`.

The R formula parser itself is not translated. Arbitrary formula behavior is
represented numerically by the `sgt_observation_model` callback interface.
