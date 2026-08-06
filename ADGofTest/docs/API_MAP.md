# API map

| Upstream R routine | Fortran routine | Notes |
|---|---|---|
| `ad.test.statistic(x)` | `ad_statistic(probabilities, ...)` | Sorts internally and computes the same statistic. |
| `ad.test.pvalue(x, n)` | `ad_distribution_cdf(statistic, n, ...)` | Despite its R name, the upstream helper returns the lower-tail CDF. |
| `1 - ad.test.pvalue(...)` | `ad_p_value(statistic, n, ...)` | Returns the upper-tail p-value. |
| `ad.test(x)` | `ad_test_uniform(x, result, ...)` | Tests values already on the uniform scale. |
| `ad.test(x, distr.fun, ...)` | `ad_test_with_cdf(x, cdf, result, ...)` or generic `ad_test` | Extra distribution parameters are supplied through a Fortran wrapper procedure. |

The Fortran result type contains `statistic`, `p_value`, `n`, `status`, and
`message` fields instead of an R `htest` object.
