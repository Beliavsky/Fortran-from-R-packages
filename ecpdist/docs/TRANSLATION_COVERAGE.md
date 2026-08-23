# Translation coverage

The upstream package exports every alphabetically named R function through
`exportPattern`. Its substantive numerical functions are all covered.

| Upstream function | Coverage |
|---|---|
| `decp` | complete |
| `pecp` | complete |
| `qecp` | complete; option semantics corrected |
| `recp` | complete |
| `secp` | complete |
| `hecp` | complete |
| `ecp_kmoment` | complete |
| `ecp_kmoment_cond` | complete |
| `ecp_mrl` | complete |
| `ecp_shape` | complete |
| `ecp_plot` | omitted: graphics only |

There are no remaining numerical algorithms in the supplied R package that
require translation for computational parity.
