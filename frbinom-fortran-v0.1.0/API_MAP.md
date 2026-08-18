# API map

`frbinom` 1.0.0 exports eight computational functions. All are translated.

| R API | Fortran API |
|---|---|
| `dfrbinom` | `dfrbinom`, `dfrbinom_vec` |
| `pfrbinom` | `pfrbinom`, `pfrbinom_vec` |
| `qfrbinom` | `qfrbinom`, `qfrbinom_vec` |
| `rfrbinom` | `rfrbinom` |
| `dfrbinom2` | `dfrbinom2`, `dfrbinom2_vec` |
| `pfrbinom2` | `pfrbinom2`, `pfrbinom2_vec` |
| `qfrbinom2` | `qfrbinom2`, `qfrbinom2_vec` |
| `rfrbinom2` | `rfrbinom2` |

The shared finite-state kernels are also public:

- `frbinom_pmf_table`
- `frbinom_cdf_table`
- `frbinom2_pmf_table`
- `frbinom2_cdf_table`

These avoid recomputing the generalized-Bernoulli-process recurrence when many
probabilities or quantiles are needed.
