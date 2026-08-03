# Translation coverage

| Upstream export | Fortran procedure | Status |
|---|---|---|
| `SV` | `sv` | Implemented |
| `SVJ` | `svj` | Implemented |
| `SV1F` | `sv1f` | Implemented |
| `SV1FJ` | `sv1fj` | Implemented with documented jump-call correction |
| `SV2F` | `sv2f` | Implemented |
| `jumptestday` | `jumptestday` | Implemented |
| `jumptestperiod` | `jumptestperiod` | Implemented |
| `pcombine` | `pcombine` | Implemented |
| `ppool` | `ppool` | All six methods implemented |

Internal R routines `bns`, `amin`, and `amed` are public in Fortran as
`bns_statistic`, `amin_statistic`, and `amed_statistic`. Native recursions `lp`,
`lp2`, `pvc`, `pvc0`, and `pv2` are represented by `lp_path`, `pvc_path`, and
`pv2_path`; zero-jump and jump variants share the same typed path routine.

No plotting or data-download code exists in the upstream package. R S4 classes,
Rcpp registration, and package-loading infrastructure are not compiled.
