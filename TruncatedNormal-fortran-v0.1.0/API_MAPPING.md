# API mapping

The computational exports from the upstream `NAMESPACE` are represented as
follows.  Internal R/C++ symbols beginning with a dot are exposed through the
corresponding public Fortran routines rather than as duplicate aliases.

| R export | Fortran |
| --- | --- |
| `.cholpermGB` | `cholperm(..., method='GB')` |
| `.cholpermGGE` | `cholperm(..., method='GGE')` |
| `.dmvnorm_arma` | `dmvnorm` |
| `.dmvt_arma` | `dmvt` |
| `Phinv` | `phinv` |
| `cholperm` | `cholperm` |
| `dtmvnorm` | `dtmvnorm` |
| `dtmvt` | `dtmvt` |
| `lnNpr` | `lnNpr` |
| `mvNcdf` | `mvncdf` |
| `mvNqmc` | `mvnqmc` |
| `mvTcdf` | `mvtcdf` |
| `mvTqmc` | `mvtqmc` |
| `mvrandn` | `mvrandn` |
| `mvrandt` | `mvrandt` |
| `norminvp` | `norminvp` |
| `pmvnorm` | `pmvnorm` |
| `pmvt` | `pmvt` |
| `ptmvnorm` | `ptmvnorm` |
| `ptmvt` | `ptmvt` |
| `qtnorm` | `qtnorm_vec` |
| `rtmvnorm` | `rtmvnorm` |
| `rtmvt` | `rtmvt` |
| `rtnorm` | `rtnorm_vec` |
| `trandn` | `trandn` |
| `tregress` | `tregress` / `tregress_result` |
