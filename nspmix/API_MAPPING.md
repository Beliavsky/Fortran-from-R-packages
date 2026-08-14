# API mapping

| R nspmix | Fortran |
|---|---|
| `disc` | `make_disc`, `type(disc_dist)` |
| `dmix` | `dmix_disc` |
| `npnorm` | `make_npnorm_data` |
| `nppois` | `make_nppois_data` |
| `npgeom` | `make_npgeom_data` |
| `npnbinom` | `make_npnbinom_data` |
| `cvps` | `cvps_from_raw`, `make_cvps_data` |
| `mlogit` | `make_mlogit_data` |
| `logd` | `logd_eval` |
| `loglik` | `loglik` |
| `gridpoints` | `gridpoints` |
| `suppspace` | `support_bounds` |
| `valid` | `valid_parameters` |
| `weight` | `data_weights` |
| `grad` | `gradient_values` |
| `maxgrad` | `maxgrad` |
| `hcnm` | `hcnm` |
| `cnm(..., model="npmle")` | `cnm` |
| `cnm(..., model="proportions")` | `cnm_proportions` |
| `cnmms` | `cnmms` |
| `cnmpl` | `cnmpl` |
| `cnmap` | `cnmap` |
| `dnpnorm`, `pnpnorm`, `rnpnorm` | same names |
| `dnppois`, `pnppois`, `rnppois` | same names |
| `dnpgeom`, `pnpgeom`, `rnpgeom` | same names |
| `dnpnbinom`, `pnpnbinom`, `rnpnbinom` | same names |
| `whist` | `weighted_histogram` |
| plotting/printing S3 methods | omitted |

`initial`, `initial0`, `llex`, and `llexdb` are represented by the typed family
default routines and built-in likelihood layer rather than R-style generic
dispatch. The upstream default `llex`/`llexdb` values are zero for all built-in
families.
