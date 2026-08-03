# API map

| R function | Fortran procedure | Status |
|---|---|---|
| `insert.at` | `insert_at` | translated for scalar insertion |
| `po.avg` | `po_avg` | LASSO, RIDGE, and EN |
| `po.grossExp` | `po_gross_exp` | NOSHORT and EQUAL |
| `po.covShrink` | `po_cov_shrink` | translated using vendored corpcor |
| `po.cols` | `po_cols` | translated |
| `po.JM` | `po_jm` | translated |
| `buh.clust` | `buh_clust` | translated |
| `po.bhu` | `po_bhu` | translated |
| `po.TZT` | `po_tzt` | translated |
| `po.SW` | `po_sw` | translated |
| `po.SW.lasso` | `po_sw_lasso` | translated |
| `prepare_data` | `prepare_data` | typed date/return interface |
| `perform_analysis` | `perform_analysis` | numerical outputs; no plots |
| `ren` | `ren_run` | translated |
| `setup_parallel` | none | R runtime infrastructure, skipped |

All public procedures are re-exported by module `ren`.
