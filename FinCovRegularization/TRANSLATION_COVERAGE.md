# Translation coverage

| Original R export | Modern Fortran procedure | Status |
|---|---|---|
| `F.norm2` | `f_norm2` | translated |
| `O.norm2` | `o_norm2` | translated |
| `FundamentalFactor.Cov` | `fundamental_factor_cov` | translated |
| `GMVP` | `gmvp` | translated |
| `Ind.Cov` | `ind_cov` | translated |
| `MacroFactor.Cov` | `macro_factor_cov` | translated |
| `RiskParity` | `risk_parity` | translated |
| `StatFactor.Cov` | `stat_factor_cov` | translated |
| `banding` | `banding` | translated |
| `banding.cv` | `banding_cv` | translated |
| `hard.thresholding` | `hard_thresholding` | translated |
| `soft.thresholding` | `soft_thresholding` | translated |
| `tapering` | `tapering` | translated |
| `tapering.cv` | `tapering_cv` | translated |
| `threshold.cv` | `threshold_cv` | translated |
| `threshold.min` | `threshold_min` | translated |

The non-computational S3 methods `plot.CovCv`, `print.CovCv`, and
`summary.CovCv` are omitted. The binary R sample dataset is preserved in
`original/data/` but no R-data reader is included.
