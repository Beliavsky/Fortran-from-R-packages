# API map

All exported computational functions in upstream `NAMESPACE` have numerical
counterparts.

| R export | Fortran counterpart |
|---|---|
| `dBER`, `pBER`, `qBER`, `rBER` | `dber`, `pber`, `qber`, `rber` |
| `dBER2`, `pBER2`, `qBER2`, `rBER2` | `dber2`, `pber2`, `qber2`, `rber2` |
| `dUHLG`, `pUHLG`, `qUHLG`, `rUHLG` | `duhlg`, `puhlg`, `quhlg`, `ruhlg` |
| `dUMB`, `pUMB`, `qUMB`, `rUMB` | `dumb`, `pumb`, `qumb`, `rumb` |
| `dUPHN`, `pUPHN`, `qUPHN`, `rUPHN` | `duphn`, `puphn`, `quphn`, `ruphn` |
| `BER()` | `zod_ber` + `zero_one_families` callbacks |
| `BER2()` | `zod_ber2` + `zero_one_families` callbacks |
| `UHLG()` | `zod_uhlg` + `zero_one_families` callbacks |
| `UMB()` | `zod_umb` + `zero_one_families` callbacks |
| `UPHN()` | `zod_uphn` + `zero_one_families` callbacks |

The R family constructors return S3 lists containing closures. In Fortran those
closures are represented by family-ID dispatch routines rather than dynamic R
objects. `fit_zero_one` supplies a direct regression interface using model
matrices rather than R formulas.

Randomized-quantile-residual plotting (`rqres.plot`) is presentation code and
is omitted. `family_cdf` exposes the numerical CDF needed to construct those
residuals externally.
