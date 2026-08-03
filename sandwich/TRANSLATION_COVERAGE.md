# Translation coverage

| Upstream area | Fortran coverage |
|---|---|
| `sandwich`, `meat`, `vcovOPG` | Complete matrix-first equivalents |
| `bread`/`estfun` generic infrastructure | Generic score/bread inputs; weighted OLS adapter included |
| `vcovHC`/`meatHC` | const and HC0-HC5, including custom omega |
| `kweights` | All five upstream kernel families |
| `vcovHAC`/`meatHAC` | Explicit weights, finite-sample adjustment, diagnostics, VAR prewhitening |
| `bwAndrews`/`weightsAndrews` | AR(1) and CSS ARMA(1,1), all upstream kernels |
| `bwNeweyWest`/`NeweyWest` | Automatic bandwidth and Bartlett lag weights |
| `weightsLumley`/`weave` | Truncate and smooth WEAVE weight construction |
| `pava.blocks`/`isoacf` | Translated |
| `lrvar` | Andrews and Newey-West long-run mean covariance |
| `vcovCL`/`meatCL` | One-way/multiway HC0-HC3, `cadjust`, `multi0`, PSD fix |
| `vcovPL`/`meatPL` | Kernel/lag selection and aggregate/nonaggregate forms |
| `vcovPC`/`meatPC` | Complete-period and pairwise covariance forms |
| `vcovBS`/`vcovJK` | General replicate covariance plus OLS cluster/bootstrap methods |
| External model adapters | Not translated; use score/bread matrices |
| R object/class/formula infrastructure | Not translated |
| Plotting | No plotting code in scope |

The complete upstream R source and tests used for the port are retained under
`original/` for license attribution and traceability.
