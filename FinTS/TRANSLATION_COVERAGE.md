# Translation coverage

FinTS 0.4-9 exports 15 names in `NAMESPACE`.

| Original export | Fortran status | Fortran entry point |
|---|---|---|
| `Acf` | Computational core translated; plotting omitted | `acf`, `cross_acf` |
| `apca` | Translated | `apca` |
| `ArchTest` | Translated | `ArchTest` |
| `ARIMA` | Translated with conditional likelihood | `ARIMA` |
| `as.yearmon2` | Numeric/integer forms translated | `as_yearmon2` |
| `AutocorTest` | Translated | `AutocorTest` |
| `compoundInterest` | Translated | `compoundInterest` |
| `findConjugates` | Translated | `findConjugates` |
| `FinTS.stats` | Translated | `FinTS_stats` |
| `package.dir` | Omitted; R package/filesystem helper | - |
| `plotArmaTrueacf` | Computational core translated; plotting omitted | `plotArmaTrueacf`, `arma_true_acf` |
| `read.yearmon` | Omitted; `read.table`/`zoo` I/O wrapper | - |
| `runscript` | Omitted; R script execution helper | - |
| `simple2logReturns` | Translated | `simple2logReturns` |
| `url2data` | Omitted; network download helper | - |

The non-exported `plot.Acf` and `plot.loadings` S3 methods are plotting-only and were omitted.

Additional internal numerical facilities were implemented to keep the project dependency-free:

- pivoted dense linear solves and least squares
- symmetric Jacobi eigendecomposition
- regularized incomplete gamma and chi-square survival probabilities
- average-tie ranks
- Durand-Kerner polynomial roots
- reflection-coefficient AR/MA transforms
- multiplicative seasonal ARMA polynomial construction
- conditional Gaussian ARIMA optimization
