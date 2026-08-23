# Translation coverage

| Upstream export | Fortran | Status |
| --- | --- | --- |
| `dChern` | `dchern` | Complete |
| `pChern` | `pchern` | Complete |
| `qChern` | `qchern` | Complete |

There are no plotting routines or additional computational exports in the
supplied package.

The R package's dependency on `gsl` is used only for Airy zeros and Airy
derivatives at those zeros. Those fixed values are embedded in the Fortran
library, so there is no runtime dependency.
