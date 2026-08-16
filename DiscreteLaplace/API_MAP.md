# API map

| R export | Fortran mapping |
|---|---|
| `Edlaplace` | `edlaplace` -> `type(edlaplace_result)` |
| `Edlaplace2` | `edlaplace2` -> `type(edlaplace2_result)` |
| `ddlaplace` | `ddlaplace` |
| `pdlaplace` | `pdlaplace` |
| `qdlaplace` | `qdlaplace` |
| `rdlaplace` | `rdlaplace` |
| `ddlaplace2` | `ddlaplace2` |
| `palaplace2` | `palaplace2` |
| `pdlaplace2` | `pdlaplace2` |
| `qdlaplace2` | `qdlaplace2` |
| `rdlaplace2` | `rdlaplace2` |
| `dlaplacelike2` | `dlaplacelike2` |
| `estdlaplace` | `estdlaplace` -> `type(estdlaplace_result)` |
| `estdlaplace2` | `estdlaplace2` |
| `iFI` | `ifi` |
| `iFI2` | `ifi2` |
| `ioFI2` | `iofi2` |
| `loss` | `loss` |

All 18 names selected by the upstream `exportPattern("^[[:alpha:]]+")` have a
computational mapping.  There is no plotting code in the supplied package.
