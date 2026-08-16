# API map

| R API | Fortran API | Status |
|---|---|---|
| `ddiweibull` | `ddiweibull` | translated |
| `pdiweibull` | `pdiweibull` | translated; optional upper/log tails added |
| `qdiweibull` | `qdiweibull` | translated; optional upper/log probabilities added |
| `rdiweibull` | `rdiweibull` | translated as array-filling subroutine |
| `Ediweibull` | `ediweibull` | translated; returns `diw_moments` |
| `hrdiweibull` | `hrdiweibull` | translated |
| `ahrdiweibull` | `ahrdiweibull` | translated |
| `loglikediw` | `loglikediw` | translated |
| `lossdiw` | `lossdiw` | translated |
| `heuristic` | `heuristic` | translated |
| `estdiweibull` | `estdiweibull` | translated for `P`, `M`, `H`, `PP` |

The R package uses `exportPattern(".")`; these are all computational functions in the supplied source tree.
