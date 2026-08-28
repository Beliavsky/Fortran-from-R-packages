# API map

| R / legacy API | Fortran API | Status |
|---|---|---|
| `pbivnorm(x,y,rho)` | elemental `pbivnorm(x,y,rho)` | complete numerical port |
| R vector recycling | `pbivnorm_recycle` | supported |
| `MVBVU` | `bvn_upper_tail` | direct modern translation |
| `.Fortran("PBIVNORM",...)` | not needed | replaced by native module API |

The package has no fitting, simulation, plotting, or other exported numerical
routines beyond `pbivnorm`.
