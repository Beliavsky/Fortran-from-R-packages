# API coverage

This translation uses the package's exported computational R functions as its coverage basis. Plotting, printing, summaries, data objects, internal helpers, and presentation-only methods are excluded.

**Status: substantial — 11 of 11 exported computational functions mapped (100%).**

The percentage counts functions with meaningful Fortran mappings; it does **not** mean exact R-interface or numerical-path compatibility. In particular, the GPCM covariance-model branch of `tclust`, R callback/formula-like scaling interfaces, parallel execution, S3/list presentation, and exact R RNG streams are not reproduced.

| R function | Fortran module | Fortran procedure(s) | Status | Main difference |
|---|---|---|---|---|
| `simula.rlg` | `tclust_mod` | `simulate_rlg` | substantial | Same scenario; different RNG stream. |
| `simula.tclust` | `tclust_mod` | `simulate_tclust` | substantial | Same type 1/2, k=3/6 scenarios; different RNG stream. |
| `rlg` | `tclust_mod` | `rlg`, `rlg_refine` | substantial | Affine/PCA trimmed clustering translated; R scaling and parallel/UI behavior omitted. |
| `tclust` | `tclust_mod` | `tclust`, `tclust_refine`, `tclust_information_criteria` | substantial | HARD/MIXT plus eigen/determinant restrictions translated; GPCM omitted. |
| `DiscrFact` | `tclust_mod` | `discr_fact` | substantial | Numerical discriminant-factor calculation translated; plots/S3 omitted. |
| `ctlcurves` | `tclust_mod` | `ctlcurves` | substantial | Numerical grid translated; R parallel/progress/dots interface omitted. |
| `tkmeans` | `tclust_mod` | `tkmeans`, `tkmeans_refine` | substantial | Trimmed k-means translated; R scaling/parallel/drop-empty interface differs. |
| `tclustIC` | `tclust_mod` | `tclust_ic` | substantial | CLA/BIC/ICL grid translated; all three are computed together. |
| `tclustICsol` | `tclust_mod` | `tclust_ic_solutions` | partial | Adjacent stability and criterion ranking translated; `findBestSolutions` heuristic simplified. |
| `randIndex` | `tclust_mod` | `rand_index_labels`, `rand_index_table` | complete | Integer labels and contingency tables supported. |
| `FowlkesMallowsIndex` | `tclust_mod` | `fowlkes_mallows_labels`, `fowlkes_mallows_table` | complete | Integer labels and contingency tables supported. |

See `fpm.toml` for machine-readable per-function mappings and notes.
