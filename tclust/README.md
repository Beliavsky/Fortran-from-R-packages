# tclust — modern Fortran translation

This package translates the computational core of the R package [`tclust`](https://github.com/valentint/tclust) to modern free-form Fortran with FPM packaging.

Implemented functionality includes robust trimmed Gaussian clustering (`tclust`), HARD and MIXT objectives, eigenvalue and determinant/shape covariance restrictions, trimmed k-means, robust clustering around affine subspaces, discriminant-factor diagnostics, CTL curves, CLA/BIC/ICL grids, simulation helpers, and Rand/Fowlkes–Mallows agreement indices.

The implementation is self-contained and does not require BLAS, LAPACK, R, Rcpp, or Armadillo.

## Build and test

With FPM and gfortran installed:

```text
fpm build
fpm test
fpm run --example tclust_example
```

The public module is `tclust_mod` and re-exports the package's single `dp` real kind.

A minimal fit looks like:

```fortran
use tclust_mod, only : dp, tclust_result, tclust

type(tclust_result) :: fit
real(dp) :: x(100, 2)

! Fill x, observations in rows.
call tclust(x, 3, fit, alpha=0.10_dp, restr_fact=12.0_dp, seed=123)
```

`fit%cluster == 0` denotes trimmed observations. Cluster centers are stored as `fit%centers(p,k)`, covariance matrices as `fit%cov(p,p,k)`, and posterior probabilities as `fit%posterior(n,k)`.

## Translation coverage

**Package status: substantial. Coverage: 11 of 11 (100%).**

The fraction measures mapped exported computational functions rather than complete R compatibility. The R package's GPCM covariance-model branch, R center/scale callbacks, parallel execution, S3/list construction, plotting/presentation, and exact R RNG streams are not reproduced. `tclustICsol` has a useful but simplified solution-ranking heuristic.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
|---|---|---|---|
| `simula.rlg` | `tclust_mod` | `simulate_rlg` | substantial |
| `simula.tclust` | `tclust_mod` | `simulate_tclust` | substantial |
| `rlg` | `tclust_mod` | `rlg`, `rlg_refine` | substantial |
| `tclust` | `tclust_mod` | `tclust`, `tclust_refine`, `tclust_information_criteria` | substantial |
| `DiscrFact` | `tclust_mod` | `discr_fact` | substantial |
| `ctlcurves` | `tclust_mod` | `ctlcurves` | substantial |
| `tkmeans` | `tclust_mod` | `tkmeans`, `tkmeans_refine` | substantial |
| `tclustIC` | `tclust_mod` | `tclust_ic` | substantial |
| `tclustICsol` | `tclust_mod` | `tclust_ic_solutions` | partial |
| `randIndex` | `tclust_mod` | `rand_index_labels`, `rand_index_table` | complete |
| `FowlkesMallowsIndex` | `tclust_mod` | `fowlkes_mallows_labels`, `fowlkes_mallows_table` | complete |

See `API_COVERAGE.md` and `[extra.translation]` in `fpm.toml` for details.

## Validation

`test/test_tclust.f90` provides deterministic regression tests. `tools/parity_driver.f90` plus `tools/check_parity.py` independently validate major numerical outputs with NumPy/SciPy. `VALIDATION.md` records the validation performed for the packaged checkpoint.

## License and provenance

The upstream R package is GPL-3. See `LICENSE`, `NOTICE.md`, and `upstream/` for license and provenance information.
