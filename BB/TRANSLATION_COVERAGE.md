# Translation coverage

Upstream: R package `BB`, version 2026.1.0.

| R source | Fortran translation | Status |
|---|---|---|
| `R/spg.R` | `src/bb_spg.f90` | Translated |
| `R/sane.R` | `src/bb_nonlinear.f90` | Translated |
| `R/dfsane.R` | `src/bb_nonlinear.f90` | Translated |
| `R/BBoptim.R` | `src/bb_drivers.f90` | Translated |
| `R/BBsolve.R` | `src/bb_drivers.f90` | Translated |
| `R/multiStart.R` | `src/bb_drivers.f90` | Translated |
| `R/project.R` | `src/bb_projection.f90` | Translated using supplied quadprog port |
| `stats::optim(..., Nelder-Mead)` | `src/bb_aux_optim.f90` | Standalone replacement |
| `stats::optim(..., L-BFGS-B)` | `src/bb_aux_optim.f90` | Standalone dense-BFGS replacement |
| `numDeriv::grad` gradient diagnostic | n/a | Omitted; not part of optimization iteration |
| R `try`, lists, attributes, printing | Fortran result/control types | Adapted |
| Vignette/demo plotting | n/a | Omitted as requested |

## API adaptations

R's dynamically typed callbacks and `...` arguments become explicit Fortran
procedure interfaces.  Additional user data can be supplied through host
association in internal procedures or through module data.

`projectLinear` keeps the R convention that constraints are rows of `A` and
are `A*x >= b`, with the first `meq` rows equalities.  Internally the matrix is
transposed before calling the supplied Fortran `quadprog` implementation.

The R `multiStart` input has one start per row; the Fortran API uses one start
per column (`starts(:,k)`) for natural contiguous Fortran storage.

## Dependency

The following source modules are copied unchanged from the attached
`quadprog-fortran` translation:

- `src/quadprog_kinds.f90`
- `src/quadprog_core.f90`
- `src/quadprog.f90`

They retain their GPL-2.0-or-later SPDX/license terms.  BB translation files
carry GPL-3.0-only SPDX identifiers, consistent with the upstream BB package.
