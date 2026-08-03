# API coverage map

## Direct high-level counterparts

| R package interface | Fortran interface |
|---|---|
| `highs_solve()` | generic `highs_solve()` |
| `highs_model()` | `highs_model_from_dense()` and `highs_model` |
| `highs_control()` | `highs_control` |
| `highs_solver()` | `highs_solver` plus persistent procedures |
| dense or `Matrix` constraint input | dense arrays, CSC/CSR, or triplets |
| `Q` quadratic objective | `q=` dense Hessian or `model%q` sparse Hessian |
| `start` | `start=` or `highs_set_start()` |

## Model construction

- `hi_new_model`: default construction of `highs_model`
- `hi_model_set_ncol`, `hi_model_set_nrow`: `num_col`, `num_row`
- objective, bounds, sense, and offset setters: corresponding `highs_model` components
- constraint matrix: `highs_sparse_matrix`
- Hessian: `highs_hessian_from_dense` or direct sparse assignment
- variable types: `integrality`

Fortran favors complete typed model construction over a long sequence of mutating setters.

## Persistent solver functions

Implemented:

- create/destroy solver;
- pass LP/MIP/QP model and Hessian;
- run and presolve;
- clear model, solver state, basis, and options;
- set Boolean, integer, double, and string options;
- obtain model status, status text, objective, solution vectors, run time, and selected information fields;
- change costs, variable bounds, row bounds, coefficients, integrality, objective sense, and offset;
- set primal/dual start values;
- get/set basis;
- get primal and dual rays;
- read and write models;
- query version and infinity.

## Adapted or omitted interfaces

The following are not currently exposed one-for-one:

- R S3/R6-like object environments, printing, `checkmate` validation, names, and list return values;
- callback logging into the R console;
- option-table enumeration and dynamic R option typing;
- adding rows/columns incrementally;
- ranged sensitivity records returned by `getRanging`;
- IIS extraction and writing;
- presolved-model extraction and explicit postsolve objects;
- sparse partial-solution setters;
- low-level row/column slice getters;
- global scheduler reset;
- R `Matrix` classes and formula/data-frame handling.

These are backend capabilities rather than missing optimization algorithms. The bridge ABI is deliberately small and can be extended without changing the high-level model API.
