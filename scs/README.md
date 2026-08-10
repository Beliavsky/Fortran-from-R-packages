# scs-fortran

`scs-fortran` is a modern Fortran translation of the computational core of the
R package **scs 3.2.7**, which embeds the **SCS (Splitting Conic Solver)**
algorithm. The project uses FPM and has no required external numerical
libraries.

Version **0.2.0** replaces the v0.1.0 dense KKT factorization with a native
Fortran sparse QDLDL backend. The public optimization API is unchanged.

The translated solver handles convex quadratic cone programs of the form

```text
minimize    (1/2) x' P x + c' x
subject to  A x + s = b
            s in K
```

where `K` may combine zero, nonnegative, box, second-order, positive
semidefinite, exponential, dual-exponential, and power cones.

## Build

```sh
fpm build
fpm test
fpm run --example basic_lp
```

The source is Fortran 2018 and is also tested directly with GNU Fortran using
`-std=f2018 -fcheck=all`.

## Main API

The principal public routine is

```fortran
call scs(data, cone, settings, solution, info)
```

from `scs_solver`.

The public types are in `scs_types`:

- `scs_matrix` -- compressed sparse column (CSC) matrix.
- `scs_data` -- `A`, optional upper-triangular `P`, `b`, and `c`.
- `scs_cone` -- cone description.
- `scs_settings` -- solver controls corresponding to the R package controls.
- `scs_solution` -- `x`, `y`, and `s`.
- `scs_info` -- status, residuals, objectives, timing, iteration statistics,
  and sparse-linear-system statistics.

Unlike the C API, CSC indices are **1-based** to follow normal Fortran
conventions. `A%p` contains column starts in the range `1 ... nnz+1`, and
`A%i` contains row numbers in `1 ... m`. `dense_to_csc` and
`dense_upper_to_csc` are supplied for convenient setup. For `P`, only the upper
triangle is stored, as in SCS.

The cone row order follows SCS:

1. zero cone (`z`)
2. nonnegative cone (`l`)
3. box cone (`bsize`, `bl`, `bu`)
4. second-order cones (`q`)
5. positive-semidefinite cones (`s`)
6. primal exponential cones (`ep`)
7. dual exponential cones (`ed`)
8. power cones (`p`)

PSD entries use the same SCS `svec` scaling convention: off-diagonal matrix
entries are multiplied by `sqrt(2)` in the packed vector.

## What was translated

The Fortran implementation includes the numerical parts of the bundled SCS
solver:

- homogeneous self-dual embedding and Douglas-Rachford iteration;
- primal/dual residuals and certificates for infeasibility/unboundedness;
- quadratic objectives;
- Ruiz/L2 data equilibration and adaptive dual scaling;
- zero, linear, box, SOC, PSD, exponential, dual-exponential, and power cone
  projections;
- Anderson acceleration (type I and type II modes);
- sparse CSC matrix-vector products;
- sparse KKT assembly and native QDLDL factorization;
- warm starts and time limits;
- solver status/objective/residual/timing reporting.

The PSD projection uses a self-contained Jacobi symmetric eigensolver. This
keeps FPM builds dependency-free.

## v0.2.0 sparse direct backend

The upstream C solver uses SuiteSparse AMD ordering followed by QDLDL. v0.2.0
ports the **QDLDL elimination-tree, numeric LDL^T factorization, and triangular
solves to native Fortran** and builds the same upper-triangular quasi-definite
KKT matrix directly in CSC form.

Important properties of the new backend:

- KKT storage is sparse CSC; no `(n+m) x (n+m)` dense matrix is formed.
- `L` is stored sparsely in CSC form.
- symbolic elimination-tree analysis is performed once and reused when SCS
  changes only the diagonal regularization during adaptive scaling.
- repeated numeric factorizations reuse all symbolic and workspace storage.
- duplicate numerical entries in a KKT column are accumulated safely by the
  translated QDLDL numeric routine.
- `scs_info` reports `kkt_nnz`, `factor_nnz`, `factorizations`, and
  `symbolic_analyses`.

The current ordering is the natural KKT ordering (`x` variables followed by
constraint variables), so `lin_sys_solver` reports
`native-sparse-qdldl-natural`.

### Remaining difference from upstream SCS

**AMD ordering is not yet translated.** Therefore v0.2.0 has the sparse memory
and arithmetic structure of QDLDL, but fill-in can be substantially larger than
upstream SCS on difficult sparse graphs. A future release can add a native AMD
(or another fill-reducing symmetric ordering) without changing the public SCS
API or the QDLDL module.

The optional MKL/GPU/indirect linear-system backends and the R
`.Call`/S3/Matrix/slam conversion layer remain omitted. CSV iteration logging
and serialized problem-data output are diagnostics/interface facilities and are
not included. There was no substantive plotting algorithm to translate.

## Tests

`test/test_scs.f90` covers:

- the package README/basic equality example;
- an unbounded problem;
- the package SOCP example;
- a convex quadratic program;
- the package's semidefinite-program regression case;
- Anderson-accelerated SOCP solving;
- nontrivial exponential, dual-exponential, primal/dual power, and box cone
  projection checks.

`test/test_sparse_ldlt.f90` additionally checks the v0.2.0 backend directly:

- sparse QDLDL solve residuals against the assembled KKT matrix;
- sparse KKT/factor nonzero accounting;
- symbolic-analysis reuse across diagonal refactorization.

The PSD regression solution and objectives are checked against the values in
`inst/tinytest/test_psd.R` from the original R package.

## Benchmark

`benchmark/benchmark_sparse_ldlt.f90` compares the v0.2.0 sparse backend with a
copy of the v0.1.0 dense LDL^T algorithm on a banded synthetic KKT problem.
It is intentionally outside FPM's automatic example/test targets.

An illustrative GNU Fortran 14.2.0 `-O2` run in the translation environment,
for KKT dimension 600, produced:

```text
KKT nnz (upper):              2094
L nnz (strict lower):         2684
Sparse QDLDL factor time:   0.000068 s
Dense LDL factor time:      0.035030 s
Approx sparse factor data:  0.035 MiB
Dense L data:               2.747 MiB
```

The timing ratio is problem- and machine-dependent; the meaningful structural
change is sparse storage/work proportional to KKT and factor nonzeros rather
than dense quadratic storage.

## Source map

| Original SCS area | Fortran module |
|---|---|
| `scs.c`, `util.c` | `scs_solver` |
| `cones.c`, `exp_cone.c` | `scs_cones` |
| `normalize.c` | `scs_normalize` |
| `linalg.c`, sparse matrix helpers | `scs_linalg`, `scs_sparse` |
| `aa.c` | `scs_acceleration` |
| QDLDL | `scs_qdldl` |
| sparse direct KKT backend | `scs_ldlt` |
| R `scs.R` data/control interface | `scs_types`, `scs_solver` |
| R sparse conversion helpers | `scs_sparse` |

See `TRANSLATION_NOTES.md` for additional implementation notes.

## Licensing and provenance

The original R package declares **GPL-3**. This translated package is provided
under GPL-3.0-only; the full text is in `LICENSE`. The bundled upstream SCS
solver is MIT-licensed, and its license is retained in
`licenses/SCS-MIT-LICENSE.txt`.

The v0.2.0 `scs_qdldl.f90` module is a Fortran translation of bundled QDLDL,
which is Apache-2.0 licensed; its license is retained in
`licenses/QDLDL-APACHE-2.0.txt`. The bundled AMD license remains in
`licenses/AMD-BSD-3-CLAUSE.txt` for provenance, although AMD itself is not yet
translated or compiled. `ORIGINAL_DESCRIPTION` preserves the R package
authorship and package metadata.
