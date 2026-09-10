# Validation record

## Current 122/122 checkpoint

This checkpoint maps **122 of 122 (100.0%)** functions in the documented
computational coverage basis. The package status remains **substantial** because
several mappings are intentionally partial and typed Fortran interfaces do not
reproduce every R/S3 option, cache, class, warning, RNG stream, or FFTW-backed
acceleration path.

The current tree contains 16 library modules and 27 maintained free-form Fortran
source/test/demo/example files (4,573 lines). The retained upstream source and
coverage/provenance documentation are also included.


## CSSA sibling-API compatibility fix

The current release fixes a compile-time mismatch with the sibling `svd`
package. `ztrlan_svd_result` intentionally returns singular values and left
singular vectors only (`d` and `u`); it has no `v` component. Complex SSA now
constructs the required right singular vectors from the conventional complex
SVD identity `v_j = A^H u_j / sigma_j`, with zero vectors used for numerically
zero singular values.

`test_decomposition` now includes a complex rank-one geometric series and checks
that a one-component CSSA decomposition reconstructs it to `1e-10`. In addition,
the affected modules were compiled and this regression was executed against a
validation `svd` stand-in whose `ztrlan_svd_result` exactly matches the current
sibling API (no `v` member). The regression reports `test_cssa: PASS`.

## Final IOSSA/EOSSA parity round

This checkpoint closes the four previously unmapped computational names:

- `iossa` -> `iossa_ssa`
- `iossa.ssa` -> `iossa_ssa`
- `eossa` -> `eossa_ssa`
- `eossa.ssa` -> `eossa_ssa`

`iossa_ssa` implements the real 1-D iterative O-SSA core. Selected eigentriples
are represented by one-based indices and integer group labels. Each iteration
reconstructs the current nested groups, recomputes their leading ordinary-SSA
triples, applies the upstream kappa separation/balance rule, performs the
low-rank oblique SVD transform, Hankelizes the resulting grouped terms, and
checks the upstream RMS-style convergence criterion. The result reports optional
iteration and convergence outputs in place of R's cached `iossa.result` list.

`eossa_ssa` implements the default real 1-D column-subspace, least-squares EOSSA
path. It forms the one-step ESPRIT shift operator, uses the current
`rfortran-linalg` general complex eigensolver, sorts roots by absolute argument,
clusters `(Re(root), abs(Im(root)))` points by deterministic complete linkage,
realifies each clustered complex eigenspace by SVD, and maps the rotated basis
back to an oblique decomposition. Multidimensional EOSSA, row-subspace mode,
TLS shift solves, and R `hclust` label identity are deliberate omissions.

A new deterministic `test_iterative_oblique` uses the analytic rank-2 series

```text
x_t = 2 * 0.82^t + 0.7 * (-0.45)^t,  t = 0,...,11.
```

For this case both I-OSSA and EOSSA recover the two exact rank-1 exponential
components. The test checks each reconstructed component, total trajectory
preservation, I-OSSA convergence, and frequency-ordered EOSSA clusters.

## Current checked and optimized validation

GNU Fortran 14.2.0 compiled the exact current source tree against external,
API-compatible numerical stand-ins for the sibling `svd` and
`rfortran-linalg` modules. The stand-ins live outside `Rssa/`, are used only to
exercise this package in an environment without FPM/network access, and are not
included in the release archive. The public interfaces used by the new EOSSA
code were also checked against the current repository `rfortran-linalg` source,
including `general_complex_eigen` and matrix `least_squares_svd`.

Checked flags:

```text
-std=f2018 -O0 -g -Wall -Wextra -Wconversion-extra -Wimplicit-interface
-Werror -Wno-maybe-uninitialized -fcheck=all -fbacktrace
-ffree-line-length-132
```

Optimized flags:

```text
-std=f2018 -O2 -Wall -Wextra -Wconversion-extra -Wimplicit-interface
-Werror -Wno-maybe-uninitialized -ffree-line-length-132
```

No fast-math, finite-math-only, or equivalent assumptions are used.
`-Wno-maybe-uninitialized` is retained only for GNU Fortran's conservative
warnings around several pre-existing allocatable-component assignments.

All eight deterministic tests pass in both checked and optimized builds:

```text
test_decomposition: PASS
test_forecast: PASS
test_fossa: PASS
test_gap_group: PASS
test_iterative_oblique: PASS
test_matrices: PASS
test_oblique: PASS
test_projection_forecast: PASS
```

The demo and both examples also compile and run in both configurations:

```text
rssa_demo: PASS
basic_ssa: PASS
forecast_gapfill: PASS
```

## Earlier bounded parity rounds retained in this tree

The current source also includes the previously validated bounded additions:

- projection-SSA recurrent/vector forecast corrections;
- diagonal weighted-oblique decomposition and ordinary/oblique correlations;
- real 1-D FOSSA with finite-gamma and filter-only modes;
- the upstream `decompose.ossa` refusal to continue an oblique decomposition.

These paths remain covered by `test_projection_forecast`, `test_oblique`, and
`test_fossa` respectively.

## Static release audit

`python scripts/audit.py` checks manifest/README coverage agreement,
ASCII/free-form source, the 132-column limit, duplicate maintained source files,
semicolon-separated executable statements, disallowed real-kind forms,
D-exponent literals, self-comparison NaN tests, and packaged build/archive
artifacts.

For this checkpoint it reports:

```text
audit: PASS (27 Fortran files, 122 mapping entries)
```

The manifest parses successfully as TOML and records:

```text
r_functions_total = 122
r_functions_translated = 122
frac_functions_translated = 1.0
untranslated_r_functions = []
```

## Recovery history

A runtime reset during the earlier 110/122 work removed the exact pre-reset
working tree. The 118/122 checkpoint was reconstructed from surviving source and
patch history, then recompiled. The current 122/122 tree has subsequently been
compiled and its complete current test suite executed as described above; the
old pre-reset results are no longer the only numerical evidence for the archive.

## FPM and formatting tools

The required FPM commands are attempted from the package root before release.
In this environment `fpm` is not installed, so they return command-not-found and
a literal FPM/Windows build cannot be claimed here. `fpm.toml` uses only sibling
path dependencies:

```toml
svd = { path = "../svd" }
rfortran-linalg = { path = "../rfortran-linalg" }
```

`fprettify` is also not installed in the current environment. The maintained
source is nevertheless audited as free-form, ASCII, and at most 132 columns.

## Dependency and archive policy

No dependency source is copied into this package. The release contains no BLAS,
LAPACK, ARPACK, FFTW, `r.f90`, `r_mod.f90`, translated dependency source,
compiler objects, module files, executables, caches, or nested ZIP archives.

The original supplied archive SHA-256 is recorded in
`provenance/source-archive.sha256`, and per-file SHA-256 hashes of the retained
upstream tree are recorded in `provenance/upstream-files.sha256`.
