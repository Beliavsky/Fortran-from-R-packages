# Build and validation record

Validation date: 2026-09-02

## Toolchain available in the translation environment

- GNU Fortran: `GNU Fortran (Debian 14.2.0-19) 14.2.0`
- `fpm`: not installed (`command -v fpm` returned no path)
- `fprettify`: not installed (`command -v fprettify` returned no path)

An attempt to obtain the current FPM release binary through the available
sandbox download path was unsuccessful. The requested commands were then
attempted explicitly and each FPM command returned `fpm: command not found`;
`fprettify --version` likewise returned `fprettify: command not found`.
Consequently, this record does **not** claim successful execution of `fpm
build`, `fpm test`, `fpm run`, `fpm clean --all`, or `fprettify`. The manifest
was parsed successfully with Python's standard TOML parser.

The intended repository-root build remains:

```text
fpm build
fpm test
fpm run --example grain_demo
fpm clean --all
```

with sibling directories `../rfortran-core` and `../gRbase` present as declared
in `fpm.toml`.

## Direct GNU Fortran validation

The maintained gRain sources were compiled from clean temporary directories at
both `-O0` and `-O2` with:

```text
-std=f2018
-pedantic
-Wall
-Wextra
-Wimplicit-interface
-Werror=implicit-interface
-fcheck=all
-fbacktrace
```

The validation compiled the actual translated gRbase source modules used by
this package (`grbase_types`, `grbase_arrays`, `grbase_sets`, `grbase_tables`,
`grbase_graphs`, and `grbase_decompositions`). The local sandbox did not contain
the transitive sibling `rfortran-core` and `igraph` directories, so two
validation-only modules were supplied outside the package tree:

- `r_kinds`, defining the same `dp = real64` interface expected from
  rfortran-core; and
- `grbase_igraph`, supplying a small maximal-clique implementation with the
  interface expected by gRbase's decomposition module.

Neither validation module is included in the archive. No dependency source is
copied into `gRain/`.

This direct compile therefore validates the maintained gRain source syntax,
module interfaces used from gRbase, strict explicit-interface behavior, and the
package's deterministic numerical tests. It is not represented as a substitute
for the unavailable literal FPM dependency-resolution run.

## Deterministic tests

At both `-O0` and `-O2`, `test/test_grain.f90` completed successfully and
printed:

```text
All gRain deterministic tests passed.
```

The suite covers:

- CPT construction and normalization;
- compilation of a binary `1 -> 2 -> 3` network into a two-clique junction
  tree;
- collect/distribute propagation;
- independently specified prior marginals and a non-clique joint;
- hard evidence and its exact evidence probability;
- posterior marginal, joint, and conditional probabilities;
- soft evidence and its exact evidence probability;
- seeded posterior simulation and evidence-respecting samples;
- CPT replacement without retriangulation;
- logical AND/OR CPTs and Mendelian segregation probability;
- CPT estimation from categorical data; and
- clique marginal -> potential -> marginal round-trip consistency.

For the reference chain the independently calculated checks include
`P(node 3 = state 1) = 0.65` and
`P(node 1 | node 3 = state 1) = (0.6923076923076923,
0.3076923076923077)`.

## Example

At both optimization levels, `example/grain_demo.f90` completed successfully:

```text
P(node 3 = state 1) =   0.650000
P(node 1 | evidence) =   0.692308  0.307692
```

## Source audit

The maintained Fortran sources, test, and example were checked mechanically for
all requested source rules. The final audit covers 12 Fortran files, 45
procedures, and 153 dummy arguments and verifies:

- every dummy argument has explicit `INTENT` or `VALUE`;
- every dummy argument has its own declaration line;
- every dummy-argument declaration has a meaningful trailing `!!` FORD
  comment;
- no maintained Fortran line exceeds 132 columns;
- no code line contains a semicolon-separated statement;
- no `double precision`, `real*8`, `kind(0.0d0)`, or `d0`/`D` exponent literal
  occurs;
- no self-comparison NaN idiom is used;
- no duplicate maintained Fortran source file is present; and
- pure procedures are declared `pure` where their interfaces and dependency
  calls permit it. Network compilation remains impure because RIP construction
  crosses the igraph-backed gRbase bridge; two status-returning functions also
  cannot be pure under their current function interfaces.

## Package-content audit and cleanup

Before ZIP creation, the package tree was checked to contain no:

- object (`.o`) or module (`.mod`/`.smod`) files;
- executables or shared/static libraries;
- caches or Python bytecode;
- nested ZIP files; or
- copied dependency/BLAS/LAPACK/ARPACK source.

Because FPM itself is unavailable, `fpm clean --all` could not be executed.
Equivalent cleanup was performed by building only in external temporary
validation directories and confirming that the package tree itself contains no
build products before archiving.
