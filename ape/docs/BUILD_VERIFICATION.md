# Build verification

## Compiler verification performed for this translation

The maintained Fortran sources, deterministic tests, and all twenty-two examples were
compiled and run from clean external build directories with GNU Fortran 14.2.0
using both of these flag sets:

```text
-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -Wsurprising -pedantic -O2
```

and

```text
-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -Wsurprising -pedantic -O0 -fcheck=all -fbacktrace
```

Both builds passed without `-ffree-line-length-none`; every maintained `.f90`
line was also checked to be at most 132 columns so normal FPM/gfortran free-form
compilation does not depend on a line-length extension. `fprettify` is not
installed in this packaging environment, so it could not be run here.
`test/test_ape.f90` printed:

```text
All ape deterministic tests passed.
```

All 22 files in `example/` compiled and ran successfully with both flag sets, including the new discrete ACE, PGLS, chronopl/chronos, comparative-model, binary-PGLMM, and reconstruct examples.

The direct compiler verification used an external sibling `rfortran-core`
`r_kinds` module. For the new `rfortran-linalg` calls, validation used an
API-compatible external scaffold backed by the container's LAPACK solely as a
test harness; it is outside the package and is not referenced by `fpm.toml`.
The shipped manifest instead points to the repository sibling `rfortran-linalg`,
whose FPM manifest pins `fortran-lapack`. All temporary dependencies and compiler
products were kept outside the `ape` directory and are not part of the ZIP.

During development, the newly translated NJ*/BIONJ*/MVR* and triangle-method
routines were additionally compared with small standalone builds of the
corresponding upstream ape C kernels. Those comparison harnesses are validation
artifacts only and are not distributed with this package. The deterministic
Fortran tests retain representative reference values from those comparisons.

FastME received a broader native-kernel regression during this parity pass. Five
modes--OLS insertion, OLS+NNI, BME insertion, BME+bNNI, and BME+bNNI+SPR--were
compared with standalone builds of ape's `me_o`/`me_b` C implementation on 750
generated complete distance matrices with 6, 7, or 8 taxa. Every canonical
split topology matched; the largest reconstructed tip-distance discrepancy was
about `1.4e-14`.


The new analytical workflows have deterministic numerical regressions as well.
PCoA checks uncorrected, Lingoes, and Cailliez eigen/correction values; continuous
ACE checks closed-form Brownian ML/REML states, variances, and likelihoods; the
skyline tests check modern/old-style population-size likelihoods and AICc; and
standard plus extended (`bd.ext`) birth-death tests check fixed deviances,
interior optima, and zero-extinction boundary behavior.

DNAbin ambiguity behavior was likewise checked directly against extracted
upstream kernels used only for development. `seg.sites` strict/non-strict
behavior matched on all 4,913 possible three-state columns over ape's 17 DNAbin
states. DNA translation matched all 29,478 combinations of 17^3 codons across
genetic codes 1--6. These comparison harnesses and upstream C sources are not
distributed in this package.

The current large-parity pass adds deterministic regressions for discrete ACE,
phylogenetic GLS/PGLS, `chronopl`, `chronos` clock/correlated/relaxed/discrete
objectives, `compar.ou`, `compar.lynch`, `corphylo`, binary PGLMM PQL/REML, and
the `reconstruct` family. A fresh clean compile of the current source and all 22
examples passed under both compiler flag sets immediately before packaging.

## Source-policy audit

The final maintained source set was checked for the requested translation
policy. The audit found no violations among all 49 library/test/example `.f90` files:

- all translated/new source uses free-form `.f90` files;
- no copied `r.f90`, `r_mod.f90`, `r_kinds.f90`, BLAS, LAPACK, ARPACK, or
  translated dependency source is present;
- no duplicate Fortran source files are present;
- no semicolon-separated Fortran statements are present;
- no self-comparison NaN idioms are present;
- no legacy `double precision`, `real*8`, `kind(0.0d0)`, or D-exponent real
  literals are present;
- every dummy argument has an explicit `INTENT` or `VALUE`, is declared alone,
  and has a trailing meaningful FORD `!!` comment;
- no build products, module files, executables, caches, or nested ZIP files are
  present;
- `fpm.toml` contains no fast-math assumptions or system BLAS/LAPACK links.

The preserved upstream `COPYING`, `DESCRIPTION`, and `inst/CITATION` files were
also byte-compared with `COPYING`, `UPSTREAM_DESCRIPTION`, and
`UPSTREAM_CITATION` in this package and were identical.

## FPM execution limitation in the packaging environment

The packaging execution environment does not contain an `fpm` executable. An
attempt to fetch the official standalone FPM release was also blocked by the
container's outbound DNS/network isolation. Fresh attempts to invoke all
requested commands returned shell status 127 with `fpm: command not found`:

```text
fpm --version
fpm build
fpm test
fpm clean --all
```

Therefore this archive does **not** claim that FPM itself was executed here.
The FPM manifest is a standard sibling-path project and the same source was
validated directly with gfortran as described above. Before merging into
`Fortran-from-R-packages`, run from the top-level `ape` directory with the
repository's sibling `rfortran-core` present:

```text
fpm build
fpm test
fpm run --example nj_example
fpm run --example incomplete_reconstruction_example
fpm run --example fastme_example
fpm run --example dna_distance_example
fpm run --example tree_statistics_example
fpm run --example consensus_example
fpm run --example diversification_example
fpm run --example chrono_mpl_example
fpm run --example pcoa_example
fpm run --example ace_likelihood_example
fpm run --example skyline_example
fpm run --example birthdeath_example
fpm run --example birthdeath_extended_example
fpm run --example discrete_ace_example
fpm run --example pgls_example
fpm run --example chronopl_example
fpm run --example chronos_clock_example
fpm run --example compar_ou_example
fpm run --example compar_lynch_example
fpm run --example corphylo_example
fpm run --example binary_pglmm_example
fpm run --example reconstruct_example
fpm clean --all
```

No system BLAS/LAPACK link flags are present in this package.
