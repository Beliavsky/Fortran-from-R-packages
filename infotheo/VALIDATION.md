# Validation

## Environment

Validated 2026-10-12 on Linux using GNU Fortran 14.2.0, FPM 0.13.0,
R 4.5.0 with the exact CRAN infotheo 1.2.0.1 source package, and sibling
rfortran-core 0.2.0 from target-repository commit
`cbd437d6d92f5481706fd8ea5a7829aaff7b6e24`.
The Fortran package needs neither R nor Python for its normal tests.
Windows portability is a design requirement; native Windows execution has
not been performed.

## Executed tests

`fpm build`, `fpm test`, and `fpm run --example` passed. Checked-debug builds,
tests, and the example also passed using:

```sh
fpm build --profile debug --flag '-O0 -g -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow -Wall -Wextra -Wimplicit-interface -Werror=implicit-interface -std=f2018 -pedantic'
fpm test --profile debug --flag '-O0 -g -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow -Wall -Wextra -Wimplicit-interface -Werror=implicit-interface -std=f2018 -pedantic'
fpm run --example --profile debug --flag '-O0 -g -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow -Wall -Wextra -Wimplicit-interface -Werror=implicit-interface -std=f2018 -pedantic'
```

The release audit additionally runs all three commands with `--profile release`.
No non-IEEE optimization flags are used.

### Deterministic analytic and edge-case suite

`test/test_infotheo.f90` covers:

- All estimators for independent nonuniform counts [1,2,7]: empirical
  0.8018185525433373; MM 0.9018185525433373; corrected SG
  0.7730536807232604; legacy SG 0.7139682539682545; shrink
  0.8625257778161765 nats. The Dirichlet references were independently
  evaluated using SciPy digamma, not copied from the Fortran implementation.
- Three-way XOR: joint entropy 2 log(2), total correlation log(2),
  co-information -log(2), zero pairwise MI, and conditional MI log(2).
- Conditional entropy and information identities, multicolumn target and
  conditioning blocks, absent optional blocks, symmetric MI matrices,
  diagonal entropy, and the scalar/matrix negative-MM distinction.
- Exact/elemental nats conversion, constants, singleton samples, no complete
  observations, empty matrices, zero-column blocks, row-count mismatch,
  invalid method, the interaction dimension guard, and 100,000-observation
  shrinkage without upstream integer-square overflow.
- Separate marginal/joint complete-case deletion under missingness.
- All discretizers, legacy remainder cut ranks, ties, descending values,
  left- versus right-closed interval boundaries, constants, default bin count,
  NaNs and entirely missing columns/matrices, bad counts/methods, infinities,
  values beyond 32-bit extrema, the full finite binary64 range, and the maximum
  default-integer equal-width bin count without bin-label overflow.

### Actual-R reference tests

`tools/generate_reference.R` runs original R infotheo 1.2.0.1 and writes
`test/r_reference_data.f90`; `test/test_r_reference.f90` checks it.
It verifies **144 information values** across empirical, MM, legacy SG,
and shrinkage on both complete and partially missing data. Each case includes
joint and marginal entropies, block conditional entropy, scalar MI, scalar
conditional MI, total correlation, co-information, and every pairwise matrix
entry. Absolute tolerance is 5e-12, accommodating the shared digamma API's
small approximation difference from upstream.

The same fixture verifies **72 bin labels** across all three discretizers on
two columns, including ties and an unequal-remainder equal-frequency case.
To regenerate it:

```sh
Rscript tools/generate_reference.R > test/r_reference_data.f90
fpm test
```

The generator checks the R package version. It uses a data frame for the MI
matrix because the original R wrapper rejects an atomic matrix before its
matrix branch. `sg_r` is deliberately used for R's `sg`; comparing corrected
`sg` against the defective R estimator would be an invalid parity test.

### Independent semantic review

A separate reviewer checked **1,200 deterministic randomized cases** against
an independent contingency-count/SciPy-digamma/combinatorial oracle across
emp/mm/sg/sg_r/shrink and all information measures. Cases included up to four
X columns, up to two Y/S columns, and missing cells. All passed; maximum
absolute error was **8.633094239485217e-13**.

The reviewer also checked **231 discretization cases** against original R,
covering all three strategies, random values, bin counts, ties, and constants;
no mismatches were found on the safe parity domain. Deliberate repairs to
unsafe or incorrect upstream cases are verified analytically rather than
requiring agreement with defective results.

## Source/package audit

The maintained-source audit found zero errors and zero warnings: all five
maintained Fortran files are free form; all 45 dummy-argument declarations
have explicit INTENT and a separate meaningful trailing FORD `!!` comment;
coverage manifest counts/paths/public names agree; no duplicate Fortran
source, copied dependency source, disallowed semicolon-separated statement,
nonstandard real-kind declaration, self-comparison NaN test, or forbidden
manifest system-library/fast-math flag was found.

The original R/C++ files in `upstream/` were compared byte-for-byte with the
official archive. Build artifacts are removed before packaging. The release
artifact contains only the `infotheo/` directory and relies on its sibling
rfortran-core; no dependency source is vendored.
