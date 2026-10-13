# Validation

## Deterministic executable tests

Normal `fpm test` requires no R installation and runs:

1. `test_oracle`: 132 actual-R reference samples; statistic, alternative zero
   convention and table p-value for every case; full predecessor/successor,
   modal interval and modal/full GCM/LCM arrays for 65 cases.
2. 987 R-generated table-calibration cases, including tied knots, interpolation
   between sample sizes, tails and the exact n=72000 boundary.
3. A 257-replicate Monte Carlo oracle using saved actual-R uniform samples and
   exact inclusive exceedance counting.
4. `test_properties`: the analytic two-point-mixture formula over all split
   counts for n=2..40 and uniform-grid closed forms under both minimum conventions.
5. Independent geometric checks that full hull segments have monotone slopes
   and lie below/above every empirical-CDF index point as appropriate.
6. Permutation, reflection and affine invariance; huge and tiny power-of-two
   scales; subnormal spacing; large offsets; overflowing signed coordinate
   ranges; tiny clusters adjacent to a huge outlier and its reflection.
7. Empty, singleton, constant, tied, nonfinite and missing-data contracts,
   sorted-input validation, direct-calibration errors, malformed null samples,
   invalid/default B, and repeatability after restoring the intrinsic RNG seed.
8. Full `dip_test` calls on 72000- and 72001-point grids verify the finite-table
   boundary and asymptotic flag without an out-of-bounds access.

The ordinary-data absolute tolerance is 3e-14 for dips, 3e-12 for p-values;
reference node arrays are checked exactly. Table probability comparisons allow
8e-12 around steep interpolation segments. Analytic extreme-scale tests use
3e-15 to 5e-15. Tests never accept NaN as a successful numerical comparison.

Both examples are executable and check status before showing results.

## Additional independent review

An independent reviewer compared 5,402 cases directly with the current upstream
C kernel. All statistic comparisons passed, including the adversarial vector
`[0:8]*1e-303, [1000:1009]*1e-303, 1e300`, whose dip is 0.2232.
This test exposed a preliminary global-rescaling approach that erased the tiny
clusters. The final implementation instead uses **local binary fraction/exponent
pairs**, and the regression is now permanently included in `fpm test`.

A further 6,300 extreme cases were compared against an independently
widened-long-double version of the corrected C kernel: scales down to 2^-1070 and
up to 2^1023, mixed 1e-300..1e300 ranges, and both minimum conventions. No dip
error exceeded 2e-13. Alternate hull choices occurred in roundoff-sensitive ties
(9 ordinary differential cases and 6 extreme cases), consistent with the
documented single-kind/long-double distinction. These external review harnesses
are not runtime or build dependencies.

## Build profiles

Default, checked-debug and optimized-release FPM build/test/all-example stages
were run successfully with **gfortran 14.2.0** and **FPM 0.13.0**. Checked debug
uses `-std=f2018 -Wall -Wextra -Wimplicit-interface -fcheck=all -fbacktrace
-ffpe-trap=invalid,zero,overflow`. Release uses FPM's normal release profile.
No fast-math flags are used. The release coordinator independently reruns the
nine stages against a frozen source snapshot.

Native Windows execution was not available. The source and FPM manifest use
portable free-form Fortran, a sibling path dependency and no system-library
links; native Windows runtime validation remains unperformed.

## Regeneration (optional, requires R)

Install the exact upstream diptest 0.77-2 package, then from this package root:

```sh
Rscript tools/generate_table.R
Rscript tools/generate_references.R
fpm test
```

The R script may report warnings that duplicate quantile abscissae are being
collapsed. These are expected and validate the upstream `approx(ties=mean)`
semantics rather than indicating failed data generation.
