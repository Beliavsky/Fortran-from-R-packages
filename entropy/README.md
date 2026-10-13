# entropy

Modern free-form Fortran/FPM translation of **R entropy 1.3.2** by Jean Hausser
and Korbinian Strimmer, licensed **GPL-3.0-or-later**. Includes entropy,
mutual information, KL/Pearson divergences, shrinkage and Dirichlet estimators,
G/Pearson test statistics and equal-width discretizers.

## Build and use

Place this directory at the root of
[Fortran-from-R-packages](https://github.com/Beliavsky/Fortran-from-R-packages),
next to its existing `rfortran-core` directory:

```text
fpm build
fpm test
fpm run --example estimate_information
```

Only the sibling FPM dependency `../rfortran-core` is needed. Its `dp`,
stable elementary transforms and chi-squared distribution are reused directly.
No dependency sources are copied, and there are no BLAS/LAPACK or other system
library links. Maintained sources use Fortran 2008 features compatible with
modern gfortran on Windows and other platforms; execution was verified on Linux
with GNU Fortran 14.2.0 and FPM 0.13.0, not on a Windows host.

```fortran
use entropy_module, only : dp, entropy, freqs_shrink, mi_empirical
real(dp) :: counts(4), joint(2, 2)
counts = [8.0_dp, 3.0_dp, 1.0_dp, 0.0_dp]
print *, entropy(counts, method='shrink', unit='log2')
print *, freqs_shrink(counts)
joint = reshape([4.0_dp, 0.0_dp, 0.0_dp, 4.0_dp], [2, 2])
print *, mi_empirical(joint, 'log2')
```

All computational procedures are pure; scalar helpers are elemental where
applicable. The public module re-exports the single real kind `dp` from `r_kinds`.
All real arrays and scalar real arguments must use that kind. Details and
calling conventions are in [API_COVERAGE.md](API_COVERAGE.md).

## Numerical and interface contract

- Units are the exact strings `log` (default), `log2`, and `log10`.
- Methods are exact, case-sensitive R names. `entropy` defaults to `ML`;
  available methods are `ML`, `MM`, `Jeffreys`, `Laplace`, `SG`, `minimax`, `CS`,
  and `shrink`. `NSB` is deliberately not approximated and returns NaN.
- Inputs must be nonempty finite nonnegative counts/weights. Positive totals
  are required unless positive Dirichlet priors supply the mass. Invalid data,
  unsupported units/methods, shape mismatches and scalar numerical errors
  produce quiet NaNs; zero-length probability arrays stay zero length.
- The probability normalization and Dirichlet prior addition are scaled to
  avoid needless overflow. Count-based corrections require a representable
  sample total and otherwise yield NaN. Chao–Shen uses stable log1p/expm1.
- Unlike permissive R arithmetic, negative/nonfinite data, negative priors,
  and supplied shrinkage intensities outside [0,1] are rejected. The empty or
  zero-count Chao–Shen estimate is NaN, rather than upstream's incidental zero.
- `freqs_empirical`, `freqs_dirichlet` and `freqs_shrink` support rank-one and
  rank-two arrays. Entropy and two-distribution divergence routines accept
  rank-one arrays; flatten matrices with `reshape`. Mutual-information and
  independence routines take rank-two arrays. No implicit R recycling occurs.
- `freqs_dirichlet` accepts a scalar prior or a same-shape prior array. The
  other Dirichlet wrappers accept scalar priors; compose the frequency and
  plugin operations for cell-specific priors.
- Zero numerator probability contributes zero to KL, including a common zero.
  Positive probability over zero reference probability yields positive
  infinity. Pearson divergence retains upstream NaN for any common-zero cell;
  this includes independence tables with a zero row or column margin.
- Shrinkage verbosity and R attributes are omitted. `shrinkage_intensity(y)`
  exposes the automatic intensity separately (flatten tables first). An all-zero
  sample has automatic intensity one but NaN estimated probabilities, as in R.
  `freqs` returns NaN arrays for `MM`, `CS`, and `NSB`, as upstream does.
- Test functions return `type(test_result)` with `stat`, integer `df`, and
  `pval`. P-values use a direct upper tail, improving cancellation-prone
  upstream `1-pchisq`. Log-unit conversion follows R even for Pearson statistics.
  Zero-degree-of-freedom and infinite statistics are handled explicitly. Negative
  statistics caused by floating-point cancellation retain their value but have
  upper-tail p-value one, matching R.
- Discretizers return allocatable integer counts and optional status (0 success,
  1 invalid input). Bins are right-closed and include the lowest endpoint.
  Explicit ranges omit nonfinite/out-of-range observations; default ranges
  reject nonfinite or empty input and constant ranges. Computed breaks must
  remain distinct at dp precision, as required by R cut. R interval labels and
  table objects are omitted. Equal-width breaks follow R's step-times-index arithmetic for finite widths
  and use convex interpolation when the endpoint difference would overflow.

Do not compile with `-ffast-math`, `-Ofast`, `-ffinite-math-only`, or equivalent
assumptions: NaN and infinity semantics are part of the API.

## Tests and validation

`test/test_r_reference.f90` contains **826 comparisons** evaluated from the
unmodified pinned R source using GNU R 4.5.0. All 33 mapped exports are exercised,
including four sample scenarios, all three units, fixed/automatic shrinkage,
Dirichlet priors, methods, statistics/p-values and both discretizers.
`test/test_properties.f90` adds **139 independent identity and edge checks**,
including missing data, zero cells, invalid dimensions, boundaries, very large
and tiny weights, unit invariance and exact mathematical identities.

To regenerate the R-based golden test (R is not required to build or run tests):

```text
Rscript --vanilla validation/generate_reference.R
```

Optional formatting afterward: `fprettify -i 3 -w 3 test/test_r_reference.f90`.
Debug and release build/test/example execution are required for releases.

## Translation coverage

Status: **substantial**. **33 of 34 (97.06%)** exported computational R functions
have mappings. The fraction measures mapped computational functions, not complete
R interface compatibility. The 34-function denominator conservatively retains
`entropy.NSB`, even though its R code only delegates to an external executable;
the NSB algorithm is not included in the upstream package. It is untranslated.
No other exports are excluded, and no internal helpers inflate the count.
Per-function statuses describe numerical scope, with interface limits above and
in `fpm.toml` applying even to complete numerical mappings.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
| --- | --- | --- | --- |
| `chi2.Dirichlet` | `entropy_module` | `chi2_dirichlet` | substantial |
| `chi2.empirical` | `entropy_module` | `chi2_empirical` | complete |
| `chi2.plugin` | `entropy_module` | `chi2_plugin` | complete |
| `chi2.shrink` | `entropy_module` | `chi2_shrink` | substantial |
| `chi2indep.Dirichlet` | `entropy_module` | `chi2indep_dirichlet` | substantial |
| `chi2indep.empirical` | `entropy_module` | `chi2indep_empirical` | complete |
| `chi2indep.plugin` | `entropy_module` | `chi2indep_plugin` | complete |
| `chi2indep.shrink` | `entropy_module` | `chi2indep_shrink` | substantial |
| `Gstat` | `entropy_module` | `gstat` | complete |
| `Gstatindep` | `entropy_module` | `gstatindep` | complete |
| `chi2stat` | `entropy_module` | `chi2stat` | complete |
| `chi2statindep` | `entropy_module` | `chi2statindep` | complete |
| `discretize` | `entropy_module` | `discretize` | substantial |
| `discretize2d` | `entropy_module` | `discretize2d` | substantial |
| `entropy` | `entropy_module` | `entropy` | substantial |
| `entropy.ChaoShen` | `entropy_module` | `entropy_chaoshen` | complete |
| `entropy.Dirichlet` | `entropy_module` | `entropy_dirichlet` | substantial |
| `entropy.empirical` | `entropy_module` | `entropy_empirical` | complete |
| `entropy.MillerMadow` | `entropy_module` | `entropy_millermadow` | complete |
| `entropy.plugin` | `entropy_module` | `entropy_plugin` | complete |
| `entropy.shrink` | `entropy_module` | `entropy_shrink` | substantial |
| `freqs` | `entropy_module` | `freqs` | substantial |
| `freqs.Dirichlet` | `entropy_module` | `freqs_dirichlet` | substantial |
| `freqs.empirical` | `entropy_module` | `freqs_empirical` | complete |
| `freqs.shrink` | `entropy_module` | `freqs_shrink` | substantial |
| `mi.Dirichlet` | `entropy_module` | `mi_dirichlet` | substantial |
| `mi.empirical` | `entropy_module` | `mi_empirical` | complete |
| `mi.plugin` | `entropy_module` | `mi_plugin` | complete |
| `mi.shrink` | `entropy_module` | `mi_shrink` | substantial |
| `KL.Dirichlet` | `entropy_module` | `kl_dirichlet` | substantial |
| `KL.empirical` | `entropy_module` | `kl_empirical` | complete |
| `KL.plugin` | `entropy_module` | `kl_plugin` | complete |
| `KL.shrink` | `entropy_module` | `kl_shrink` | substantial |
| `entropy.NSB` | — | — | untranslated |

See [NOTICE.md](NOTICE.md) for exact upstream provenance and copyrights,
[LICENSE](LICENSE) for the complete GPL v3 text, and [API_COVERAGE.md](API_COVERAGE.md) for call shapes.
