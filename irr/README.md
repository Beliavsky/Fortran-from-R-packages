# irr: agreement and reliability in modern Fortran

A modern free-form Fortran translation of the computational code of **irr 0.85**,
by Matthias Gamer, Jim Lemon, Ian Fellows and Puspendra Singh. GPL version 2 or
later; see COPYING and NOTICE.md. All maintained arithmetic uses `dp` from
`r_kinds`. Numerical procedures are pure except the LP-backed sample-size call,
whose shared solver is impure.

## Build and run

Place `irr` at the repository root beside `rfortran-core`, `rfortran-linalg`, and
`lpSolve`. No dependencies are copied into this package.

```sh
fpm build
fpm test
fpm run --example agreement_demo
```

A current gfortran/FPM installation is sufficient, including on Windows. There
are no system BLAS/LAPACK links. The sibling `rfortran-linalg` provides the pure
`solve_system` interface and brings its tested, pinned Fortran LAPACK dependency.
The first build compiles that large dependency and takes longer. `lpSolve` is
referenced under its actual FPM package name `lpsolve-fortran`.

Do not compile with `-ffast-math`, `-Ofast`, `-ffinite-math-only`, or equivalent
non-IEEE assumptions: NaNs represent missing data and unavailable statistics.

## Quick start

```fortran
use irr, only: dp, irr_result, irr_ok, kappa2_table
real(dp) :: counts(2,2)
type(irr_result) :: result
counts = reshape([20.0_dp,10.0_dp,5.0_dp,15.0_dp], [2,2])
result = kappa2_table(counts)
if (result%status /= irr_ok) error stop 'Kappa failed'
print *, result%value ! 0.4
```

The public facade `irr` re-exports every mapped procedure and `dp`. See the
trailing FORD comments on each procedure argument and `example/agreement_demo.f90`.

## Data and result contract

- Most procedures accept a real(dp) **subjects × raters** matrix. `kripp_alpha`
  intentionally follows R's opposite orientation: **raters × subjects**. `iota`
  takes **subjects × raters × variables**.
- Categorical labels are finite numbers, compared exactly and sorted numerically.
  Encode strings or factors as ordered numeric labels before calling. R factor
  level metadata, unused levels, character ordering, partial string matching,
  matrix/data-frame coercion, printing and S3 objects are not reproduced.
- Most functions omit any subject containing a NaN (listwise). `iota` applies
  this across every variable, an explicit policy replacing ambiguous R list
  handling. `kripp_alpha` uses all available rating pairs within each subject;
  subjects with fewer than two available ratings contribute no coincidences.
  Infinities are invalid observations. `rel_inter_intra` requires a complete,
  balanced finite matrix and never silently drops cells.
- `irr_result` provides `value`, `statistic`, `p_value`, `df1`, `df2`, `lbound`,
  `ubound`, `subjects`, `raters`, `status`, and `dropped`. Inapplicable numeric
  fields are quiet NaNs. `levels` provides category order where relevant;
  Fleiss `detail(:,1:3)` is kappa/z/p without R's presentation-only rounding;
  alpha `coincidence` is the complete coincidence matrix.
- Status constants: `irr_ok=0`, `irr_invalid=1`, `irr_degenerate=2`,
  `irr_resource=3`. A degenerate inferential calculation can preserve a valid
  coefficient (notably Light's kappa). Check each needed output, not just NaN
  propagation. Zero standard errors may produce valid infinite statistics.
- Table APIs require explicit square nonnegative counts; fractional counts are
  accepted. Subject count is the nearest integer table total. Score APIs never
  guess that a square score matrix is a contingency table. Singular covariance
  systems return `irr_degenerate` rather than aborting.
- General bounds: 10,000,000 input cells, 1,000 raters and 256 distinct categories.
  `iota`'s flattened rater-by-variable width is also limited to 1,000. These are
  safety ceilings, not performance guarantees: Light's kappa and alpha include
  pairwise-rater work. N2 retains R's stricter 2–10-category domain.

## Numerical compatibility and deliberate choices

1. `kripp_alpha(..., legacy_complete_weights=.true.)` defaults to R 0.85's
   behavior: if no value is missing anywhere, each within-subject pair has
   weight one. With three or more raters this differs from normalized
   Krippendorff coincidences. Pass `.false.` to always divide pairs by
   the available-rater count minus one. The corrected mode is invariant to
   appending an entirely missing subject; the R-compatible mode need not be.
2. `stuart_maxwell_mh(..., legacy_equal_marginals=.true.)` defaults to R's
   deletion of any category whose two marginal totals agree. Equal marginals
   need not mean perfect agreement. Pass `.false.` to retain discordant
   categories and remove only isolated perfect-agreement categories. The
   nonlegacy mode uses the usual reduced nonsingular covariance calculation;
   disconnected/singular covariance structures still return degenerate status.
3. Cohen numeric `costs` are normalized and assembled exactly as R, including
   its unusual asymmetric matrix for unsorted custom costs. Ordinary equal and
   squared weights remain symmetric. Constant cost vectors are rejected.
4. Light's upstream inferential variance can be negative even with a valid mean
   kappa. Its coefficient is returned with NaN inference and degenerate status.
   Fisher averaging drops exact perfect correlations; their number is `dropped`.
   Inadequate sample size or all-perfect pairs yield unavailable inference.
5. ICC preserves the upstream McGraw–Wong confidence-bound formulas, including
   its average-agreement degrees-of-freedom expression. Perfect nonconstant
   agreement receives the limiting [1,1] interval; all-constant data is
   degenerate. `kind` is the Fortran argument corresponding to R's `type`.
6. `n_cohen_kappa` and `n2_cohen_kappa` return the upward-rounded sample size as
   real(dp), with optional `raw_n`, avoiding default-integer overflow. Invalid
   or infeasible probability designs are rejected rather than using negative
   cell probabilities. N2 accepts a marginal sum within 64 machine epsilons of
   one instead of R's exact equality check. Zero marginal categories are allowed.
7. `rel_inter_intra(x, nrater, result, ...)` requires at least two subjects,
   two raters and two repeats per rater. Columns are rater-major, with contiguous
   repeats. The result's `rohat`, `fs`, `pvalue`, `lowvalue` and `sem` vectors
   preserve upstream ordering; see API_COVERAGE.md for the exact layout. Null
   reliabilities must be in [0,1), and the confidence level in (0,1).
8. Upper-tail probabilities are computed directly through shared distribution
   routines, avoiding the cancellation in R code written as `1-pf` or `1-pnorm`.
   Warnings, labels and other presentation data are replaced by numeric status.

## Validation

Deterministic `fpm test` tests compare all 19 computational functions with actual
R 4.5.0 / irr 0.85 / lpSolve 5.6.23 values. Reproducible oracle scripts are in
`validation/`; their generated constants and repeated-ANOVA references ship
with the tests, so R is not required to build or test. Independent checks cover
binary contingency formulas, normalized alpha coincidence counts, ANOVA sums
of squares, transformation invariance, asymmetric scores, weights, missingness,
constant/perfect/zero agreement, infeasible LPs and resource limits. The LP
translation also matched 150 deterministic randomized 2–10-category R designs.

## Translation coverage

Package status: **substantial**. **19 of 19 (100%)** exported computational R
functions are mapped. This fraction measures mapped computational functions,
not complete R compatibility. Every mapping is substantial because numeric
interfaces replace R objects/coercions and explicit validation differs at edges.
The inventory was obtained from the complete 0.85 `NAMESPACE` and `R/` tree;
`print.irrlist` and `print.icclist` are excluded presentation methods. Internal
nested helpers and imported dependency functions are excluded.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
|---|---|---|---|
| `agree` | `irr_agreement` | `agree` | substantial |
| `bhapkar` | `irr_marginal` | `bhapkar`, `bhapkar_table` | substantial |
| `finn` | `irr_agreement` | `finn` | substantial |
| `icc` | `irr_icc` | `icc` | substantial |
| `iota` | `irr_agreement` | `iota` | substantial |
| `kappa2` | `irr_kappa` | `kappa2`, `kappa2_table` | substantial |
| `kappam.fleiss` | `irr_kappa` | `kappam_fleiss` | substantial |
| `kappam.light` | `irr_kappa` | `kappam_light` | substantial |
| `kendall` | `irr_agreement` | `kendall` | substantial |
| `kripp.alpha` | `irr_agreement` | `kripp_alpha` | substantial |
| `maxwell` | `irr_agreement` | `maxwell` | substantial |
| `meancor` | `irr_agreement` | `meancor` | substantial |
| `meanrho` | `irr_agreement` | `meanrho` | substantial |
| `N.cohen.kappa` | `irr_designs` | `n_cohen_kappa` | substantial |
| `N2.cohen.kappa` | `irr_designs` | `n2_cohen_kappa` | substantial |
| `rater.bias` | `irr_marginal` | `rater_bias`, `rater_bias_ratings` | substantial |
| `relInterIntra` | `irr_designs` | `rel_inter_intra` | substantial |
| `robinson` | `irr_agreement` | `robinson` | substantial |
| `stuart.maxwell.mh` | `irr_marginal` | `stuart_maxwell_mh`, `stuart_maxwell_ratings` | substantial |
