# Reference fixtures and data attribution

These text fixtures were generated on 2026-10-12 by running the **actual R
package diptest 0.77-2**, including its corrected compiled C kernel, under
**R 4.5.0 (2025-04-11)**. The generation script is
`tools/generate_references.R`; it checks the installed package version.

- `oracle_cases.txt`: 132 sample cases, each with observed dip, `min.is.0` dip and
  the tabulated p-value. Sixty-five cases also contain all original `mn`, `mj`,
  modal and full-sample hull node arrays and modal endpoints. Diagnostic checks
  deliberately use cases with stable hull choices rather than claiming universal
  cross-platform tie equivalence.
- `pvalue_cases.txt`: 987 independent R `approx` results covering small n,
  original table nodes, between-row values, repeated quantile knots, tails,
  n=71999/72000/72001 and n=100000. Small-n boundary behavior comes from the
  upstream wrapper. The interpolation oracle follows the original R formula
  directly and does not invoke the Fortran implementation.
- `monte_carlo.txt`: the actual R-generated 10-by-257 uniform draws,
  observed data and R exceedance proportion. R seed: 20261012. The exact sample
  matrix, rather than a claim of equal RNG algorithms, makes the test portable.

The other generated cases use R seed 217. Reference values and inputs are
recorded with 17 significant decimal digits.

## Dataset sources

- **statfaculty** (63 observations): distributed in diptest; its exact source is
  retained at `upstream/data/statfaculty.R`. This is the historical Hartigan
  statistics-faculty age example and the deprecated `exHartigan` data alias.
  Licensed with diptest under GPL (>=2).
- **faithful$eruptions** (272 observations): from R's datasets package, Old
  Faithful eruption durations in minutes. Source credited by R: **W. Härdle**;
  see Härdle (1991), *Smoothing Techniques with Implementation in S*, Springer,
  and Azzalini and Bowman (1990), *Applied Statistics* 39, 357–365,
  DOI 10.2307/2347385.
- **precip** (70 observations): from R's datasets package, annual precipitation
  in inches. Source credited by R: *Statistical Abstracts of the United States*,
  1975; reference McNeil (1977), *Interactive Data Analysis*, Wiley.

R's datasets DESCRIPTION identifies **R Core Team and contributors worldwide**
and says `License: Part of R 4.5.0`. R's overall distribution notice permits
**GPL version 2 or version 3**. Those terms are preserved for the mixed
`oracle_cases.txt` fixture containing these datasets. The original metadata is
in `R-datasets-DESCRIPTION`, the unchanged distribution-wide notice in
`../../licenses/R-COPYRIGHTS.txt`, and the GPL texts in `../../LICENSE` and
`../../licenses/GPL-3.txt`. The R runtime and package implementation source are
not included; these fixtures preserve data and attributed computed reference
outputs only.

The Mechler examples and the 1994 termination / Yong Lu regression data follow
upstream diptest tests; attribution is retained in NOTICE.md.
