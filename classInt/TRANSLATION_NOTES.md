# Translation notes

## Upstream

- R package: `classInt`
- translated upstream version: 0.4-11
- upstream license: GPL (>= 2)
- upstream snapshot: `upstream/`

The objective is computational parity for numeric univariate class-interval
selection, not reproduction of R's object system or graphical presentation.

## Computational coverage

| Upstream functionality | Fortran status | Notes |
| --- | --- | --- |
| `classIntervals(..., style="fixed")` | translated | sorted caller-supplied breaks |
| `sd` | translated | sample standardization plus pretty/custom standardized breaks |
| `equal` | translated | equal-width classes |
| `pretty` | translated | native default pretty-break calculation |
| `quantile` | translated | R quantile types 1-9 |
| `kmeans` | translated | native 1-D Lloyd k-means with deterministic starts |
| `hclust` | translated | native 1-D hierarchy and recutting |
| `bclust` | translated | uses vendored `e1071-fortran` |
| `fisher` | translated | exact dynamic programming based on upstream `fish1.f` |
| `jenks` | translated | R dynamic-programming structure and right-closure semantics |
| `dpih` | translated | KernSmooth direct plug-in recursion and power-of-two FFT convolution |
| `headtails` | translated | iterative mean splitting and threshold stop rule |
| `maximum` | translated | largest adjacent gaps, including tied-gap behavior |
| `box` | translated | modern/legacy endpoint behavior and configurable IQR multiplier |
| `findCols` | translated | left/right interval closure and non-finite sentinel handling |
| `classify_intervals` | translated | existing-fit and fit-and-classify overloads |
| `getHclustClassIntervals` | translated | recuts retained hierarchy |
| `getBclustClassIntervals` | translated | recuts retained bagged hierarchy |
| `jenks.tests` | translated | GVF, TAI, optional area-weighted OAI |
| `logLik.classIntervals` | translated | Gaussian within-class likelihood; discrete class contribution 0 |
| `nPartitions` | translated | returns positive infinity above 170 unique observations |
| plotting, colors, printed tables | omitted | presentation-only |
| `classIntervals2shingle` | omitted | R `class`/shingle object integration |
| Date/POSIX/`units` metadata | omitted | R-specific class metadata |

## Fisher kernel parity

The upstream package contains the original fixed-form Fortran `fish1.f`.
During translation it was compiled separately and compared with
`fisher_exact` on:

```text
1, 2, 2, 4, 5, 8, 9, 10   with k = 3
```

Both implementations choose the same classes.  The upstream high-to-low class
statistics are reproduced to roundoff:

```text
min  max  mean                 population SD
8    10   9                    0.8164965809277289
4     5   4.5                  0.5
1     2   1.6666666666666667   0.4714045207910313
```

The only observed differences were at roughly 1e-15 in the standard deviations
because the modern code evaluates the equivalent variance expression in a
slightly different order.

## Jenks behavior

The R Jenks implementation deliberately differs from Fisher in its boundary
representation.  The Fortran port preserves that distinction.  For
`[0,0,0,1,2,50]` and three classes:

```text
Fisher breaks: 0, 0.5, 26, 50
Jenks breaks:  0, 0,   2,  50
```

Jenks forces right-closed intervals as upstream does.

## `dpih`

The `KernSmooth::dpih` algorithm is translated through its linear binning and
recursive `bkfe` kernel-functional estimates.  The default grid size is 401 and
the supported recursion levels are 0 through 5.  The kernel functional now uses
the same FFT layout as KernSmooth: the symmetric Hermite/Gaussian weights are
wrapped into a vector of length equal to the smallest power of two at least
`M + L + 1`, the binned counts are zero padded to that length, both sequences
are transformed, and the inverse transform of their product supplies the
finite convolution used in the weighted functional sum.

A retained direct-convolution implementation is used only as a validation
oracle.  Tests compare FFT and direct `dpih` results at all plug-in levels 1-5
on the default 401-point grid to a relative tolerance of `2e-12`, and also
exercise an 8193-point grid.  In a development benchmark using the strict-build
objects, the FFT path was about 3x faster at 401 points, 25x at 4097 points, and
33x at 8193 points.  Those timings are illustrative rather than contractual.

KernSmooth's source explicitly permits unlimited use and distribution.  The
algorithm is attributed in `NOTICE.md`.

## Clustering adaptations

`kmeans` and `bclust` use the deterministic RNG in `e1071-fortran`, not R's RNG,
so random initialization/resampling streams are intentionally not bit-identical.
The 1-D hclust implementation supports the most useful linkage methods for this
package: `complete`, `single`, `average`, `centroid`, `ward.d`, and `ward.d2`.
Exact tie ordering and less common R hclust methods are not claimed.

## Validation

`tools/run_strict_tests.sh` builds proxy, e1071, and classInt from source and
runs every test and example with GNU Fortran strict diagnostics.  It also runs
the source-rule checker independently for all three maintained source trees.

The test suite covers:

- upstream Fisher/Jenks small fixtures;
- direct original-`fish1.f` statistics;
- saved upstream Jenks log-likelihood result;
- fixed/equal/pretty/quantile/SD styles;
- k-means, hclust recutting, and bclust;
- head/tails, maximum, and box styles;
- `dpih` bandwidth, FFT/direct parity, large-grid FFT, and interval generation;
- GVF/TAI/OAI and AIC;
- fit-and-classify convenience API.
