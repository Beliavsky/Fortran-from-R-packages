# API reference

Import from `diptest_api`. `dp` is the sibling `r_kinds` kind.

## Status codes

| Name | Value | Meaning |
|---|---:|---|
| `dip_ok` | 0 | Successful computation |
| `dip_empty_sample` | 1 | No observations remain |
| `dip_unsorted` | 2 | A decreasing pair violated `sorted=.true.` |
| `dip_nonfinite` | 3 | NaN/Infinity for `dip`; Infinity after NaN removal for a test |
| `dip_invalid_simulation` | 4 | Wrong matrix shape, no replicates, nonfinite/out-of-range null draws, or nonpositive B |
| `dip_invalid_argument` | 5 | Invalid statistic or sample size supplied directly to calibration |

Invalid computations return a quiet NaN statistic or p-value as appropriate.
For `dip_compute`, `n` and `status` are always defined; invalid-result arrays are
unallocated and modal indices remain zero. For an invalid null matrix the
observed statistic remains available in `dip_test_result`.

## `dip` and `dip_compute`

`dip(x, min_is_zero)` is a pure scalar function. `x` is rank-one `real(dp)` of
positive length and finite values. The optional logical `min_is_zero` defaults
to false. Input is never modified. Error reporting for this short form is NaN;
use `dip_compute` for explicit status.

`call dip_compute(x, result, min_is_zero, sorted)` additionally returns complete
`type(dip_result)` diagnostics. Optional `sorted` defaults false. When true,
nondecreasing input is required; no silent sorting occurs. Ties are valid.

| Component | Meaning |
|---|---|
| `statistic` | Hartigan's dip |
| `n`, `status` | Number of input observations and status code |
| `x(:)` | Sorted original-coordinate copy of the data |
| `low`, `high` | Inclusive 1-based modal interval endpoints |
| `lower`, `upper` | `x(low)`, `x(high)` |
| `gcm(:)`, `lcm(:)` | Final modal-interval convex/concave hull nodes |
| `full_gcm(:)`, `full_lcm(:)` | Full-sample convex/concave hull nodes |
| `mn(:)`, `mj(:)` | Full-sample convex predecessor / concave successor arrays |

All node indices refer to sorted `x`. GCM runs high-to-low, LCM low-to-high.
For a singleton each hull is `[1]`. For a constant sample of size n>1 each GCM
is `[n,1]`, each LCM `[1,n]`, `mn(:)=1`, and `mj(:)=n`.

## `dip_test`

`call dip_test(x, result, uniform_samples)` is pure. The result has type
`dip_test_result`, with real components `statistic` and `p_value`, integer
components `n`, `status`, `replicates`, and logical `simulated`, `asymptotic`.
`n` counts observations after NaN omission. A supplied null sample matrix must
have shape `(n,B)` for positive B; each column is one sample and all entries must
be finite in `[0,1]`. Endpoints are accepted for deterministic fixtures even
though an ideal continuous generator hits them with probability zero.

Without the optional matrix the original table is used, `replicates=0`, and
`simulated=.false.`. With valid null samples, `simulated=.true.` and
`replicates=B`. An invalid null matrix yields `dip_invalid_simulation`; its
requested simulation metadata are retained. If `n<=3`, p=1 and all simulation
options are ignored; metadata keep their default false/zero values.

## `dip_test_simulate`

`call dip_test_simulate(x, result, replicates)` uses the intrinsic RNG and is
impure. `replicates` is an optional positive integer, default 2000. It accepts the
same observed data as `dip_test`; seed externally with `random_seed`. Null
samples are streamed one column at a time. There is no global RNG maintained by
this package. A nonpositive B is ignored only when n<=3, matching R's small-n
short circuit. No cross-compiler or R-seed equivalence is claimed.

## `dip_p_value`

`call dip_p_value(statistic, n, p_value, status, asymptotic)` is pure elemental. A finite
statistic in `[0,0.5]` and positive integer n are required; the wider bound
accommodates the n=1 classical floor. Optional `status` and `asymptotic` are
outputs. `asymptotic` is true only for `n>72000`. Calibration is available
independently of computing the statistic; no significance level is required.

## Complexity

Sorting reuses `rfortran-core`. The corrected hull algorithm uses O(n)
auxiliary storage and data-dependent iterative hull passes. A caller-provided
Monte Carlo matrix occupies O(nB) caller storage; the intrinsic generator
requires O(n) temporary storage. Calibration has a fixed 21-by-26 table.
