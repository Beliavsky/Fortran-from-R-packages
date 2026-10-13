# API coverage and omissions

Package status: substantial. Coverage basis: exported computational R functions.
The pinned NAMESPACE exports every named function through `exportPattern(".")`.
There are 14 distinct computational functions: all 14 are mapped (100%). The
fraction is not a claim of complete R compatibility. Presentation, plotting and
the sweetpotato dataset are excluded. `randtests.aux` is included conservatively
because it is exported, despite its internal documentation keyword.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
|---|---|---|---|
| `bartels.rank.test` | `randtests` | `bartels_rank_test` | substantial |
| `cox.stuart.test` | `randtests` | `cox_stuart_test` | substantial |
| `difference.sign.test` | `randtests` | `difference_sign_test` | substantial |
| `rank.test` | `randtests` | `rank_test` | substantial |
| `runs.test` | `randtests` | `runs_test` | substantial |
| `turning.point.test` | `randtests` | `turning_point_test` | substantial |
| `druns` | `randtests_distributions` | `druns` | substantial |
| `pruns` | `randtests_distributions` | `pruns` | substantial |
| `qruns` | `randtests_distributions` | `qruns` | substantial |
| `rruns` | `randtests_distributions` | `rruns` | substantial |
| `dbartelsrank` | `randtests_distributions` | `dbartelsrank` | substantial |
| `pbartelsrank` | `randtests_distributions` | `pbartelsrank` | substantial |
| `permut` | `randtests_permutations` | `permut`, `permut_apply` | partial |
| `randtests.aux` | `randtests_permutations` | `randtests_aux` | substantial |

## Mapping notes

- `bartels.rank.test`: Normal, beta, auto, and corrected exact p-values; exact mode requires untied n <= 9. R htest metadata is omitted.
- `cox.stuart.test`: Exact binomial tails, odd-middle deletion, and pair-tie removal; two-sided p-values capped at one. NaNs are omitted.
- `difference.sign.test`: Normal approximation with consecutive equal differences removed; typed result replaces R htest.
- `rank.test`: Pair-count rank statistic and original untied normal variance; no tie correction or R htest metadata.
- `runs.test`: Normal and exact methods; plotting and htest metadata omitted. Degenerate exact case returns p=1, two-sided values are capped at one.
- `turning.point.test`: Consecutive ties collapsed, original normal approximation; NaNs are consistently omitted, R htest omitted.
- `druns`: Elemental mass and log mass; positive integer counts with total <= 100000. Conformable arrays replace R recycling.
- `pruns`: Correct log CDF and direct upper tail; cutoff is floored rather than upstream rounding. Counts total <= 100000.
- `qruns`: Quantile inversion honors lower-tail and log-probability flags and endpoints; corrected upstream ignored flags.
- `rruns`: Inverse-CDF sampling with Fortran random_number; caller controls random_seed, no R RNG sequence parity.
- `dbartelsrank`: Exact permutation distribution for n=1..9, including corrected n=1 and n=2 behavior; no factorial-sized matrix.
- `pbartelsrank`: Exact permutation tails and log probabilities for n=1..9, corrected small-n callback behavior.
- `permut`: Finite real source vectors and optional pure scalar-real callback pointer; arbitrary R types, vector callback results and variadic arguments omitted. One million matrix entries maximum.
- `randtests.aux`: Bartels squared-difference and difference-sign increase statistics; invalid methods, nonfinite data or unrepresentable squared sums return NaN.

Every mapping names public implementing procedures. `randtests` re-exports the
submodule APIs for convenience. No formula/S3 interface, R object construction,
recycling or R RNG implementation is included. See README for error contracts,
resource limits, numerical caveats, and corrected upstream behavior.
