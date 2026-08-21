# Translation coverage

BGFD has eight baseline distributions, each wrapped by the Bell-G transform and
the complementary Bell-G transform. The Fortran port implements the complete
numerical API for all sixteen resulting families.

| Baseline | Bell-G | Complementary Bell-G | Parameters before `lambda` |
| --- | --- | --- | --- |
| Exponential | BellE | CBellE | `alpha` |
| Exponentiated exponential | BellEE | CBellEE | `alpha, beta` |
| Weibull | BellW | CBellW | `alpha, beta` |
| Exponentiated Weibull | BellEW | CBellEW | `alpha, beta, theta` |
| Fisk / log-logistic | BellF | CBellF | `a, b` |
| Lomax | BellL | CBellL | `b, q` |
| Burr XII | BellB | CBellB | `a, b, k` |
| Burr X | BellBX | CBellBX | `a` |

For every family the following operations are translated:

- density and log-density (`d_*`)
- CDF, upper tail, and log probabilities (`p_*`)
- quantile with lower/upper-tail and log-probability options (`q_*`)
- random generation (`r_*`)
- survival/log-survival (`s_*`)
- hazard/log-hazard (`h_*`)
- MLE and goodness-of-fit wrapper (`m_*`)

The generic `bgfd_core` API also exposes a family-ID based interface so code can
select distributions at run time without a large `select case` in user code.

## Intentionally omitted

The upstream package imports graphics because `AdequacyModel` can plot TTT and
other diagnostics, but BGFD itself contains no substantive custom plotting
algorithm. No graphics/device code is included in the Fortran library.

R-specific list/matrix formatting in the `m*` return values is replaced by the
`bgfd_fit_result` derived type.
