# Translation notes

## Upstream target

- R package: `multcomp`
- upstream version: 1.4-32
- license: GPL-2
- source snapshot: `upstream/`
- numerical dependency: vendored `mvtnorm-fortran` 1.4.2

The upstream package has no compiled source. Its main R implementation delegates
multivariate normal/t probabilities and quantiles to `mvtnorm`; this translation
preserves the same separation.

## Computational coverage

Translated directly or represented by native typed APIs:

- `parm` / `modelparm` numerical representation;
- general matrix `glht` hypotheses and identity hypotheses;
- `univariate`, `adjusted`, `Ftest`, and `Chisqtest` calculations;
- adjusted and univariate confidence critical values;
- single-step, free, Shaffer, Westfall, and every upstream `p.adjust` method;
- `contrMat` families;
- the `maxsets` closed-testing combinatorics;
- compact-letter display insert-absorb calculations;
- multiple marginal model covariance assembly and block-diagonal linear
  functions;
- coefficient-subset testing corresponding to `cftest`.

## Deliberate interface adaptations

The following are R object-system/model-description facilities rather than
standalone numerical algorithms, and are not reproduced literally:

- formula, expression, and character parsing for hypotheses;
- `mcp()` inspection of R factor contrasts and design-matrix `assign` metadata;
- automatic coefficient/covariance extraction from `lm`, `glm`, `coxph`, `lme`,
  `merMod`, `fixest`, `glmmTMB`, `polr`, and similar model objects;
- S3 `print`, `plot`, `summary`, `coef`, `vcov`, and `confint` dispatch;
- response-data extraction used only to order/display CLD letters.

The native Fortran API accepts the corresponding numerical objects directly.
This makes the statistical core usable with coefficients and covariance matrices
from any Fortran model-fitting library.

## Singular contrast correlations

All-pair contrast families such as Tukey are linearly dependent, so their
correlation matrices are positive semidefinite rather than positive definite.
The earlier `mvtnorm-fortran` backend used a strict Cholesky integrator. Before
calling that integrator, `multcomp-fortran` detects rank deficiency and applies
a `1e-10` diagonal-nugget-equivalent shrinkage to off-diagonal correlations.
The covariance/correlation stored and returned to callers is not changed.

This is a numerical integration stabilization only. It reproduces the upstream
cholesterol Shaffer/Westfall and warpbreaks Tukey fixtures within the printed R
precision. The tiny perturbation means singular max-t probabilities are not
claimed bit-for-bit identical to R `mvtnorm` in the last digits.

## mvtnorm dependency safety patch

The vendored `mvtnorm-fortran` quantile objective used `huge(real(dp))` as a
stand-in for one-sided infinite integration bounds. Under strict floating-point
traps, downstream Student-t calculations can overflow while squaring such a
number. Those internal bounds are changed to `+/-1e100_dp`, which is effectively
infinite for the supported probability calculations without triggering overflow.

## Closed-testing size

The native `maxsets` enumerator uses a signed 64-bit bit mask and accepts at most
60 hypotheses. Closed-testing set enumeration is exponential in the number of
hypotheses in any case. Single-step, free, marginal, and ordinary `p.adjust`
methods do not have this limit.

## Reproducibility

The translated `mvtnorm` randomized quasi-Monte-Carlo backend has an explicit
seed in `probability_control` (default 12345), so default Fortran results are
repeatable. They need not follow R's RNG stream.

## Source conventions

Maintained `multcomp-fortran` source uses one public `dp = real64` kind. Every
dummy argument has explicit `INTENT`/`VALUE`, is declared on its own line, and
has a meaningful trailing FORD `!!` documentation comment. Maintained source is
kept within 132 columns and is written to be compatible with `fprettify`.
