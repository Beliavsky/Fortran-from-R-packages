# Translation notes

## Upstream

- Package: RMKdiscrete
- Version: 0.1
- Date in DESCRIPTION: 2014-10-17
- Author/Maintainer: Robert M. Kirkpatrick
- Upstream license: GPL (>= 2)
- Upstream implementation: approximately 1,000 lines of R and 800 lines of C

The translation preserves the upstream GPL-2.0-or-later licensing and attribution.

## API coverage

| Upstream R function | Fortran API |
|---|---|
| `LGPMVP` | `lgp_from_mu_sigma2`, `lgp_from_theta_lambda`, `lgp_from_mu_theta`, `lgp_from_sigma2_lambda`, `lgp_from_sigma2_theta`, `lgp_from_mu_lambda` |
| `LGP.findmax` | `lgp_findmax` |
| `LGP.get.nc` | `lgp_get_nc` |
| `dLGP` | `dlgp` |
| `pLGP` | `plgp` |
| `qLGP` | `qlgp` |
| `rLGP` | `rlgp`, `rlgp_sample` |
| `sLGP` | `slgp`, returning `type(lgp_summary)` |
| `dbiLGP` | `dbilgp` |
| `biLGP.logMV` | `bilgp_logmv` |
| `rbiLGP` | `rbilgp`, `rbilgp_sample` |
| `dnegbin` | `dnegbin` |
| `negbinMVP` | `negbin_from_nu_p`, `negbin_from_mu_sigma2`, `negbin_from_mu_nu`, `negbin_from_mu_p`, `negbin_from_sigma2_p` |
| `dbinegbin` | `dbinegbin` |
| `binegbin.logMV` | `binegbin_logmv` |
| `rbinegbin` | `rbinegbin`, `rbinegbin_sample` |
| `dmanaclash.dmg` | `dmanaclash_dmg` |
| `dmanaclash.xyN` | `dmanaclash_xyn` |
| `dmanaclash.net` | `dmanaclash_net` |
| `rmanaclash` | `rmanaclash`, `rmanaclash_sample` |

Internal R helpers used only for vector recycling, row sums, or `.C` dispatch are
not reproduced as public procedures.

## Numerical implementation

### Lagrangian Poisson

The upstream PMF is translated directly:

`theta * (theta + lambda*x)^(x-1) * exp(-theta-lambda*x) / x!`

For negative `lambda`, the support is finite and the same upstream support rule is
used. The upstream C code optionally accumulates probabilities in 21 magnitude
bins (`carefulprobsum`). The Fortran translation instead uses Kahan compensated
summation throughout. This is simpler and at least as accurate for binary64.

The upstream C implementation may stop normalizer accumulation early for extremely
large finite supports. The Fortran implementation sums the complete finite support
unless the caller supplies a positive tolerance for an exceptionally large support.

CDF and quantile calculations use direct accumulated PMFs. The R package contains
optimized special paths for vectors of probabilities/quantiles; those are R-level
performance optimizations rather than distinct statistical algorithms, so the
Fortran API is scalar and callers can use ordinary Fortran loops or array helpers.

`qlgp` uses conventional lower-tail quantile semantics consistently. This also
matches the upstream vector-specialized path. The upstream scalar `lambda=0` R
branch appears to call `qpois(..., lower.tail=FALSE)` after already transforming
`p`; that inconsistent special-case behavior is not reproduced.

### Random generation

The LGP generators retain the algorithms used by the R frontend:

- Poisson special case for `lambda=0`
- normal approximation in the same high-`theta` regions
- finite-support inversion for negative `lambda`
- branching-process generation for positive `lambda`

Poisson, gamma, negative-binomial, normal, and binomial RNG primitives are
implemented locally so the package has no Rmath or external-library dependency.
Negative-binomial RNG uses the exact gamma-Poisson mixture.

### Bivariate laws

The bivariate LGP and negative-binomial PMFs are direct convolutions over the
shared latent component, as in the upstream C code. Equal-parameter convolution
shortcuts are retained where appropriate. Log-transformed means, variances, and
covariances use the same convergence idea as the upstream routines.

### Mana Clash

The trivariate `(x,y,N)` PMF, conditional/marginal damage PMFs, net-damage PMF,
and simulation logic are direct translations of the R formulas. As upstream,
the density functions require strictly positive four-outcome weights; the random
generator permits zero weights for the first three outcomes provided the stopping
outcome has positive weight.

## Fortran-specific API choices

- Counts are integers, so invalid non-integer count inputs are prevented by type.
- R's vector recycling is intentionally not reproduced.
- Optional `log`, `lower.tail`, and `log.p` behavior is represented by optional
  logical arguments.
- The overloaded R parameter-conversion helpers are split into explicit Fortran
  procedures, avoiding ambiguous optional-output interfaces.
- `dp` is defined as `kind(1.0d0)`.
- Source is free-form Fortran 2018, uses `implicit none`, stays within 132 columns,
  and avoids multiple statements per line.

## Validation

The included tests cover:

1. positive- and negative-`lambda` LGP PMF/CDF/quantile values;
2. finite-support normalization and all LGP parameter-conversion directions;
3. negative-binomial PMFs in all three upstream parameterizations and conversion helpers;
4. bivariate LGP and bivariate negative-binomial PMFs, including independent special cases;
5. log-moment covariance in independent bivariate cases;
6. Mana Clash joint, conditional, marginal, and net probabilities;
7. RNG support and structural constraints.

Fixed reference probabilities were independently computed from the defining PMFs,
not from the Fortran routines under test.
