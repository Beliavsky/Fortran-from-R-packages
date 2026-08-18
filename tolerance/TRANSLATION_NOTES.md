# Translation notes

## Upstream

- R package: `tolerance`
- Upstream version: 3.0.0 (2024-04-18)
- Authors: Derek S. Young and Kedai Cheng
- Upstream license: GPL (>= 2)
- Upstream `NeedsCompilation`: no

The source package is pure R. The Fortran translation therefore ports the statistical algorithms themselves rather than wrapping an existing C/C++ backend.

## Coverage map

### Probability/distribution helpers

| R routine/file | Fortran equivalent |
|---|---|
| `F1` | `appell_f1` |
| `d2exp/p2exp/q2exp/r2exp` | same names in `tolerance_distributions` |
| `ddpareto/pdpareto/qdpareto/rdpareto` | same names |
| `dpoislind/ppoislind/qpoislind/rpoislind` | same names |
| `dnhyper/pnhyper/qnhyper/rnhyper` | same names |
| `dzipfman/pzipfman/qzipfman/rzipfman`, `zeta.fun` | `dzipfman/...`, `zeta_fun` |
| `ddiffprop/pdiffprop/qdiffprop/rdiffprop` | same names |
| `rwishart` | `rwishart_identity` |

For the beta-difference distribution, density/CDF evaluation uses the mathematically equivalent beta-convolution integral rather than relying on the Appell-F1 expression used in R. `appell_f1` remains available independently.

### Normal-family methods

| R routine | Fortran equivalent |
|---|---|
| `K.factor` | `k_factor` |
| `K.factor.sim` | `k_factor_sim` |
| `K.table` | `k_table` |
| `normtol.int` | `normtol_int` |
| `bayesnormtol.int` | `bayesnormtol_int` |
| `simnormtol.int` | `simnormtol_int` |
| `diffnormtol.int` | `diffnormtol_int` |
| numerical part of `norm.OC` | `norm_oc_content`, `norm_oc_alpha` |
| `norm.ss` | `norm_ss_fw`, `norm_ss_dir`, `norm_ss_ygzo` |

The HE/HE2/WBE/ELL/KM/EXACT/OCT K-factor paths are represented where applicable.

### Discrete tolerance intervals

- `bintol.int` -> `bintol_int`
- `poistol.int` -> `poistol_int`
- `negbintol.int` -> `negbintol_int`
- `hypertol.int` -> `hypertol_int`
- `neghypertol.int` -> `neghypertol_int`
- `umatol.int` -> `uma_upper`
- `acc.samp` -> `acceptance_sampling`

The method selectors from the R package are retained in the Fortran routines (for example CP/WS/AC/etc. for binomial intervals and the corresponding Poisson/negative-binomial method sets).

### Continuous tolerance intervals

- `exptol.int` -> `exptol_int`
- `exp2tol.int` -> `exp2tol_int`
- `uniftol.int` -> `uniftol_int`
- `laptol.int` -> `laptol_int`
- `gamtol.int` -> `gamtol_int`
- `logistol.int` -> `logistol_int`
- `cautol.int` -> `cautol_int`
- `paretotol.int` -> `paretotol_int`
- `exttol.int` -> `exttol_int`

Where upstream R delegated to `nlm`, `mle`, or related machinery, the translation uses the local Nelder-Mead/golden-section/Hessian implementation. The target likelihoods and parameterizations are retained.

### Fitted custom discrete models

- `dpareto.ll`, `dparetotol.int` -> `dpareto_mle`, `dparetotol_int`
- `poislind.ll`, `poislindtol.int` -> `poislind_mle`, `poislindtol_int`
- `zm.ll`, `zipftol.int` -> `zipf_mle`, `zipfman_mle`, `zeta_mle`, `zipftol_int`

### Nonparametric methods

- `npbetol.int` -> `npbeta_tol_int`
- `nptol.int` -> `nptol_int`
- exact two-row Hahn-Meeker cases -> `nptol_hm_options`
- Young-Mathew OS/FOS interpolation/extrapolation -> `nptol_ym_options`
- `np.order` -> `np_order`
- `distfree.est` -> scalar building blocks `distfree_sample_size`, `distfree_alpha`, `distfree_confidence`, `distfree_content`
- `interp`, `extrap`, `two.sided` -> implemented through `nptol_ym_options`

The scalar Fortran routines replace R's vectorized table-building behavior; callers can loop over parameter vectors directly.

### Regression and ANOVA

- `regtol.int` -> `regtol_int`
- `nonparregtolint` / `npregtol.int` -> `npregtol_int`
- `nonlinregtolint` / `nlregtol.int` -> `nlregtol_int`
- `anovatol.int` -> `anova_group_tol`
- `mvregtol.region` -> `mvregtol_region`

R formula objects and `lm`/`nls` objects are replaced by design matrices and, for nonlinear regression, a user-supplied model callback plus a numerical Jacobian.

### Multivariate tolerance regions

- `mvtol.region` -> `mvtol_region`
- KM/AM/GM/HM/MHM/V11/HM.V11/MC factors -> `mvtol_factor`
- `npmvtol.region` -> `npmvtol_region`

`npmvtol_region` accepts precomputed data-depth values rather than an R function object. Central and one-sided semispace coordinate types and floor/ceiling adjustment modes are preserved.

### Fiducial and semicontinuous procedures

- `fidbintol.int` -> `fidbintol_int`
- `fidpoistol.int` -> `fidpoistol_int`
- `fidnegbintol.int` -> `fidnegbintol_int`
- `semiconttol.int` -> `semiconttol_int`

The fiducial two-sample routines accept a scalar user callback in place of an R function closure.

### Bonferroni wrapper

- `bonftol.int` -> `bonf_tol_int` for a callback-driven equivalent
- `bonf_combine` combines two already-computed one-sided intervals.

### Deliberately omitted R-only code

- `plotly_anovatol`, `plotly_controltol`, `plotly_histtol`, `plotly_multitol`, `plotly_normOC`, `plotly_npmvtol`, `plotly_regtol`
- `plottol`
- package startup `.onAttach`
- `rFUN`, which deparses and rewrites R function bodies
- S3/data-frame printing, names, row names, and formula parsing

## Numerical implementation differences

1. The package is self-contained: regularized beta/gamma functions, central/noncentral CDFs and quantiles, integration, root-finding, small dense linear algebra, and RNGs are implemented in Fortran rather than delegated to R/MASS/stats4.
2. R's optimization functions are replaced by local Nelder-Mead, golden-section, numerical-Hessian, and explicit Newton/root iterations.
3. Regression APIs are numerical rather than formula-based.
4. Stochastic routines target the same probability laws but do not reproduce R's exact RNG stream.
5. Monte Carlo methods naturally vary from run to run; callers can increase their simulation counts where exposed.

## Validation

The source was compiled with GNU Fortran 14.2 in Fortran 2018 mode using `-Wall -Wextra -Wimplicit-interface -fcheck=all`. Ten test programs exercise the probability core, continuous/discrete/custom families, nonparametric methods, linear/nonlinear regression, multivariate regions, fiducial/semi-continuous routines, K-factor/OC calculations, Bonferroni callbacks, and Young-Mathew/Hahn-Meeker edge cases.

Independent numerical reference checks include:

- one-sided normal K factor: `K(20, alpha=.05, P=.95) = 2.3960016837521687`;
- multivariate AM factor for `n=20, p=2, alpha=.05, P=.90`: `7.382647443824329`;
- one-sided exponential tolerance-limit formulas checked against independent chi-square quantiles;
- polygamma values checked against standard special-function reference values.

FPM itself was not installed in the validation environment, so `fpm test` could not be executed there; `fpm.toml` is included and the same sources/tests were compiled and run directly with gfortran.
