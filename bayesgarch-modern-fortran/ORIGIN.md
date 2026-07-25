# Origin and translation scope

This project translates the computational code in the attached R package:

- Package: `bayesGARCH`
- Version: 2.1.10
- Date: 2021-05-16
- Original author and copyright holder: David Ardia
- Original package license field: `GPL (>= 2)`
- Original computational sources: `R/sampler.R`, `R/functions.R`, and `src/fnGarchC.c`

The package license field means GPL version 2 or later, represented by the SPDX
identifier `GPL-2.0-or-later`. The complete GPL version 2 text from the original
package is retained as `LICENSE`. Every Fortran source file carries the original
copyright notice, SPDX identifier, and version 2-or-later grant.

## Mapping

| Original routine | Fortran implementation |
| --- | --- |
| `bayesGARCH`, `fn.bayesGARCH`, `fn.block` | `run_bayesgarch` in `bayesgarch_sampler` |
| `cvGARCH`, `fn.cvGARCH` | `garch11_filter` |
| `fn.alpha.full` | internal `draw_alpha_mh` |
| `fn.beta.full` | internal `draw_beta_mh` |
| `fn.w.full` | `draw_latent_scales` |
| `fn.nu.full` | internal `draw_nu_conditional` |
| `fn.post.garch` | `augmented_log_posterior` |
| `fn.neg.alpha` | `nu_root_function` |
| `fn.accept.proba` | `nu_log_acceptance` |
| `fn.filterAlpha` / `fnFilterAlphaC` | `filter_alpha` |
| `fnFilterAlphaAsymC` | `filter_alpha_asymmetric` |
| `fn.filterW` / `fnFilterWC` | `filter_w` |
| `fnFilterWAsymC` | `filter_w_asymmetric` |
| `fn.W` / `fnQDiffC` | `quasi_difference` |
| `fn.Dd.D` | `regression_posterior_1` and `regression_posterior_2` |
| `fnGarchC` filtering and simulation | general GARCH and threshold-GARCH filter/simulation routines |
| `formSmpl` | `form_posterior_sample` |

## Intentional exclusions

R-specific `mcmc`, `mcmc.list`, S3 behavior, `coda` integration, R list and
matrix metadata, console progress formatting, package registration, and plotting
used in examples are not translated. The supplied data set is not needed by the
numerical library and is not redistributed; a generated example CSV is supplied
instead.
