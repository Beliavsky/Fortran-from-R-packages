# Changelog

## 0.2.0

- Added `MCMChierEI` and `MCMCdynamicEI` computational samplers.
- Added `MCMCdynamicIRT1d` and robust multidimensional `MCMCirtKdRob`.
- Added truncated-DP `MCMCpaircompare2dDP`.
- Added `SSVSquantreg` asymmetric-Laplace variable-selection sampler.
- Added `HMMpanelFE` and `HMMpanelRE` panel hidden-Markov samplers.
- Added finite-truncation `HDPHMMpoisson`, `HDPHMMnegbin`, and explicit-duration `HDPHSMMnegbin`.
- Extended `MCMCirtHier1d` with parameter expansion and Chib-style level-2 marginal-likelihood output.
- Added five specialized-family test programs; complete suite is now 22 tests.
- Wrapped all Fortran source to the standard free-form 132-column limit.

## 0.1.0

- Initial translation of the core MCMCpack computational routines and major standalone samplers.
