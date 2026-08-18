# MCMCpack-fortran

Modern Fortran/FPM translation of computational routines from the R package **MCMCpack 1.7-1**.

The upstream package is by Andrew D. Martin, Kevin M. Quinn, Jong Hee Park, and contributors. This translation preserves the upstream GPL-3.0 license and attribution. The supplied Fortran translations of `coda`, `mcmc`, and `quantreg` are vendored under `vendor/` and retain their own notices and license files.

This is a computational translation: plotting, R formula parsing, S3 methods, interactive printing, R object serialization, and other presentation/R-runtime glue are intentionally omitted. Fortran interfaces use explicit arrays, priors, starting values, and control parameters.

## Implemented computational areas

The v0.2.0 tree contains the v0.1.0 translation plus the specialized samplers that were previously listed as gaps.

- Common random-number, probability, linear-algebra, truncated-normal, Wishart/inverse-Wishart, Dirichlet, inverse-gamma, inverse-Gaussian, and noncentral-hypergeometric utilities
- Conjugate Monte Carlo routines: binomial-beta, Poisson-gamma, normal-normal, and multinomial-Dirichlet
- Basic samplers: `MCMCregress`, `MCMCprobit`, `MCMClogit`, `MCMCpoisson`, `MCMCtobit`, `MCMCquantreg`, and generic `MCMCmetrop1R`
- Multinomial and count models: `MCMCmnl`, `MCMCnegbin`, and `MCMCnegbinChange`
- Ordered changepoint models: `MCMCbinaryChange`, `MCMCprobitChange`, `MCMCregressChange`, `MCMCpoissonChange`, `MCMCnegbinChange`, and `MCMCoprobitChange`
- Residual-break analysis computational model
- Hierarchical regression models: `MCMChregress`, `MCMChlogit`, and `MCMChpoisson`
- IRT and latent-variable models: `MCMCirt1d`, parameter-expanded/non-PX `MCMCirtHier1d` with Chib-style level-2 marginal likelihood, `MCMCirtKd`, `MCMCirtKdRob`, and `MCMCdynamicIRT1d`
- Ecological inference: `MCMChierEI` and `MCMCdynamicEI`
- Paired comparisons: `MCMCpaircompare`, `MCMCpaircompare2d`, and `MCMCpaircompare2dDP`
- Factor models: `MCMCfactanal`, `MCMCordfactanal`, and `MCMCmixfactanal`
- Ordinal probit and SVD regression: `MCMCoprobit` and `MCMCSVDreg`
- Variable-selection quantile regression: `SSVSquantreg`
- Panel hidden-Markov models: `HMMpanelFE` and `HMMpanelRE`
- Finite-truncation hierarchical/sticky count-state models: `HDPHMMpoisson`, `HDPHMMnegbin`, and explicit-duration `HDPHSMMnegbin`
- Model utilities: `BayesFactor` arithmetic, `PostProbMod`, break-list and transition-prior helpers, agreement matrices, SSVS marginal/top-model summaries, `vech`, `xpnd`, Procrustes alignment, and WAIC
- Conversion of sampler matrices to the translated `coda` chain type
- Quantile-regression starting values through the supplied `quantreg-fortran`
- Re-export of the supplied `mcmc-fortran` generic Metropolis interfaces

See `TRANSLATION_NOTES.md` for entry-point mappings and deliberate transition-kernel differences.

## Build with FPM

```sh
fpm build
fpm test
fpm run --example regression_demo
```

The package uses local path dependencies under `vendor/`, so no network access is required.

## Compiler validation

The complete source tree, including the three vendored dependencies, compiles with GNU Fortran 14.2.0 under ordinary free-form Fortran 2018 line-length rules. The validation build used runtime checks:

```sh
gfortran -std=f2018 -fcheck=all ...
```

All **22** supplied test programs pass in that build. The source does not require `-ffree-line-length-none`.

An `fpm` executable was not installed in the translation environment, so `fpm build`/`fpm test` could not be executed there; the FPM manifest and local dependency declarations are included and the equivalent full-tree GNU Fortran build was exercised directly.

## Minimal example

```fortran
program demo
   use mcmcpack, only : dp, set_seed, mcmc_result, mcmc_regress
   implicit none
   real(dp) :: x(6,2), y(6), beta_start(2), b0(2), b0prec(2,2)
   type(mcmc_result) :: fit

   x(:,1) = 1.0_dp
   x(:,2) = [-2.0_dp,-1.0_dp,0.0_dp,1.0_dp,2.0_dp,3.0_dp]
   y = 1.0_dp + 0.75_dp*x(:,2)
   beta_start = 0.0_dp
   b0 = 0.0_dp
   b0prec = 0.0_dp
   b0prec(1,1) = 0.01_dp
   b0prec(2,2) = 0.01_dp

   call set_seed(12345)
   fit = mcmc_regress(y,x,beta_start,b0,b0prec,2.0_dp,1.0_dp,100,1000,2)
   print *, sum(fit%draws,dim=1)/real(size(fit%draws,1),dp)
end program demo
```

The runnable version is `example/regression_demo.f90`.
