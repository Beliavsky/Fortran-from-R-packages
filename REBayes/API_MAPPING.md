# API mapping

Upstream R name | Fortran entry point | v0.1.0 status
---|---|---
`KWDual` | `kw_dual`, `kw_fit` | implemented for NPMLE (`alpha=0`)
`KWPrimal` | `kw_primal`, `kw_fit` | implemented
`GLmix` | `glmix` | implemented
`Bmix` | `bmix` | implemented; R histogram/collapse is unnecessary for correctness
`B2mix` | `b2mix` | implemented; also returns marginal posterior means
`BPmix` | `bpmix` | implemented
`Pmix` | `pmix` | implemented
`NPmix` | `npmix` | implemented
`GVmix` | `gvmix` | implemented
`Gammamix` | `gammamix` | implemented
`Weibullmix` | `weibullmix` | implemented, censoring supported
`Gompertzmix` | `gompertzmix` | implemented
`TLmix` | `tlmix` | implemented
`Tncpmix` | `tncpmix` | implemented using deterministic quadrature
`Umix` | `umix` | implemented
`HLmix` | `hlmix` | implemented
`Cosslett` | `cosslett` | implemented
`GLVmix` | `glvmix` | implemented
`WGVmix` | `wgvmix` | implemented
`WGLVmix` | `wglvmix` | implemented
`WLVmix` | `wlvmix` | implemented
`WTLVmix` | `wtlvmix` | implemented
`predict.*mix` | `posterior_summary` plus fitted kernel/grid | computational core implemented
`Lfdr.*` | `lfdr_1d`; bivariate fitted `A` is exposed | computational core implemented
`KWsmooth` | `kw_smooth` | implemented
`KW2smooth` | use product smoothing on returned bivariate grid | low-level pieces implemented
`qKW` | `kw_quantiles` | implemented
`qKW2` | `kw2_marginal_quantiles` | implemented
`rKW` | `kw_random_sample` | implemented
`bwKW` | `bw_kw` | implemented
`bwKW2` | `bw_kw2` | implemented
`traprule` | `traprule_values` | implemented
`L1norm` | `l1_step_distance` | implemented from knot/value arrays
`dhuber` | `huber_pdf` | implemented
`hubereps` | `huber_eps` | implemented analytically
`ThreshFDR` | `thresh_fdr` | implemented
`Finv` | `finv` | implemented using procedure callback and bracket expansion
`RLR` | `rlr_fit` | implemented with ADMM/Newton, including general penalty matrix D
`medde` | `medde_fit` | implemented for Dorder 1--3 and Renyi/Shannon/log objectives
`qmedde` | `medde_quantiles` | implemented numerically
`rmedde` | `medde_random` | implemented (without R's optional extra smoothing noise)
`BDGLmix` | -- | deferred
`HuberSpline` | -- | deferred
`HodgesLehmann` | -- | deferred
`Rxiv` | -- | intentionally omitted, nonnumerical file-management utility
plot methods | -- | intentionally omitted

`mixture_fit` returns the support grid, probability masses, fitted mixture
density at observations, posterior mean, log likelihood, status, iteration
count and KKT gap.  `bivariate_mixture_fit` additionally returns both support
grids, flattened joint masses, posterior marginal means, posterior product
mean, and the likelihood matrix used by posterior/FDR calculations.
