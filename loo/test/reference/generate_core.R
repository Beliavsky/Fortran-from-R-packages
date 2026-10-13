# Run with R >= 3.5 and exact CRAN loo 2.10.1, posterior 1.7.0, matrixStats 1.5.0.
# From loo package root: Rscript test/reference/generate_core.R
stopifnot(packageVersion('loo') == '2.10.1')
suppressPackageStartupMessages(library(loo))
dir.create('test/reference/data', showWarnings=FALSE)
put <- function(name, x) write(ifelse(is.na(as.vector(x)), 'NaN', sprintf('%.17g', as.vector(x))), paste0('test/reference/data/', name, '.txt'), ncolumns=1)
i <- 1:200
ll <- cbind(-.3-.5*sin(i*.17)^2, -1-.7*cos(i*.113)-.002*i, -.015*i-.0004*i^2, -2+.05*sin(i*.031))
eff <- c(.5,1,.8,1.2)
for (method in c('psis','sis','tis')) {
  obj <- suppressWarnings(do.call(method, list(-ll, r_eff=eff)))
  put(paste0(method,'_lw'), weights(obj))
  put(paste0(method,'_k'), obj$diagnostics$pareto_k)
  put(paste0(method,'_ess'), obj$diagnostics$n_eff)
  res <- suppressWarnings(loo(ll,r_eff=eff,is_method=method))
  put(paste0(method,'_pointwise'),res$pointwise)
  put(paste0(method,'_estimates'),res$estimates)
}
obj <- suppressWarnings(psis(-ll, r_eff=eff))
put('waic_pointwise',suppressWarnings(waic(ll))$pointwise)
put('waic_estimates',suppressWarnings(waic(ll))$estimates)
put('elpd_pointwise',elpd(ll)$pointwise)
put('elpd_estimates',elpd(ll)$estimates)
old <- suppressWarnings(psislw(-ll))
put('old_lw',old$lw_smooth)
put('old_k',old$pareto_k)
logp <- -.002*i
logg <- -.003*i+.05*sin(i*.31)
ap <- suppressWarnings(loo_approximate_posterior(ll,log_p=logp,log_g=logg))
put('ap_pointwise',ap$pointwise)
put('ap_estimates',ap$estimates)
put('eff',relative_eff(exp(array(ll,c(50,4,4)))))
x <- matrix(sin(i*.17),200,4)+rep(0:3,each=200)
for (mode in c('mean','variance','sd','quantile')) {
   e <- E_loo(x,obj,type=mode,probs=c(.1,.5,.9),log_ratios=-ll)
   put(paste0('expect_',mode),e$value)
   put(paste0('expect_k_',mode),e$pareto_k)
}
put('mcse',mcse_loo(suppressWarnings(loo(ll,r_eff=eff)), threshold=10))
# Comparisons use exact same observations with three differing model pointwise scores.
models <- lapply(c(0,.1,.3), function(a) suppressWarnings(loo(ll+a*cos(i*.29), r_eff=eff)))
comp <- loo_compare(models)
put('compare_diff', comp$elpd_diff)
put('compare_se', comp$se_diff)
put('compare_p', comp$p_worse)
put('thin', loo:::.thin_draws(matrix(1:10,ncol=1),3))
cat('Core fixture files regenerated.\n')
perm <- c(10,1,3,7,2,9,5,4,8,6)
controlled <- function(fn, perm) {
  environment(fn) <- list2env(list(sample.int=function(n) perm,
                                 sample=function(x) perm[perm %in% x]), parent=environment(fn))
  fn
}
put('fold_random',controlled(loo::kfold_split_random,perm)(3,10))
put('fold_stratified',controlled(loo::kfold_split_stratified,perm)(3,c(-1,-1,8,8,-1,3,8,3,3,3)))
put('fold_grouped',controlled(loo::kfold_split_grouped,c(5,2,1,4,3))(3,c(8,-1,2,8,10,-1,2,1,10,1)))
