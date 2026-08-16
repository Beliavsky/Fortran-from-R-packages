# Rfast API map

Upstream version: **2.1.5.2**. The supplied NAMESPACE contains **444 unique exported names**. Rfast is library-scale; v0.3.0 ports the reusable computational core and records remaining advanced routines as pending rather than claiming false parity.

Status counts: **257 ported/consolidated**, **1 partial**, **13 direct Fortran-intrinsic equivalents**, **28 R-only omissions**, **145 pending computational exports**.

| R export | v0.3 status | Fortran mapping / note |
|---|---|---|
| `Elem<-` | r-only | `omitted: R object/package infrastructure` |
| `AddToNamespace` | r-only | `omitted: R object/package infrastructure` |
| `Choose` | ported | `choose_real` |
| `Crossprod` | ported | `crossprod` |
| `Diag.fill` | pending | `not in v0.3.0 computational-core port` |
| `Diag.matrix` | pending | `not in v0.3.0 computational-core port` |
| `Digamma` | ported | `digamma_r` |
| `Dist` | ported | `dist_matrix / vecdist_upper / total_dist` |
| `Elem` | r-only | `omitted: R object/package infrastructure` |
| `Hash` | r-only | `omitted: R object/package infrastructure` |
| `Hash.key.multi` | r-only | `omitted: R object/package infrastructure` |
| `Lbeta` | ported | `log_beta` |
| `Lchoose` | ported | `log_choose` |
| `Lgamma` | ported | `Fortran intrinsic log_gamma` |
| `Log` | ported | `Fortran intrinsic log` |
| `Mad` | ported | `mad_r` |
| `Match` | pending | `not in v0.3.0 computational-core port` |
| `Median` | ported | `median_r` |
| `Norm` | pending | `not in v0.3.0 computational-core port` |
| `Order` | ported | `order_real` |
| `Outer` | ported | `outer_vec` |
| `Pmax` | ported | `pmax_vec` |
| `Pmin` | ported | `pmin_vec` |
| `Pmin_Pmax` | ported | `pmin_pmax_vec` |
| `Rank` | ported | `rank_average` |
| `RemoveFromNamespace` | r-only | `omitted: R object/package infrastructure` |
| `Rnorm` | ported | `rnorm_vector` |
| `Round` | intrinsic | `nint / anint` |
| `Sort` | ported | `sort_real` |
| `Sort.int` | ported | `sort_integer` |
| `Stack` | r-only | `omitted: R object/package infrastructure` |
| `Table` | pending | `not in v0.3.0 computational-core port` |
| `Table.sign` | pending | `not in v0.3.0 computational-core port` |
| `Tcrossprod` | ported | `tcrossprod` |
| `Trigamma` | ported | `trigamma_r` |
| `Unique` | ported | `unique_real / unique_int` |
| `Var` | ported | `variance_r` |
| `XopY.sum` | pending | `not in v0.3.0 computational-core port` |
| `acg.mle` | pending | `not in v0.3.0 computational-core port` |
| `allbetas` | pending | `not in v0.3.0 computational-core port` |
| `allttests` | pending | `not in v0.3.0 computational-core port` |
| `ancova1` | pending | `not in v0.3.0 computational-core port` |
| `ancovas` | pending | `not in v0.3.0 computational-core port` |
| `anova1` | ported | `one_way_anova` |
| `anova_propreg` | pending | `not in v0.3.0 computational-core port` |
| `anova_qpois.reg` | pending | `not in v0.3.0 computational-core port` |
| `anova_quasipois.reg` | pending | `not in v0.3.0 computational-core port` |
| `anovas` | ported | `one_way_anovas` |
| `apply.condition` | r-only | `omitted: R object/package infrastructure` |
| `ar1` | ported | `ar1_fit` |
| `as.Rfast.function` | r-only | `omitted: R object/package infrastructure` |
| `as_integer` | intrinsic | `int intrinsic` |
| `auc` | ported | `auc_score` |
| `bc` | pending | `not in v0.3.0 computational-core port` |
| `bcdcor` | pending | `not in v0.3.0 computational-core port` |
| `beta.mle` | ported | `beta_mle` |
| `betabinom.mle` | pending | `not in v0.3.0 computational-core port` |
| `betageom.mle` | pending | `not in v0.3.0 computational-core port` |
| `betaprime.mle` | ported | `betaprime_mle` |
| `bic.corfsreg` | ported | `bic_corfsreg` |
| `bic.fs.reg` | pending | `not in v0.3.0 computational-core port` |
| `binary_search` | ported | `binary_search_real` |
| `bincomb` | pending | `not in v0.3.0 computational-core port` |
| `binom.mle` | ported | `binomial_mle (known N)` |
| `block.anova` | pending | `not in v0.3.0 computational-core port` |
| `block.anovas` | pending | `not in v0.3.0 computational-core port` |
| `boot.ttest2` | pending | `not in v0.3.0 computational-core port` |
| `borel.mle` | ported | `borel_mle` |
| `bs.reg` | pending | `not in v0.3.0 computational-core port` |
| `btmprobs` | pending | `not in v0.3.0 computational-core port` |
| `cat.goftests` | pending | `not in v0.3.0 computational-core port` |
| `cauchy.mle` | ported | `cauchy_mle` |
| `checkAliases` | r-only | `omitted: R object/package infrastructure` |
| `checkExamples` | r-only | `omitted: R object/package infrastructure` |
| `checkNamespace` | r-only | `omitted: R object/package infrastructure` |
| `checkTF` | r-only | `omitted: R object/package infrastructure` |
| `checkUsage` | r-only | `omitted: R object/package infrastructure` |
| `check_data` | pending | `not in v0.3.0 computational-core port` |
| `chi2Test` | pending | `not in v0.3.0 computational-core port` |
| `chi2Test_univariate` | pending | `not in v0.3.0 computational-core port` |
| `chi2tests` | pending | `not in v0.3.0 computational-core port` |
| `chisq.mle` | ported | `chisq_mle` |
| `cholesky` | ported | `cholesky_lower` |
| `circlin.cor` | ported | `circlin_cor` |
| `coeff` | pending | `not in v0.3.0 computational-core port` |
| `col.coxpoisrat` | pending | `not in v0.3.0 computational-core port` |
| `col.yule` | pending | `not in v0.3.0 computational-core port` |
| `colAll` | ported | `colall` |
| `colAny` | ported | `colany` |
| `colCountValues` | pending | `not in v0.3.0 computational-core port` |
| `colCumMaxs` | ported | `cumulative_max per column` |
| `colCumMins` | ported | `cumulative_min per column` |
| `colCumProds` | ported | `cumulative_prod per column` |
| `colCumSums` | ported | `cumulative_sum per column` |
| `colFalse` | ported | `colfalse` |
| `colMads` | ported | `colmads` |
| `colMaxs` | ported | `colmaxs` |
| `colMedians` | ported | `colmedians` |
| `colMins` | ported | `colmins` |
| `colMinsMaxs` | pending | `not in v0.3.0 computational-core port` |
| `colOrder` | ported | `colorder` |
| `colPmax` | pending | `not in v0.3.0 computational-core port` |
| `colPmin` | pending | `not in v0.3.0 computational-core port` |
| `colRanks` | ported | `colranks` |
| `colShuffle` | pending | `not in v0.3.0 computational-core port` |
| `colSort` | ported | `colsort` |
| `colTabulate` | pending | `not in v0.3.0 computational-core port` |
| `colTrue` | ported | `coltrue` |
| `colTrueFalse` | ported | `coltruefalse` |
| `colVars` | ported | `colvars` |
| `colanovas` | ported | `one_way_anovas` |
| `colar1` | ported | `ar1_fit per column` |
| `colaucs` | ported | `column_aucs` |
| `colcvs` | ported | `colcvs` |
| `coldiffs` | pending | `not in v0.3.0 computational-core port` |
| `colexp2.mle` | ported | `colexponential2_mle` |
| `colexpmle` | ported | `colexponential_mle` |
| `colgammamle` | ported | `colgamma_mle` |
| `colgeom.mle` | ported | `colgeom_mle` |
| `colhameans` | ported | `colhameans` |
| `colinvgauss.mle` | ported | `colinvgauss_mle` |
| `colkurtosis` | ported | `colkurtosis` |
| `collaplace.mle` | ported | `collaplace_mle` |
| `collindley.mle` | ported | `collindley_mle` |
| `colmaxboltz.mle` | ported | `colmaxboltz_mle` |
| `colmeans` | ported | `colmeans` |
| `colnormal.mle` | ported | `colnormal_mle` |
| `colnormlog.mle` | ported | `collognormal_mle` |
| `colnth` | ported | `colnth` |
| `colpareto.mle` | ported | `colpareto_mle` |
| `colpois.tests` | pending | `not in v0.3.0 computational-core port` |
| `colpoisdisp.tests` | pending | `not in v0.3.0 computational-core port` |
| `colpoisson.anovas` | pending | `not in v0.3.0 computational-core port` |
| `colpoisson.mle` | ported | `colpoisson_mle` |
| `colprods` | ported | `colprods` |
| `colquasipoisson.anovas` | pending | `not in v0.3.0 computational-core port` |
| `colrange` | ported | `colranges` |
| `colrayleigh.mle` | ported | `colrayleigh_mle` |
| `colrint.regbx` | ported | `colrint_regbx` |
| `colrow.value` | pending | `not in v0.3.0 computational-core port` |
| `colskewness` | ported | `colskewness` |
| `colsums` | ported | `colsums` |
| `columns` | intrinsic | `array section x(:,idx)` |
| `colvarcomps.mle` | ported | `colvarcomps_mle (balanced and unbalanced groups)` |
| `colvarcomps.mom` | ported | `varcomps_mom` |
| `colvm.mle` | pending | `not in v0.3.0 computational-core port` |
| `colwatsons` | pending | `not in v0.3.0 computational-core port` |
| `colweibull.mle` | ported | `colweibull_mle` |
| `comb_n` | ported | `combination_count` |
| `cor.fbed` | pending | `not in v0.3.0 computational-core port` |
| `cor.fsreg` | pending | `not in v0.3.0 computational-core port` |
| `cora` | ported | `correlation_matrix` |
| `corpairs` | ported | `correlation_pairs` |
| `correls` | ported | `column_correlations` |
| `count_value` | ported | `count_value_real / count_value_int` |
| `cova` | ported | `covariance_matrix` |
| `cox.poisrat` | pending | `not in v0.3.0 computational-core port` |
| `cqtest` | pending | `not in v0.3.0 computational-core port` |
| `cqtests` | pending | `not in v0.3.0 computational-core port` |
| `ct.mle` | pending | `not in v0.3.0 computational-core port` |
| `data.frame.to_matrix` | r-only | `omitted: R object/package infrastructure` |
| `dcor` | ported | `distance_correlation` |
| `dcor.ttest` | pending | `not in v0.3.0 computational-core port` |
| `dcov` | ported | `distance_covariance` |
| `design_matrix` | pending | `not in v0.3.0 computational-core port` |
| `diri.nr2` | ported | `dirichlet_mle` |
| `dirimultinom.mle` | pending | `not in v0.3.0 computational-core port` |
| `dirknn` | pending | `not in v0.3.0 computational-core port` |
| `dirknn.cv` | pending | `not in v0.3.0 computational-core port` |
| `dista` | ported | `dista_matrix` |
| `dmvnorm` | ported | `dmvnorm` |
| `dmvt` | ported | `dmvt` |
| `dvar` | ported | `distance_variance` |
| `eachcol.apply` | r-only | `omitted: R object/package infrastructure` |
| `eachrow` | r-only | `omitted: R object/package infrastructure` |
| `edist` | ported | `energy_distance` |
| `eel.test1` | ported | `eel_test1` |
| `eel.test2` | ported | `eel_test2` |
| `eigen.sym` | ported | `eigen_sym_jacobi` |
| `el.test1` | ported | `el_test1` |
| `el.test2` | ported | `el_test2` |
| `env.copy` | r-only | `omitted: R object/package infrastructure` |
| `exact.ttest2` | pending | `not in v0.3.0 computational-core port` |
| `exp2.mle` | ported | `exponential2_mle` |
| `expmle` | ported | `exponential_mle` |
| `expregs` | pending | `not in v0.3.0 computational-core port` |
| `fish.kent` | pending | `not in v0.3.0 computational-core port` |
| `floyd` | ported | `floyd_warshall` |
| `foldnorm.mle` | pending | `not in v0.3.0 computational-core port` |
| `freq.max` | pending | `not in v0.3.0 computational-core port` |
| `freq.min` | pending | `not in v0.3.0 computational-core port` |
| `fs.reg` | pending | `not in v0.3.0 computational-core port` |
| `ftest` | ported | `ftest_variance` |
| `ftests` | ported | `column_ftests` |
| `g2Test` | pending | `not in v0.3.0 computational-core port` |
| `g2Test_perm` | pending | `not in v0.3.0 computational-core port` |
| `g2Test_univariate` | pending | `not in v0.3.0 computational-core port` |
| `g2Test_univariate_perm` | pending | `not in v0.3.0 computational-core port` |
| `g2tests` | pending | `not in v0.3.0 computational-core port` |
| `g2tests_perm` | pending | `not in v0.3.0 computational-core port` |
| `gammacon` | pending | `not in v0.3.0 computational-core port` |
| `gammamle` | ported | `gamma_mle` |
| `gammanb` | ported | `fit_gamma_nb` |
| `gammanb.pred` | ported | `predict_gamma_nb` |
| `gammareg` | ported | `gamma_regression` |
| `gammaregs` | ported | `gamma_regs` |
| `gaussian.nb` | ported | `fit_gaussian_nb` |
| `gaussiannb.pred` | ported | `predict_gaussian_nb` |
| `gchi2Test` | pending | `not in v0.3.0 computational-core port` |
| `geom.anova` | pending | `not in v0.3.0 computational-core port` |
| `geom.anovas` | pending | `not in v0.3.0 computational-core port` |
| `geom.mle` | ported | `geometric_mle` |
| `geom.nb` | ported | `fit_geometric_nb` |
| `geom.regs` | pending | `not in v0.3.0 computational-core port` |
| `geomnb.pred` | ported | `predict_geometric_nb` |
| `gini` | ported | `gini_r` |
| `ginis` | pending | `not in v0.3.0 computational-core port` |
| `glm_logistic` | ported | `glm_logistic` |
| `glm_poisson` | ported | `glm_poisson` |
| `group` | pending | `not in v0.3.0 computational-core port` |
| `group.sum` | pending | `not in v0.3.0 computational-core port` |
| `groupcorrels` | pending | `not in v0.3.0 computational-core port` |
| `gumbel.mle` | ported | `gumbel_mle` |
| `halfnorm.mle` | ported | `halfnormal_mle` |
| `hash.find` | r-only | `omitted: R object/package infrastructure` |
| `hash.list` | r-only | `omitted: R object/package infrastructure` |
| `hash2list` | r-only | `omitted: R object/package infrastructure` |
| `hd.eigen` | pending | `not in v0.3.0 computational-core port` |
| `hsecant01.mle` | pending | `not in v0.3.0 computational-core port` |
| `iag.mle` | pending | `not in v0.3.0 computational-core port` |
| `ibeta.mle` | pending | `not in v0.3.0 computational-core port` |
| `invdir.mle` | ported | `inverse_dirichlet_mle` |
| `invgauss.mle` | ported | `invgauss_mle` |
| `invgauss.reg` | ported | `invgauss_regression` |
| `invgauss.regs` | ported | `invgauss_regs` |
| `is.symmetric` | ported | `is_symmetric_matrix` |
| `is_element` | ported | `is_element_real` |
| `is_integer` | ported | `is_integer_value` |
| `iterator` | r-only | `omitted: R object/package infrastructure` |
| `james` | pending | `not in v0.3.0 computational-core port` |
| `knn` | ported | `knn_classify / knn_regress` |
| `knn.cv` | pending | `not in v0.3.0 computational-core port` |
| `kruskaltest` | ported | `kruskal_test` |
| `kruskaltests` | pending | `not in v0.3.0 computational-core port` |
| `kuiper` | ported | `kuiper_test` |
| `kurt` | ported | `kurtosis_r` |
| `kurt.test2` | pending | `not in v0.3.0 computational-core port` |
| `laplace.mle` | ported | `laplace_mle` |
| `lindley.mle` | ported | `lindley_mle` |
| `list.ftests` | pending | `not in v0.3.0 computational-core port` |
| `lmfit` | ported | `lmfit` |
| `logcauchy.mle` | ported | `logcauchy_mle` |
| `logistic.cat1` | pending | `not in v0.3.0 computational-core port` |
| `logistic.mle` | ported | `logistic_mle` |
| `logistic_only` | ported | `logistic_only` |
| `logitnorm.mle` | pending | `not in v0.3.0 computational-core port` |
| `loglogistic.mle` | ported | `loglogistic_mle` |
| `lognorm.mle` | ported | `lognormal_mle` |
| `logseries.mle` | ported | `logseries_mle` |
| `lomax.mle` | ported | `lomax_mle` |
| `lower_tri` | ported | `lower_tri_values` |
| `lower_tri.assign` | pending | `not in v0.3.0 computational-core port` |
| `mahala` | ported | `mahalanobis_sq` |
| `mat.mat` | intrinsic | `matmul intrinsic` |
| `mat.mult` | intrinsic | `matmul intrinsic / crossprod helpers` |
| `match.coefs` | pending | `not in v0.3.0 computational-core port` |
| `matrnorm` | pending | `not in v0.3.0 computational-core port` |
| `maxboltz.mle` | ported | `maxboltz_mle` |
| `mcnemar` | ported | `mcnemar_test` |
| `mcnemars` | pending | `not in v0.3.0 computational-core port` |
| `med` | ported | `median_r` |
| `mediandir` | ported | `spatial_median` |
| `min_max` | intrinsic | `minval / maxval` |
| `multinom.mle` | pending | `not in v0.3.0 computational-core port` |
| `multinom.nb` | ported | `fit_multinomial_nb` |
| `multinom.reg` | ported | `multinomial_regression` |
| `multinom.regs` | ported | `multinomial_regs` |
| `multinomnb.pred` | ported | `predict_multinomial_nb` |
| `multivmf.mle` | pending | `not in v0.3.0 computational-core port` |
| `mv.eeltest1` | ported | `mv_eeltest1` |
| `mv.eeltest2` | ported | `mv_eeltest2 (standard chi-square calibration)` |
| `mvbetas` | pending | `not in v0.3.0 computational-core port` |
| `mvkurtosis` | pending | `not in v0.3.0 computational-core port` |
| `mvlnorm.mle` | ported | `mvlognormal_mle` |
| `mvnorm.mle` | ported | `mvnorm_mle` |
| `mvt.mle` | ported | `mvt_mle` |
| `negative` | intrinsic | `PACK(x,x<0)` |
| `negbin.mle` | ported | `negbin_mle` |
| `normal.mle` | ported | `normal_mle` |
| `normlog.mle` | ported | `lognormal_mle` |
| `normlog.reg` | ported | `normlog_regression` |
| `normlog.regs` | ported | `normlog_regs` |
| `nth` | ported | `nth_value` |
| `odds` | pending | `not in v0.3.0 computational-core port` |
| `odds.ratio` | ported | `odds_ratio_2x2` |
| `omp` | ported | `omp_glm / omp_multivariate / omp_multinomial` (all upstream response-family branches) |
| `ompr` | ported | `ompr` |
| `ordinal.mle` | ported | `ordinal_mle` |
| `pareto.mle` | ported | `pareto_mle` |
| `pc.skel` | partial | `pc_skeleton` (Pearson/Spearman analytic modes; categorical and permutation modes pending) |
| `pdcor` | ported | `partial_distance_correlation` |
| `percent.ttest` | pending | `not in v0.3.0 computational-core port` |
| `percent.ttests` | pending | `not in v0.3.0 computational-core port` |
| `permcor` | pending | `not in v0.3.0 computational-core port` |
| `permutation` | pending | `not in v0.3.0 computational-core port` |
| `permutation.next` | ported | `next_permutation` |
| `permutation.prev` | ported | `prev_permutation` |
| `pois.test` | pending | `not in v0.3.0 computational-core port` |
| `poisdisp.test` | ported | `poisson_dispersion_test` |
| `poisson.anova` | pending | `not in v0.3.0 computational-core port` |
| `poisson.anovas` | pending | `not in v0.3.0 computational-core port` |
| `poisson.cat1` | pending | `not in v0.3.0 computational-core port` |
| `poisson.mle` | ported | `poisson_mle` |
| `poisson.nb` | ported | `fit_poisson_nb` |
| `poisson_only` | ported | `poisson_only` |
| `poissonnb.pred` | ported | `predict_poisson_nb` |
| `poly.cor` | pending | `not in v0.3.0 computational-core port` |
| `pooled.cov` | ported | `pooled_covariance` |
| `positive` | intrinsic | `PACK(x,x>0)` |
| `positive.negative` | intrinsic | `COUNT/PACK expressions` |
| `prop.reg` | ported | `proportion_regression` |
| `prop.regs` | ported | `proportion_regs` |
| `proptest` | ported | `proportion_test` |
| `proptests` | pending | `not in v0.3.0 computational-core port` |
| `qpois.reg` | ported | `quasipoisson_regression` |
| `qpois.regs` | ported | `quasipoisson_regs` |
| `quasi.poisson_only` | pending | `not in v0.3.0 computational-core port` |
| `quasipoisson.anova` | pending | `not in v0.3.0 computational-core port` |
| `quasipoisson.anovas` | pending | `not in v0.3.0 computational-core port` |
| `racg` | pending | `not in v0.3.0 computational-core port` |
| `rayleigh.mle` | ported | `rayleigh_mle` |
| `rbing` | pending | `not in v0.3.0 computational-core port` |
| `rbingham` | pending | `not in v0.3.0 computational-core port` |
| `read.directory` | r-only | `omitted: R object/package infrastructure` |
| `read.examples` | r-only | `omitted: R object/package infrastructure` |
| `regression` | pending | `not in v0.3.0 computational-core port` |
| `rel.risk` | ported | `relative_risk_2x2` |
| `rep_col` | ported | `rep_col` |
| `rep_row` | ported | `rep_row` |
| `rint.mle` | ported | `rint_mle` |
| `rint.reg` | ported | `rint_reg` |
| `rint.regbx` | ported | `rint_regbx` |
| `rint.regs` | ported | `rint_regs` |
| `rm.anova` | ported | `rm_anova` |
| `rm.anovas` | ported | `rm_anovas` |
| `rm.lines` | ported | `rm_lines` |
| `rmdp` | pending | `not in v0.3.0 computational-core port` |
| `rmvlaplace` | ported | `rmvlaplace` |
| `rmvnorm` | ported | `rmvnorm` |
| `rmvt` | ported | `rmvt` |
| `rowAll` | ported | `rowall` |
| `rowAny` | ported | `rowany` |
| `rowCountValues` | pending | `not in v0.3.0 computational-core port` |
| `rowFalse` | ported | `rowfalse` |
| `rowMads` | ported | `rowmads` |
| `rowMaxs` | ported | `rowmaxs` |
| `rowMedians` | ported | `rowmedians` |
| `rowMins` | ported | `rowmins` |
| `rowMinsMaxs` | pending | `not in v0.3.0 computational-core port` |
| `rowOrder` | ported | `roworder` |
| `rowRanks` | ported | `rowranks` |
| `rowShuffle` | pending | `not in v0.3.0 computational-core port` |
| `rowSort` | ported | `rowsort` |
| `rowTabulate` | pending | `not in v0.3.0 computational-core port` |
| `rowTrue` | ported | `rowtrue` |
| `rowTrueFalse` | ported | `rowtruefalse` |
| `rowVars` | ported | `rowvars` |
| `rowcvs` | ported | `rowcvs` |
| `rowhameans` | ported | `rowhameans` |
| `rowmeans` | ported | `rowmeans` |
| `rownth` | ported | `rownth` |
| `rowprods` | ported | `rowprods` |
| `rowrange` | ported | `rowranges` |
| `rows` | intrinsic | `array section x(idx,:)` |
| `rowsums` | ported | `rowsums` |
| `rvmf` | pending | `not in v0.3.0 computational-core port` |
| `rvonmises` | ported | `rvonmises` |
| `score.betaregs` | ported | `score_betaregs` |
| `score.expregs` | ported | `score_expregs` |
| `score.gammaregs` | ported | `score_gammaregs` |
| `score.geomregs` | ported | `score_geomregs` |
| `score.glms` | ported | `score_glms` |
| `score.invgaussregs` | ported | `score_invgaussregs` |
| `score.multinomregs` | ported | `score_multinomregs` |
| `score.negbinregs` | ported | `score_negbinregs` |
| `score.weibregs` | ported | `score_weibregs` |
| `score.ztpregs` | ported | `score_ztpregs` |
| `sftest` | pending | `not in v0.3.0 computational-core port` |
| `sftests` | pending | `not in v0.3.0 computational-core port` |
| `skew` | ported | `skewness_r` |
| `skew.test2` | pending | `not in v0.3.0 computational-core port` |
| `sort_cor_vectors` | pending | `not in v0.3.0 computational-core port` |
| `sort_mat` | pending | `not in v0.3.0 computational-core port` |
| `sort_unique` | pending | `not in v0.3.0 computational-core port` |
| `sort_unique.length` | pending | `not in v0.3.0 computational-core port` |
| `sourceR` | r-only | `omitted: R object/package infrastructure` |
| `sourceRd` | r-only | `omitted: R object/package infrastructure` |
| `spat.med` | ported | `spatial_median` |
| `spatmed.reg` | ported | `spatial_median_regression` |
| `spdinv` | ported | `spd_inverse` |
| `spml.mle` | pending | `not in v0.3.0 computational-core port` |
| `spml.reg` | pending | `not in v0.3.0 computational-core port` |
| `spml.regs` | pending | `not in v0.3.0 computational-core port` |
| `squareform` | ported | `squareform_from_vector` |
| `sscov` | pending | `not in v0.3.0 computational-core port` |
| `standardise` | ported | `standardise_cols / standardise_vector` |
| `submatrix` | intrinsic | `array section` |
| `tmle` | pending | `not in v0.3.0 computational-core port` |
| `tobit.mle` | ported | `tobit_mle` |
| `topological_sort` | ported | `topological_sort` |
| `total.dist` | ported | `total_dist` |
| `total.dista` | ported | `dista_matrix + sum` |
| `transpose` | intrinsic | `transpose intrinsic` |
| `ttest` | ported | `ttest1 / ttest2` |
| `ttest1` | ported | `ttest1` |
| `ttest2` | ported | `ttest2` |
| `ttests` | ported | `column_ttests` |
| `ttests.pairs` | pending | `not in v0.3.0 computational-core port` |
| `twoway.anova` | pending | `not in v0.3.0 computational-core port` |
| `twoway.anovas` | pending | `not in v0.3.0 computational-core port` |
| `ufactor` | r-only | `omitted: R object/package infrastructure` |
| `univglms` | ported | `univglms` (numeric-matrix families) |
| `univglms2` | pending | `not in v0.3.0 computational-core port` |
| `upper_tri` | ported | `upper_tri_values` |
| `upper_tri.assign` | pending | `not in v0.3.0 computational-core port` |
| `var2test` | pending | `not in v0.3.0 computational-core port` |
| `var2tests` | pending | `not in v0.3.0 computational-core port` |
| `varcomps.mle` | ported | `varcomps_mle` |
| `varcomps.mom` | ported | `varcomps_mom (use a one-column matrix for scalar input)` |
| `vartest` | ported | `vartest_chisq` |
| `vartests` | ported | `column_vartests` |
| `vecdist` | ported | `vecdist_upper` |
| `vm.mle` | ported | `vm_mle` |
| `vmf.mle` | ported | `vmf_mle` |
| `watson` | ported | `watson_test` |
| `weib.reg` | ported | `weibull_regression` |
| `weibull.mle` | ported | `weibull_mle` |
| `which.is` | intrinsic | `Fortran type/kind inquiry differs; no runtime R-type analogue` |
| `wigner.mle` | ported | `wigner_mle` |
| `wrapcauchy.mle` | pending | `not in v0.3.0 computational-core port` |
| `yule` | pending | `not in v0.3.0 computational-core port` |
| `zip.mle` | ported | `zip_mle` |
| `ztp.mle` | ported | `ztp_mle` |
| `.onAttach` | r-only | `omitted: R object/package infrastructure` |
