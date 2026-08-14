# sna-fortran

Modern Fortran translation of the computational portions of the R package
`sna` (Tools for Social Network Analysis), packaged for the Fortran Package
Manager (FPM).

Upstream package: `sna` 2.8 by Carter T. Butts.  The upstream package is
licensed under GPL version 2 or later.  This translation preserves that
license; see `COPYING` and `upstream/DESCRIPTION`.

This is a source-level numerical port, not an R compatibility layer.  Graphs
are represented directly as `real(dp)` adjacency matrices and graph stacks as
rank-3 arrays.  R S3 classes, printing, plotting, interactive graphics, and
`network` object conversion are intentionally not reproduced.

## Build

```text
fpm build
fpm test
fpm run --example basic_network
```

The sources use Fortran 2018 syntax and no external numerical library.
`dp` is `real64` from `iso_fortran_env`.  IEEE NaN is used where the R code
uses numeric `NA`.

## Modules

- `sna_kinds` - kinds, NaN/infinity helpers.
- `sna_types` - result types for distances, components, censuses, regressions,
  permutation tests, and brokerage.
- `sna_linalg` - small self-contained linear algebra used by the statistical
  routines.
- `sna_prep` - graph preparation, symmetrization, dichotomization,
  vectorization, interval graphs, CSS conversion, and log-domain helpers.
- `sna_graph` - connectivity, distances, components, cut points, cores,
  flow, path/cycle/clique and dyad/triad census, and graph-level indices.
- `sna_centrality` - node centrality and prestige measures.
- `sna_random` - random graph generators and rewiring routines.
- `sna_permutation` - constrained permutations and label optimization.
- `sna_multivariate` - covariance/correlation, Hamming and structural
  distances, central graphs, and distance matrices.
- `sna_roles` - structural/regular equivalence and blockmodel operations.
- `sna_bn_triad` - exact Skvoretz biased-net triadic pseudolikelihood kernel translated from the upstream C source.
- `sna_models` - Bayesian network accuracy, biased-net models, brokerage,
  consensus, network autocorrelation, network regression/logit, L-NAM, and
  P-star-style logistic change-statistic fitting.
- `sna_testing` - CUG and QAP Monte Carlo tests.
- `sna` - umbrella module that re-exports the public API.

## Main translated routines

R names are converted to conventional Fortran snake_case names where needed.
The major mappings include:

| R `sna` routine | Fortran routine |
| --- | --- |
| `add.isolates` | `add_isolates` |
| `bbnam.probtie` | `bbnam_probtie` |
| `bbnam.fixed` | `bbnam_fixed_posterior`, `bbnam_fixed_draws` |
| `bbnam.pooled` | `bbnam_pooled` |
| `bbnam.actor` | `bbnam_actor` |
| `bbnam.jntlik` | `bbnam_joint_loglik` |
| `bbnam.bf` | `bbnam_bayes_factor` |
| `betweenness` | `betweenness` |
| `bicomponent.dist` | `bicomponent_dist` |
| `blockmodel` | `blockmodel` |
| `blockmodel.expand` | `blockmodel_expand_density` |
| `bn` | `bn_fit` |
| `bn.nlpl.dyad` | `bn_nll_dyad` |
| `bn.nlpl.edge` | `bn_nll_edge` |
| `bn` triadic pseudolikelihood kernel | `bn_lpt_triad`, `bn_nll_triad` |
| biased-net triad probability kernel | `bn_ptriad` |
| `bonpow` | `bonpow` |
| `brokerage` | `brokerage` |
| `centralgraph` | `centralgraph` |
| `centralization` core | `centralization_from_scores` |
| `clique.census` | `clique_census` |
| `closeness` | `closeness` |
| `component.dist` | `component_dist` |
| `component.largest` | `component_largest_mask` |
| `component.size.byvertex` | `component_size_byvertex` |
| `components` | `components` |
| `connectedness` | `connectedness` |
| `consensus` | `consensus` |
| `cug.test` / `cugtest` | `cug_test` |
| `cutpoints` | `cutpoints` |
| `degree` | `degree` |
| `diag.remove` | `diag_remove` |
| `dyad.census` | `dyad_census` |
| `efficiency` | `efficiency` |
| `ego.extract` numerical selection | `ego_extract_mask` |
| `eval.edgeperturbation` scalar form | `eval_edgeperturbation` |
| `evcent` | `evcent` |
| `event2dichot` | `event2dichot` |
| `flowbet` | `flowbet` |
| `gcor` | `graph_correlation` |
| `gcov` | `graph_covariance` |
| `gden` | `gden` |
| `geodist` | `geodist` |
| `gilschmidt` | `gilschmidt` |
| `grecip` | `grecip` |
| `gscor` | `gscor` |
| `gscov` | `gscov` |
| `gt` | `graph_transpose` |
| `gtrans` | `gtrans` |
| `gvectorize` | `gvectorize` |
| `hdist` | `hdist` |
| `hierarchy` | `hierarchy` |
| `infocent` | `infocent` |
| `interval.graph` | `interval_graph` |
| `is.connected` | `is_connected` |
| `is.isolate` | `is_isolate` |
| `isolates` | `isolates` |
| `kcores` | `kcores` |
| `kcycle.census` | `kcycle_census` |
| `kpath.census` | `kpath_census` |
| `lab.optimize.anneal` | `lab_optimize_anneal` |
| `lab.optimize.exhaustive` | `lab_optimize_exhaustive` |
| `lab.optimize.gumbel` | `lab_optimize_gumbel` |
| `lab.optimize.hillclimb` | `lab_optimize_hillclimb` |
| `lab.optimize.mc` | `lab_optimize_mc` |
| `lnam` | `lnam` |
| `loadcent` | `loadcent` |
| `logMean` | `log_mean` |
| `logSub` | `log_sub` |
| `logSum` | `log_sum` |
| `lower.tri.remove` | `lower_tri_remove` |
| `lubness` | `lubness` |
| `make.stochastic` | `make_stochastic` |
| `maxflow` | `maxflow` |
| `mutuality` | `mutuality` |
| `nacf` | `nacf` |
| `netcancor` | `netcancor` |
| `neighborhood` | `neighborhood` |
| `netlm` | `netlm` |
| `netlogit` | `netlogit` |
| `npostpred` scalar form | `npostpred_scalar` |
| `nties` | `nties` |
| `numperm` | `numperm` |
| `prestige` | `prestige` |
| `pstar` | `pstar` |
| `pstar` scalar-effects convenience form | `pstar_basic` |
| `qaptest` | `qap_test` |
| `reachability` | `reachability` |
| `redist` | `redist` |
| `rewire.ud` | `rewire_ud` |
| `rewire.ws` | `rewire_ws` |
| `rgbn` MCMC kernel | `rgbn_mcmc` |
| `rgnm` | `rgnm` |
| `rgnmix` | `rgnmix_probability`, `rgnmix_exact` |
| `rgraph` | `rgraph` |
| `rguman` | `rguman` |
| `rgws` | `rgws` |
| `rmperm` | `rmperm` |
| `rperm` | `rperm` |
| `sdmat` | `sdmat` |
| `sedist` | `sedist` |
| `simmelian` | `simmelian` |
| `sr2css` | `sr2css` |
| `stackcount` | `stackcount` |
| `stresscent` | `stresscent` |
| `structdist` | `structdist` |
| `structure.statistics` | `structure_statistics` |
| `symmetrize` | `symmetrize` |
| `triad.census` | `triad_census` |
| `triad.classify` | `triad_classify` |
| `upper.tri.remove` | `upper_tri_remove` |

## Deliberately omitted R-only/presentation functionality

The following are not part of the computational Fortran library:

- all `gplot*`, `plot.*`, layout/graphics and `visualization.R` functionality;
- `print.*`, `summary.*`, `coef.*`, and R class construction;
- `as.sociomatrix.sna`, `as.edgelist.sna`, and `network` object dispatch;
- `read.dot`, `read.nos`, `write.dl`, and `write.nos` file-format wrappers;
- R formula, `...`, and arbitrary `match.fun` dispatch wrappers such as the
  fully generic forms of `gapply`/`gliop`.

The numerical operations behind those wrappers are exposed directly where
practical.

## Known API/algorithm adaptations

This Fortran port preserves the computational intent while replacing R runtime
services with explicit typed APIs.  The main adaptations are:

1. `make_stochastic(..., mode='rowcol')` uses deterministic iterative
   proportional fitting (Sinkhorn scaling) rather than the upstream stochastic
   annealing procedure.  It targets the same row/column normalization fixed
   point when one is feasible.
2. `bicomponent_dist` exposes per-vertex bicomponent membership counts rather
   than reproducing every component of the R list object.
3. Path/cycle census results use explicit Fortran result arrays.  Core counts
   and vertex/dyad participation data are present, while R-specific list layouts
   and every optional copath-comembership object are not duplicated.
4. `bn_fit` includes the exact specialized Skvoretz triadic pseudolikelihood
   kernel from the upstream C implementation and defaults to `mple.triad`, but
   uses a self-contained bounded coordinate/pattern search instead of R's
   `optim`/BFGS driver.  `mple.dyad`, `mple.edge`, and `mtle` are also available.
5. `netlm` and `netlogit` expose classical and permutation-based network fits
   directly.  R aliases whose only purpose is dynamic callback construction are
   consolidated into the corresponding explicit Fortran paths.
6. `lnam` contains the likelihood, conditional beta fit, autocorrelation
   parameter optimization, fitted values, residuals, disturbances, and the
   upstream-style numerical-Hessian asymptotic covariance/standard errors.
   Optimization uses the self-contained coordinate search rather than R's BFGS.
7. `pstar` supports scalar graph effects, node-valued outdegree/indegree/
   betweenness/closeness effects, centralization effects, continuous attribute
   differences, and categorical membership-similarity columns.  R formula/S3
   construction is intentionally replaced by explicit effect-name and array
   arguments.
8. `brokerage` includes Gould-Fernandez raw role counts, totals, null
   expectations, standard deviations, and z-scores, returned in a typed result
   instead of an R summary object.
9. `bbnam_bayes_factor` implements the Monte Carlo model-integration wrapper;
   its return value is a typed result rather than the upstream R object.
10. `netcancor` uses self-contained canonical-correlation linear algebra and
    implements QAP, CUG, density-conditioned CUG, and tie-resampling Monte Carlo
    nulls without relying on R's `cancor`.
11. `equiv.clust` is not reproduced because its final product is an R `hclust`
    object.  The computational dissimilarities `sedist` and `redist` are
    translated.
12. Arbitrary R callbacks (`match.fun`, `...`) and `network`/S3 coercion are not
    emulated.  Computational operations are exposed as typed procedures and,
    where useful, Fortran procedure arguments.

These adaptations are explicit so that the package does not claim R object or
dynamic-dispatch compatibility that a standalone Fortran library cannot have.

## Validation

`test/test_core.f90` exercises representative graph, census, centrality,
multivariate, random-graph, biased-network, Bayesian-network-accuracy,
brokerage, and regression routines.  The translated 64-state biased-net
triadic pseudolikelihood kernel was additionally compared directly with the
upstream C implementation for all 64 edge states at a nontrivial parameter
point; all double-precision values matched exactly.  The translation has also been compiled
from a clean module directory with GNU Fortran using:

```text
gfortran -std=f2018 -Werror=implicit-interface
```

All source lines are within the standard 132-column free-form limit, so no
compiler-specific unlimited-line-length flag is required.  The package contains
no external binary dependencies.
