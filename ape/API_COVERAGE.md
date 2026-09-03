# API coverage

This document records computational parity against upstream ape 5.8-1. It is
intentionally explicit: functionality not listed as translated should not be
inferred from the package name.

## Translated computational surface

### Tree representation and structural helpers

| Fortran API | Upstream concept |
|---|---|
| `phylo_tree`, `make_phylo_tree` | core `phylo` edge/length representation |
| `parent_vector`, `child_counts`, `edge_index_to_child` | tree traversal helpers |
| `is_binary_tree` | `is.binary` / binary topology check |
| `is_ultrametric_tree` | `is.ultrametric` |
| `mrca`, `mrca_many` | MRCA calculations |
| `node_path` | `nodepath` |
| `is_monophyletic`, `clade_tips`, `tip_descendant_matrix` | clade/descendant operations |
| `topological_distance_ph85` | `dist.topo(..., method="PH85")` split distance |
| `branch_score_distance` | `dist.topo(..., method="score")` branch-score distance |

### Tree editing

| Fortran API | Upstream concept / current scope |
|---|---|
| `has_singles` | `has.singles` |
| `is_rooted_tree` | `is.rooted` structural root-degree convention |
| `collapse_singles` | `collapse.singles`; current tree type has no separate `root.edge` field |
| `drop_tips` | numeric-index computational core of `drop.tip`, returning old-tip mapping |
| `keep_tips` | numeric-index computational core of `keep.tip`, returning old-tip mapping |
| `extract_clade` | `extract.clade` for a specified node number |
| `reroot_node` | `root(..., node=...)` at an existing internal node |
| `root_outgroup` | numeric-tip outgroup rooting, including `resolve.root=TRUE` zero-length basal insertion |
| `unroot_tree` | standard rooted-binary to unrooted conversion, preserving path lengths |
| `di2multi` | `di2multi(..., tip2root=FALSE)` edge contraction |
| `multi2di` | deterministic `multi2di(..., random=FALSE)` zero-edge resolution |

Label/S3 plumbing, arbitrary edge-position rooting, `root.edge` special cases, and
`multiPhylo` wrappers are not represented by these typed routines.

### Tree metrics and comparative calculations

| Fortran API | Upstream concept |
|---|---|
| `node_depth_edgelength` | `node.depth.edgelength` |
| `node_depth_count` | `node.depth` methods 1 and 2 |
| `dist_nodes` | `dist.nodes` / patristic node distances |
| `descendant_tip_counts` | descendant-tip counts |
| `balance_counts` | binary-tree balance counts |
| `cherry_count` | cherry count |
| `branching_times` | `branching.times` |
| `phylogenetic_vcv` | `vcv.phylo`, covariance or correlation form |
| `pic` | `pic`, including optional scaled contrasts, rescaled tree, and ancestral values |
| `ace_pic` | `ace(..., type="continuous", model="BM", method="pic")` ancestral estimates/contrast variances |
| `ace_continuous_ml` | continuous Brownian `ace(..., method="ML")`, including sigma2 and Hessian-based SEs |
| `ace_continuous_reml` | continuous Brownian `ace(..., method="REML")`, including PIC-root variance fit and residual likelihood |
| `ace_continuous_gls` | continuous `ace(..., method="GLS")` under Brownian, Martins, or Grafen correlation structures |
| `chrono_mpl` | `chronoMPL`, including dated edges, internal-node SEs, and clock-test p-values |
| `compute_brtime` | deterministic numeric-age branch of `compute.brtime(..., force.positive=FALSE)` |

### Discrete ancestral states, phylogenetic GLS, and comparative models

| Fortran API | Upstream concept / current scope |
|---|---|
| `ace_rate_index_matrix` | discrete `ace` ER/SYM/ARD transition-index construction |
| `ace_discrete_likelihood` | discrete-state continuous-time Markov pruning plus joint/marginal ancestral likelihoods |
| `ace_discrete_fit` | fitted ER/SYM/ARD or custom indexed discrete `ace`, including rate SEs |
| `cor_brownian`, `cor_martins`, `cor_grafen`, `cor_pagel`, `cor_blomberg` | ape phylogenetic correlation structures used by PGLS/nlme workflows |
| `pgls_fit` | fixed-covariance ML/REML generalized least squares |
| `pgls_fit_model` | Brownian, Martins, Grafen, Pagel, or Blomberg GLS with optional one-parameter profiling |
| `compar_ou_fit`, `compar_ou_likelihood` | numerical core of `compar.ou`, including regime shifts and profiled/fixed alpha |
| `compar_lynch_fit` | numerical multivariate Lynch comparative method |
| `corphylo_fit`, `corphylo_objective` | multivariate `corphylo` likelihood/fit with optional measurement error and covariates |
| `binary_pglmm_fit`, `binary_pglmm_reml_objective` | binary phylogenetic GLMM PQL/REML numerical core |
| `reconstruct_fit` | deterministic `reconstruct` ML, REML, GLS, GLS_ABM, GLS_OUS, and GLS_OU methods |
| `reconstruct_gls_bm`, `reconstruct_gls_abm`, `reconstruct_gls_ou_stationary`, `reconstruct_gls_ou` | direct ancestral-reconstruction entry points |

R formula/model-frame construction, `nlme` object methods, printed summaries, and
stochastic optimizer/simulation wrappers are intentionally replaced by typed
arrays and deterministic numerical optimization.

### Penalized-likelihood molecular dating

| Fortran API | Upstream concept / current scope |
|---|---|
| `chronopl_fit`, `chronopl_objective` | `chronopl` penalized likelihood with fixed/interval node calibrations |
| `chronos_clock_fit` | `chronos(..., model="clock")` |
| `chronos_fit` | `chronos` correlated, relaxed, discrete, and clock rate models |
| `chronos_objective` | fixed-parameter Poisson likelihood plus model-specific rate penalty |

The dating routines use deterministic bounded optimization, accept fixed or
interval calibration ages, and return dated trees, fitted rates, likelihoods,
and PHIIC summaries. The R `chronopl` leave-one-tip-out cross-validation wrapper
and plotting/object-printing layers are not translated.

### Phylogeny reconstruction and incomplete distances

| Fortran API | Upstream source/function |
|---|---|
| `nj` | `nj`, `src/nj.c` |
| `bionj` | `bionj`, `src/BIONJ.c` |
| `mvr` | `mvr`, `src/mvr.c` |
| `njs` | `njs`, `src/njs.c`, including `fs` shortlist/tie-breaking behavior |
| `bionjs` | `bionjs`, `src/bionjs.c` |
| `mvrs` | `mvrs`, `src/mvrs.c` |
| `triang_mtd` | `triangMtd`, `src/triangMtd.c` |
| `triang_mtds` | `triangMtds`, `src/triangMtds.c` |
| `additive_completion` | `src/additive.c` |
| `ultrametric_completion` | `src/ultrametric.c` |
| `fastme_ols` | `fastme.ols`, greedy OLS minimum evolution with optional OLS-NNI |
| `fastme_bal` | `fastme.bal`, greedy balanced minimum evolution with optional bNNI and SPR |

`fastme_ols` and `fastme_bal` use a modern explicit unrooted edge graph rather
than ape/FastME's pointer-heavy C structures, while evaluating the same OLS/BME
branch-length and search objectives. During development all five search modes
were regression-compared with ape's native kernels on 750 generated complete
distance matrices for 6--8 taxa with no topology mismatches.

The starred/sparse methods accept NaN or negative symmetric off-diagonal
entries as missing. The NJ* family preserves ape's transition from incomplete
NJ* selection to ordinary NJ selection once the reduced matrix becomes
complete. The reconstruction APIs accept full symmetric matrices rather than R
`dist` objects and return `phylo_tree` directly.

### Split collections, clade support, and consensus

| Fortran API | Upstream concept / current scope |
|---|---|
| `split_collection`, `prop_part` | rooted clade collection/frequencies from `prop.part` |
| `bitsplits`, `tree_bipartitions` | `bitsplits` and one-wise canonical edge bipartitions |
| `count_bipartitions` | `countBipartitions` |
| `prop_clades` | `prop.clades`, including SHORTwise unrooted treatment |
| `consensus_tree` | strict and majority-rule `consensus`, `p` in [0.5,1], rooted or unrooted treatment |

The Fortran interfaces assume identical numeric tip ordering instead of R label
matching/compression. `consensus_tree` returns node-support proportions as a
separate array because `phylo_tree` intentionally has no R `node.label` field.
Bootstrap resampling itself remains deferred because it requires an RNG/callback
policy.

### Phylogenetic/statistical utilities

| Fortran API | Upstream concept |
|---|---|
| `gamma_stat` | Pybus-Harvey `gammaStat` |
| `yule_fit` | standard `yule` pure-birth rate/log-likelihood fit |
| `coalescent_intervals` | lineage/interval extraction used by coalescent summaries |
| `collapsed_intervals` | `collapsed.intervals` deterministic interval grouping |
| `skyline_from_intervals`, `skyline_tree` | skyline population-size estimates, likelihood, and AICc |
| `find_skyline_epsilon` | deterministic `find.skyline.epsilon` grid search |
| `birthdeath_fit`, `birthdeath_from_times` | standard `birthdeath` likelihood fit, Hessian SEs, and ape fixed-profile intervals |
| `birthdeath_deviance` | exact standard birth-death deviance for fixed `(d/b, b-d)` |
| `birthdeath_extended_fit`, `birthdeath_extended_from_data` | `bd.ext` conditional/unconditional taxonomic-richness likelihood fits |
| `birthdeath_extended_deviance` | exact `bd.ext` deviance for fixed `(d/b, b-d)` and richness data |
| `ltt_coordinates` | numerical core of `ltt.plot.coords`, without plotting |
| `moran_i` | `Moran.I` observed, expectation, and SD |
| `minimum_spanning_tree` | `mst` computational kernel |
| `delta_plot_statistics` | quartet delta statistic/binning from `delta.plot` kernel |
| `howmanytrees` | `howmanytrees`, including overflow-safe scientific representation |
| `gene_distance_matrix` | `dist.gene` count/percentage distances and variance |
| `splits_compatible`, `all_splits_compatible` | split compatibility tests (`arecompatible`/bitsplit logic) |
| `diversification_gof` | `diversi.gof` Cramer-von Mises and Anderson-Darling statistics |
| `diversification_time` | complete deterministic three-model `diversi.time` survival analysis |
| `slowinski_guyer_test` | Slowinski-Guyer sister-clade test |
| `mcconway_sims_test` | McConway-Sims sister-clade test |
| `diversity_contrasts` | deterministic contrast calculation used by `diversity.contrast.test` |

The standard `birthdeath` and extended `bd.ext` likelihood/deviance equations are
translated directly. R's `nlm` driver is replaced by a deterministic bounded
profile/golden-section search over `d/b` and positive `b-d`; local standard
errors use the numerical Hessian of the same deviance. The standard model also
reproduces ape's fixed-other-parameter profile-interval stepping rule.

`diversity_contrasts` intentionally omits permutation/Wilcoxon inference and
RNG-driven wrappers. `chi_square_survival` is exposed as a supporting numerical
helper used by the translated sister-clade tests.

### Ordination and decomposition

| Fortran API | Upstream concept / current scope |
|---|---|
| `pcoa` | `pcoa`, including ordinary eigen diagnostics and Lingoes/Cailliez corrections |
| `pcoa_result` | eigenvalues, relative/cumulative summaries, broken-stick expectations, and coordinates |

PCoA delegates symmetric and general eigenvalue calculations to the sibling
`rfortran-linalg` package instead of carrying private LAPACK wrappers. Row names,
printing, and `biplot.pcoa` plotting are intentionally omitted.

### DNA and sequence computations

`dna_distance_matrix` exposes all model names accepted by upstream `dist.dna`:

- RAW
- JC69
- K80
- F81
- K81
- F84
- T92
- TN93
- GG95
- LOGDET
- BH87 (directional matrix)
- PARALIN
- N
- TS
- TV
- INDEL
- INDELBLOCK

`dna_distance` provides pair-level models that do not require a whole alignment
for model estimation/directional output. GG95 and BH87 are handled by the
matrix API. Gamma correction is translated for JC69, K80, F81, and TN93.

`dna_distance_with_variance` and `dna_distance_matrix_with_variance` expose the
analytical sampling variances actually computed by upstream `dist_dna.c` for
JC69, K80, F81, K81, F84, T92, TN93, GG95, LOGDET, and PARALIN. BH87 variance
is unavailable upstream, and upstream allocates but does not populate variance
for RAW/N/TS/TV/INDEL/INDELBLOCK; the Fortran API returns status 6 for those
cases instead of exposing uninitialized data. The LogDet/ParaLin implementation
reproduces the upstream `dgesv` overwrite semantics with an internal 4-by-4
partial-pivot LU solve and therefore adds no BLAS/LAPACK dependency.

Additional sequence APIs: `dna_base_frequencies`, `dna_base_proportions`
(all 17 DNAbin states), `dna_leading_trailing_gaps_to_n`,
`dna_global_deletion_mask`, ambiguity-aware `dna_segregating_sites`,
`dna_contingency_table`, `dna_pattern_positions`, and `translate_dna`.
The public encoding now represents A/C/G/T plus R/M/W/S/K/Y/V/H/D/B/N, gap,
and unknown. `translate_dna` reproduces ape's historical DNAbin lookup behavior
for all six genetic codes, including ambiguous codons.

## Deliberately deferred computational areas

These are meaningful upstream computations, but are not claimed as translated
in this archive:

- decomposition-heavy ordinations other than the translated PCoA workflow;
- broader multi-tree operations not covered by the translated split collection,
  clade-support, and consensus APIs, especially bootstrap resampling/callback
  workflows;
- arbitrary edge-position and label-aware `root.phylo` behavior beyond numeric
  outgroup rooting and rerooting at an existing internal node;
- stochastic simulation, random trees/traits, bootstrap/permutation wrappers,
  stochastic SANN/MCMC workflows, and other routines requiring a package-level
  RNG policy;
- R formula/model-frame, `nlme` object, and S3 plumbing around the translated
  GLS/PGLS/comparative-model numerical cores;
- `chronopl` cross-validation and R-side `chronos` convenience/object methods
  beyond the translated clock/correlated/relaxed/discrete numerical fits;
- remaining specialized comparative methods not covered by `ace`, PGLS,
  `compar.ou`, `compar.lynch`, `corphylo`, `binaryPGLMM`, and `reconstruct`;
- remaining DNAbin/list-oriented sequence operations not covered by the
  translated 17-state ambiguity semantics, six-code translation, and distance utilities;
- sequence/tree parsers and writers, XML/BioConductor bridges, compression,
  and external executable orchestration.

Matrix-exponential functionality is not copied from ape because a top-level
`expm` translation already exists in `Fortran-from-R-packages`; an API that
needs it should use that sibling package rather than vendoring code. Likewise,
future linear algebra should use `rfortran-linalg`, and iterative spectral work
should use `rfortran-arpack`.

## Intentional interface differences

- R `dist`, `phylo`, `DNAbin`, lists, attributes, factors, environments, and S3
  methods are replaced by explicit typed Fortran arguments/results.
- Functions do not plot, prompt, start viewers, invoke R, or launch external
  phylogenetics applications.
- Error handling uses `info` status codes and IEEE NaNs where the mathematical
  result is undefined.
