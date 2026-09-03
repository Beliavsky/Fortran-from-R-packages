# ape

Modern free-form Fortran translation of a substantial computational core of the
R package **ape** 5.8-1 (Analyses of Phylogenetics and Evolution).

The goal is a reusable numerical/phylogenetic library rather than an emulation
of R.  Plotting, interactive interfaces, S3 object plumbing, file-format UI,
Rcpp registration, and launchers for external programs are excluded.  See
`API_COVERAGE.md` for the exact translated/deferred surface, `NOTICE.md` for
license/provenance information, and `docs/BUILD_VERIFICATION.md` for the
compiler/FPM verification status of this archive.

## Dependency layout

This package is intended to live as the top-level directory `ape` inside
`Fortran-from-R-packages`, next to the shared packages `rfortran-core` and `rfortran-linalg`:

```text
Fortran-from-R-packages/
  ape/
  rfortran-core/
  rfortran-linalg/
```

`fpm.toml` therefore uses:

```toml
[dependencies]
rfortran-core = { path = "../rfortran-core" }
rfortran-linalg = { path = "../rfortran-linalg" }
```

No system BLAS/LAPACK link is used. `rfortran-linalg` supplies the eigenvalue,
SPD-solve, and Cholesky/inverse operations needed by PCoA and continuous ACE;
that shared package in turn pins the pure-Fortran `fortran-lapack` FPM
dependency. No dependency sources are vendored in `ape`.

## Build and test

From the extracted `ape` directory:

```text
fpm build
fpm test
fpm run --example nj_example
```

The maintained Fortran is free form and uses the shared `dp` real kind from
`rfortran-core` (`r_kinds`).  NaN-aware algorithms rely on IEEE behavior and must not be built with
`-ffast-math`, `-Ofast`, `-ffinite-math-only`, or equivalent assumptions.

## Main API

Import the umbrella module:

```fortran
use ape
```

Major groups are:

- `phylo_tree` construction, validation, parent/child helpers;
- `nj`, `bionj`, `mvr`, and FastME OLS/BME phylogeny reconstruction with NNI/SPR search;
- incomplete-distance `njs`, `bionjs`, `mvrs`, `triang_mtd`, and `triang_mtds`;
- additive and ultrametric missing-distance completion;
- node depths, paths, MRCA, patristic distances, branching times, balances,
  cherries, independent contrasts, Brownian PIC/ML/REML/GLS ancestral-state estimates,
  discrete ER/SYM/ARD ancestral-state likelihood fitting, `chronoMPL` dating,
  numeric branching-time assignment, monophyly and ultrametric checks;
- tree editing: tip pruning/keeping, clade extraction, singleton collapse,
  rerooting at an internal node, numeric outgroup rooting, unrooting, and
  deterministic multi/di conversion;
- PH85 topological and branch-score distances between trees;
- phylogenetic variance/covariance matrices;
- gamma statistic, Yule-rate fit, standard and extended (`bd.ext`) birth-death likelihood fitting,
  coalescent intervals, skyline estimation/AICc epsilon search, lineage-through-time
  coordinates, Moran's I, minimum spanning tree, quartet delta statistics,
  tree-count combinatorics, `dist.gene`, split compatibility, diversification
  GOF, `diversi.time` survival models, and sister-clade tests;
- split/clade collections, `prop_clades` support counting, bipartition counts,
  and rooted/unrooted strict or majority-rule consensus trees;
- PCoA with uncorrected negative-eigenvalue diagnostics and Lingoes or
  Cailliez correction, using the shared `rfortran-linalg` eigensolvers;
- PGLS/GLS covariance structures and fitting for Brownian, Martins, Grafen,
  Pagel, and Blomberg models; `compar.ou`, `compar.lynch`, `corphylo`, and the
  numerical PQL/REML core of `binaryPGLMM`;
- `chronopl` plus `chronos` clock/correlated/relaxed/discrete penalized-likelihood
  dating with fixed or interval calibrations;
- deterministic `reconstruct` ML, REML, GLS, GLS_ABM, GLS_OUS, and GLS_OU
  ancestral-reconstruction methods;
- `dist.dna`-style nucleotide distances, including the 17 upstream model names,
  with matrix-only dispatch for GG95 and BH87 and analytical variance for the
  ten upstream models that actually populate it;
- 17-state DNAbin base proportions, global deletion, ambiguity-aware
  segregating sites, terminal-gap conversion, contingency tables, pattern location,
  and ambiguity-aware translation under ape genetic codes 1--6.

All public numerical procedures report errors through an integer `info` or a
well-defined result rather than invoking R errors/warnings.

## Tree representation

`type(phylo_tree)` stores:

- `n_tip` tips numbered `1:n_tip`;
- internal nodes numbered above the tips;
- an `(nedge,2)` parent-child edge matrix;
- optional `edge_length`.

Unrooted trees returned by NJ/BIONJ/MVR and their incomplete-distance variants
use ape's rooted-storage convention with a trifurcating root. Most traversal
routines are independent of edge order.

## DNA encoding

The Fortran API uses these public constants:

```text
dna_unknown = 0
dna_a       = 1
dna_c       = 2
dna_g       = 3
dna_t       = 4
dna_gap     = 5
dna_r       = 6
dna_m       = 7
dna_w       = 8
dna_s       = 9
dna_k       = 10
dna_y       = 11
dna_v       = 12
dna_h       = 13
dna_d       = 14
dna_b       = 15
dna_n       = 16
```

The matrix orientation is `(taxon, site)`.  `dna_distance_matrix` defaults to
ape-style global deletion, except `INDEL` and `INDELBLOCK`, for which pairwise
gap handling is forced as in the R wrapper.  Set `pairwise_deletion=.true.` for
pairwise missing-data deletion in other models.

Supported names are `RAW`, `JC69`, `K80`, `F81`, `K81`, `F84`, `T92`, `TN93`,
`GG95`, `LOGDET`, `BH87`, `PARALIN`, `N`, `TS`, `TV`, `INDEL`, and
`INDELBLOCK`.  Optional gamma-shape correction is implemented for JC69, K80,
F81, and TN93, matching the upstream availability.

Analytical sampling variance is available through `dna_distance_with_variance`
and `dna_distance_matrix_with_variance` for JC69, K80, F81, K81, F84, T92,
TN93, GG95, LOGDET, and PARALIN. Upstream does not provide BH87 variance and
does not populate the allocated variance buffer for RAW/N/TS/TV/INDEL/
INDELBLOCK; the Fortran API reports status 6 for those models instead of
returning undefined memory. Corrected-distance models retain ape's known-base
deletion rules. Ambiguity-aware sequence utilities reproduce the historical
DNAbin bit semantics through the public language-neutral constants above.
`translate_dna` supports ape genetic codes 1--6 and preserves upstream results
for ambiguous codons.

## Example

The package currently ships 22 deterministic examples. In addition to the
reconstruction, DNA, consensus, diversification, ordination, skyline, and
birth-death examples from earlier parity passes, the current archive includes
`discrete_ace_example`, `pgls_example`, `chronopl_example`,
`chronos_clock_example`, `compar_ou_example`, `compar_lynch_example`,
`corphylo_example`, `binary_pglmm_example`, and `reconstruct_example`.

## License

GPL-2.0-only OR GPL-3.0-only, matching upstream ape.  See `COPYING`,
`LICENSE.GPL-3`, and `NOTICE.md`.
