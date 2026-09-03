# Validation notes

The deterministic test program exercises tree construction and validation,
MRCA/path/depth calculations, phylogenetic covariance, independent contrasts,
PH85 and branch-score tree distances, tree editing, split/clade collections,
consensus construction, complete/incomplete/FastME reconstruction,
missing-distance completion, branching times, `chronoMPL`, numeric branch-time
assignment, PCoA, continuous ACE ML/REML/GLS, discrete ACE, PGLS, `chronopl`/`chronos`,
OU/Lynch/corphylo/binary-PGLMM comparative fits, ancestral reconstruction,
lineage-through-time coordinates, skyline and standard/extended birth-death likelihoods, gamma/Yule/coalescent summaries, Moran's I, MST, quartet delta statistics,
tree-count combinatorics, discrete genetic distances, split compatibility,
diversification tests, and representative DNA models/utilities.

Key numerical invariants used in the tests include:

- NJ and BIONJ exactly reproduce the classic five-taxon additive matrix in
  patristic distances.
- MVR with constant pairwise variances reduces to NJ and reproduces the same
  additive matrix.
- NJ* reproduces ordinary NJ on complete input and recovers a removed entry of
  the classic additive matrix.
- BIONJ* and MVR* are checked against branch lengths emitted by ape's original
  C kernels for a deterministic non-additive six-taxon matrix with two missing
  distances. The Fortran values agree to floating-point precision.
- `triang_mtd` and `triang_mtds` are checked against ape's original
  `triangMtd.c`/`triangMtds.c` behavior and reproduce the additive reference
  tree with complete and incomplete input.
- FastME OLS insertion, OLS+NNI, BME insertion, BME+bNNI, and BME+bNNI+SPR
  have deterministic package regressions. During development those five modes
  also matched ape's native FastME C topology on 750 generated matrices with
  6--8 taxa; reconstructed tip-distance differences stayed near roundoff.
- tip pruning, keeping, clade extraction, singleton collapse, rerooting,
  single- and multi-tip outgroup rooting, unrooting, and deterministic
  `multi2di`/`di2multi` preserve the expected patristic distances.
- split collection tests check one-wise bipartitions, rooted `prop.part`
  frequencies, rooted/unrooted `prop_clades`, and strict/majority consensus. A
  three-tree quartet sample retains the split occurring in 2/3 of trees at
  `p=0.5`, while strict consensus collapses it to a star.
- alternate four-tip unrooted topologies have PH85 distance 2 and branch-score
  distance `sqrt(5)` for internal-edge lengths 1 and 2.
- `howmanytrees` covers rooted/unrooted, binary/nonbinary, and labeled/unlabeled
  reference counts; `dist.gene` tests count/percentage forms and variances.
- split compatibility and the deterministic diversification/Slowinski-Guyer/
  McConway-Sims calculations are checked against closed-form reference values.
- `diversification_time` is regression-tested for all three upstream models,
  including the Weibull shape solve, standard rate estimates, likelihoods, and
  both one-degree-of-freedom likelihood-ratio p-values.
- lineage-through-time coordinates are checked in backward/forward and upper/
  lower step conventions on a balanced ultrametric tree.
- a balanced four-tip unit-branch tree has the expected covariance matrix,
  independent contrasts, PIC ancestral-state estimates, branching times, gamma
  statistic, and Yule estimates. A non-clocklike four-tip tree checks
  `chronoMPL` ages, standard errors, and molecular-clock p-values; numeric
  `compute_brtime` is checked edge-by-edge.
- PCoA regression values cover the ordinary negative eigenvalue, Lingoes
  constant/eigenvalues, and Cailliez correction constant/eigenvalues on a
  non-Euclidean four-object distance matrix.
- continuous Brownian `ace` ML and REML are checked on a balanced four-tip tree
  for ancestral states, sigma-squared estimates, Hessian standard errors, and
  likelihoods; continuous GLS is checked under phylogenetic covariance models.
- discrete `ace` checks ER/SYM/ARD indexing, fixed-rate pruning/smoothing, fitted
  ER rate/SE/log likelihood, and ancestral state probabilities.
- PGLS checks Brownian, Martins, Grafen, Pagel, and Blomberg matrices, fixed GLS,
  and profiled Pagel lambda. `chronopl` and `chronos` check fixed objectives,
  calibrated ages, equal-rate clock/correlated fits, relaxed/discrete objectives,
  and PHIIC summaries.
- `compar.ou`, `compar.lynch`, `corphylo`, and `binaryPGLMM` have deterministic
  parameter/objective regressions, and `reconstruct` checks its ML/REML/GLS/OU
  ancestral-reconstruction paths.
- skyline tests cover deterministic interval collapse, modern and old-style
  population-size estimates, likelihood, and AICc. Standard `birthdeath` tests
  cover both an interior optimum and the zero-extinction boundary; `bd.ext`
  checks conditional/unconditional fixed deviances, an interior conditional
  optimum, an unconditional boundary fit, and the tree-to-data entry path.
- DNA tests cover simple closed-form RAW/JC69/K80/F81/K81/F84/T92/TN93 values,
  BH87, GG95, LOGDET/PARALIN domains, TS/TV counts, gap preservation for INDEL,
  base frequencies, 17-state base proportions, deletion masks, ambiguity-aware
  segregating sites, terminal-gap conversion, pattern search, and translation
  under ape genetic codes 1--6. Analytical-variance regression values are checked for JC69, K80,
  F81, K81, F84, T92, TN93, GG95, LOGDET, and PARALIN, including gamma cases.
  The LogDet/ParaLin constants retain upstream `dgesv`-overwrite semantics.

During development, the newly translated NJ*/BIONJ*/MVR*, triangle-method,
and FastME routines were compared with small standalone builds of the
corresponding upstream C source. DNAbin `seg.sites` was exhaustively compared
for all 4,913 three-state columns, and translation was exhaustively compared for
29,478 state/code combinations. Those validation harnesses are not distributed
in the package.

The source is audited before packaging for duplicate Fortran files, disallowed
statement semicolons, self-comparison NaN idioms, legacy real-kind forms,
missing dummy `INTENT`/FORD annotations, vendored dependency source, and build
products.
