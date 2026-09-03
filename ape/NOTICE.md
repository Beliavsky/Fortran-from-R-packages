# NOTICE and provenance

This directory is an independent modern-Fortran translation of computational
algorithms from the R package **ape** (Analyses of Phylogenetics and Evolution),
version **5.8-1**, dated 2024-12-10.

Upstream project: https://github.com/emmanuelparadis/ape

The upstream `DESCRIPTION` and `inst/CITATION` files are preserved verbatim as
`UPSTREAM_DESCRIPTION` and `UPSTREAM_CITATION`.  The upstream GNU GPL version 2
text is preserved verbatim as `COPYING`; the GNU GPL version 3 text is included
as `LICENSE.GPL-3`.  The translated source is distributed under the same choice
of **GPL-2.0-only OR GPL-3.0-only**.

## Upstream attribution retained for translated native kernels

- FastME OLS/BME insertion and topology search: ape `src/me.c`, `src/me_ols.c`,
  `src/me_balanced.c`, `src/NNI.c`, `src/bNNI.c`, and `src/SPR.c`; Copyright
  2007-2008 Olivier Gascuel, Rick Desper, Vincent Lefort, with modifications by
  Emmanuel Paradis; SPR code Copyright 2009 Richard Desper.
- `nj`: ape `src/nj.c`, Copyright 2006-2023 Emmanuel Paradis.
- `bionj`: ape `src/BIONJ.c`, Copyright 2007-2008 Olivier Gascuel and
  Hoa Sien Cuong; R port by Vincent Lefort and Emmanuel Paradis.
- `mvr`: ape `src/mvr.c`, Copyright 2011-2012 Andrei-Alin Popescu.
- `njs`: ape `src/njs.c`, Copyright 2011-2013 Andrei-Alin Popescu.
- `bionjs`: ape `src/bionjs.c`, Copyright 2011-2014 Andrei-Alin Popescu.
- `mvrs`: ape `src/mvrs.c`, Copyright 2011-2012 Andrei-Alin Popescu.
- `triang_mtd` and `triang_mtds`: ape `src/triangMtd.c` and
  `src/triangMtds.c`, Copyright 2011-2012 Andrei-Alin Popescu.
- `additive_completion` and `ultrametric_completion`: ape `src/additive.c` and
  `src/ultrametric.c`, Copyright 2011 Andrei-Alin Popescu.
- continuous and discrete `ace`: ape `src/pic.c`, Copyright 2006-2017 Emmanuel Paradis;
  `R/ace.R`, Copyright 2005-2024 Emmanuel Paradis and 2005 Ben Bolker.
- phylogenetic correlation/PGLS numerical models: ape `R/PGLS.R`, Copyright
  2004-2021 Julien Dutheil and 2006-2017 Emmanuel Paradis.
- penalized molecular dating: ape `R/chronopl.R` and the corresponding
  `chronos` numerical workflow, Copyright Emmanuel Paradis as recorded upstream.
- OU, Lynch, multivariate correlation, binary phylogenetic GLMM, and ancestral
  reconstruction workflows: ape `R/compar.ou.R`, `R/compar.lynch.R`,
  `R/corphylo.R`, `R/binaryPGLMM.R`, and `R/reconstruct.R`; complete upstream
  author/copyright attribution remains in `UPSTREAM_DESCRIPTION` and source headers.
- molecular dating helpers: `R/chronoMPL.R`, Copyright 2007-2017 Emmanuel
  Paradis; `R/compute.brtime.R`, Copyright 2011-2012 Emmanuel Paradis.
- DNA-distance and sequence kernels: ape `R/DNA.R`, Copyright 2002-2023
  Emmanuel Paradis, 2015 Klaus Schliep, 2017 Franz Krah; and `src/dist_dna.c`,
  Copyright 2005-2020 Emmanuel Paradis.
- tree depths/distances/path calculations: `src/dist_nodes.c`, Copyright
  2012-2023 Emmanuel Paradis; `src/plot_phylo.c`, Copyright 2004-2017 Emmanuel
  Paradis; `R/nodepath.R`, Copyright 2014 Emmanuel Paradis; and
  `R/branching.times.R`, Copyright 2002-2018 Emmanuel Paradis.
- topology/comparative helpers: `R/is.binary.tree.R`, Copyright 2016-2023
  Emmanuel Paradis; `R/is.ultrametric.R`, Copyright 2003-2016 Emmanuel Paradis;
  `R/vcv.phylo.R`, Copyright 2002-2012 Emmanuel Paradis; `R/is.monophyletic.R`,
  Copyright 2009-2022 Johan Nylander and Emmanuel Paradis; and `R/which.edge.R`,
  Copyright 2004-2017 Emmanuel Paradis, 2017 Joseph W. Brown, 2017 Klaus
  Schliep.
- tree comparison, split collections, clade support, and consensus:
  `R/dist.topo.R`, Copyright 2005-2023 Emmanuel Paradis, 2016-2021 Klaus
  Schliep; `R/as.bitsplits.R`, Copyright 2011-2024 Emmanuel Paradis, 2019 Klaus
  Schliep; `src/bitsplits.c`, Copyright 2005-2024 Emmanuel Paradis; and
  `src/prop_part.cpp`, Copyright 2017 Klaus Schliep.
- tree balance/cherries: `R/balance.R`, Copyright 2002-2015 Emmanuel Paradis,
  2022 Klaus Schliep; and `R/cherry.R`, Copyright 2002-2009 Emmanuel Paradis.
- tree editing: `R/drop.tip.R`, Copyright 2003-2023 Emmanuel Paradis,
  2017-2023 Klaus Schliep, 2018 Joseph Brown; `R/collapse.singles.R`, Copyright
  2015 Emmanuel Paradis, 2017 Klaus Schliep; `R/multi2di.R`, Copyright
  2005-2021 Emmanuel Paradis, 2018-2022 Klaus Schliep; and `R/root.R`, Copyright
  2004-2024 Emmanuel Paradis.
- diversification/coalescent/LTT statistics, including standard and extended birth-death likelihoods: `R/gammaStat.R`, Copyright 2002-2009
  Emmanuel Paradis; `R/yule.R`, Copyright 2003-2011 Emmanuel Paradis;
  `R/birthdeath.R`, Copyright 2002-2022 Emmanuel Paradis;
  `R/coalescent.intervals.R`, `R/collapsed.intervals.R`, and `R/skyline.R`,
  Copyright 2002 Korbinian Strimmer; `R/diversi.time.R`, Copyright 2002-2007
  Emmanuel Paradis; and `R/ltt.plot.R`, Copyright 2002-2021 Emmanuel Paradis.
- principal-coordinate analysis: `R/pcoa.R`, authored by Pierre Legendre
  (October 2007), with literature attribution retained in the translated source
  documentation.
- Moran's I: `R/MoranI.R`, Copyright 2004 Julien Dutheil, 2007-2008 Emmanuel
  Paradis.
- minimum spanning tree: `R/mst.R`, Copyright 2002-2006 Yvonnick Noel, Julien
  Claude, and Emmanuel Paradis.
- quartet delta statistic: `src/delta_plot.c`, Copyright 2010-2011 Emmanuel
  Paradis.
- tree-count and discrete genetic-distance utilities: `R/howmanytrees.R`,
  Copyright 2004-2022 Emmanuel Paradis; and `R/dist.gene.R`, Copyright 2002-2012
  Emmanuel Paradis.
- split compatibility: `R/is.compatible.R`, Copyright 2011 Andrei-Alin
  Popescu.
- diversification goodness-of-fit and sister-clade tests: `R/diversi.gof.R`,
  Copyright 2002-2006 Emmanuel Paradis; and `R/SlowinskiGuyer.R`, Copyright
  2011-2016 Emmanuel Paradis.

The complete upstream author/copyright-holder list is in
`UPSTREAM_DESCRIPTION` and remains controlling attribution for the upstream
work.

## Translation-specific design choices

The R `phylo` list representation is replaced by `type(phylo_tree)` with the
same one-based tip/internal-node numbering convention. R's externally visible
`DNAbin` raw bytes are replaced by language-neutral integer constants for the
17 upstream states (A/C/G/T, IUPAC ambiguity states, N, gap, and unknown); the
translation maps these states back to the historical bit semantics internally
where exact upstream behavior depends on them.
Plotting, GUI/interactive code, R S3/S4 dispatch, Rcpp glue, serialization, and
external-program launchers are intentionally not translated.

No source from `rfortran-core`, `rfortran-linalg`, BLAS, LAPACK, ARPACK, or
another translated R package is copied into this directory. `rfortran-core` is
a sibling FPM path dependency for the repository-wide `dp` kind, and
`rfortran-linalg` is a sibling FPM path dependency for shared eigenvalue and SPD
linear-algebra operations.
