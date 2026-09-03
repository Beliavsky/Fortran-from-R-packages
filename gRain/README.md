# gRain — modern Fortran computational translation

This directory is a modern free-form Fortran translation of the portable
computational core of the R package **gRain 1.4.6**, a package for probability
propagation in Bayesian networks.

The translation preserves the Bayesian-network and junction-tree numerical
algorithms while replacing R objects, factor labels, formulas, S3 methods,
Rcpp, Eigen, and Armadillo interfaces with a small native Fortran API. Nodes
are numbered `1:n` and discrete states are numbered `1:cardinality(node)`.
Multidimensional probability tables use the same column-major ordering as R
and Fortran.

## Implemented computational areas

- conditional-probability-table construction and first-dimension normalization;
- Boolean AND/OR CPT helpers and Mendelian segregation probabilities;
- DAG validation, moralization, optional root completion, weighted minimal
  triangulation, and running-intersection/junction-tree construction;
- clique-potential construction and CPT insertion;
- Lauritzen-Spiegelhalter collect/distribute propagation;
- hard and soft evidence insertion, replacement, retraction, and absorption;
- probability of evidence;
- marginal, joint, and conditional queries, including non-clique joint
  reconstruction from calibrated clique/separator potentials;
- forward simulation from a calibrated network, including evidence;
- empirical discrete tables and CPT estimation from integer categorical data;
- clique-marginal estimation and marginal/potential conversion; and
- CPT replacement without rebuilding the triangulation when graph structure is
  unchanged.

See `API_COVERAGE.md` for a detailed upstream-to-Fortran mapping and explicit
parity boundaries.

## Dependencies

The package is intended to live at repository root beside the shared packages:

```text
Fortran-from-R-packages/
  gRain/
  gRbase/
  rfortran-core/
  ...
```

`fpm.toml` uses sibling path dependencies:

```toml
rfortran-core = { path = "../rfortran-core" }
gRbase = { path = "../gRbase" }
```

The `gRbase` translation supplies column-major probability-table operations,
graph algorithms, triangulation/RIP support, and its own `igraph` and
`rfortran-linalg` dependencies. No dependency source, BLAS, LAPACK, ARPACK,
Rcpp, Eigen, or Armadillo source is copied into this package.

## Build and test

From this directory, with the sibling dependencies present:

```text
fpm build
fpm test
fpm run --example grain_demo
```

The maintained Fortran uses a single working precision `dp` imported from
`r_kinds` through `rfortran-core`. The code is free-form Fortran and uses
explicit interfaces and explicit dummy-argument intents.

The validation environment used to create this archive did not contain the
`fpm` or `fprettify` executables. `BUILD_VALIDATION.md` records the strict
GNU Fortran validation that was performed instead and the exact limitation.

## Minimal example

The example `example/grain_demo.f90` constructs the binary chain
`1 -> 2 -> 3`, propagates it, conditions on node 3 being in state 1, and
queries the posterior distribution of node 1.

The CPT convention is child variable first. For example,

```fortran
cpts(2) = make_cpt([2, 1], [2, 2], &
                   [0.7_dp, 0.3_dp, 0.2_dp, 0.8_dp])
```

represents `P(node 2 | node 1)`. Each consecutive block of length two is a
conditional distribution over node 2.

## Random simulation

The R package draws from R's RNG. This standalone Fortran translation instead
provides two portable alternatives:

- `simulate_from_uniforms`, which accepts caller-supplied uniforms and is the
  preferred route when an external RNG or exact reproducibility is required;
- `simulate_network`, which uses a deterministic portable Park-Miller stream
  and accepts an optional integer seed.

The RNG stream is therefore intentionally not bit-for-bit compatible with R.
The conditional sampling probabilities are the translated gRain computation.

## License and provenance

The upstream package declares `GPL (>= 2)`. `LICENSE`, `COPYING`, `NOTICE.md`,
`upstream/DESCRIPTION`, and `upstream/CITATION` preserve the applicable
licensing, attribution, and package provenance. This translation is distributed
under GPL-2.0-or-later consistently with upstream.
