# scip 1.10.0-3

- Upgrade to SCIP 10.0.2, SoPlex 8.0.2, PaPILO 3.0.0.
- Enable OpenMP thread pool interface (TPI=omp) when the platform
  supports it, giving SCIP parallel branch-and-bound. Falls back
  gracefully to TPI=none when OpenMP is unavailable.
- Use `SHLIB_OPENMP_CXXFLAGS` in both `PKG_CXXFLAGS` and `PKG_LIBS`
  per R-exts §1.2.1.1.
- Drop all tinycthread patches (no longer compiled with TPI=omp/none).
  Reduces R-specific patch burden from 14 to 10 across submodules.

# scip 1.10.0-1

- Switched build system from hand-maintained Makevars.in (472 lines) to
  CMake-based build at install time, following the highs R package pattern.
- Solver sources (SCIP, SoPlex) moved from `src/` to `inst/` and deleted
  after compilation when installing from tarball, reducing installed size
  from ~258 MB to ~10 MB.
- Version now tracks SCIP Optimization Suite (10.0.x).

# scip 0.0.2

Initial CRAN submission.

- One-shot solver interface (`scip_solve`) and incremental
  model-building API (`scip_model`, `scip_add_var`,
  `scip_add_linear_cons`, `scip_add_quadratic_cons`,
  `scip_add_sos1_cons`, `scip_add_sos2_cons`,
  `scip_add_indicator_cons`).
- Solver control parameters (`scip_control`).
- Sparse matrix support (dgCMatrix, simple_triplet_matrix).
- Vignette with LP, MIP, quadratic, and indicator constraint examples.
- Builds on macOS, Linux, and Windows using vendored SCIP 10.0.1 and
  SoPlex 8.0.1 sources.
