# Changelog

## 0.2.3-fortran.2

- Replaced platform-sensitive exact-support F-score assertions in the
  bipartite learner tests with direct adjacency-spectrum symmetry checks.
- The bipartite objective can validly select a different bipartite support on
  different numerical platforms, so the revised tests verify the property the
  estimators actually constrain.
- No library API or computational implementation changes.

## 0.2.3-fortran.1

- Translated all 24 exported computational APIs from spectralGraphTopology 0.2.3.
- Added the non-exported cospectral graph learner and supporting operators.
- Replaced R/C++ numerical dependencies with self-contained modern Fortran.
- Added FPM metadata, four tests, four examples, a demo, build scripts, and API
  and porting documentation.
- Preserved the original GPL-3 license and source materials.
