# Changelog

## 0.1.0

- Initial modern Fortran/FPM translation of the computational interface in
  R package `scip` 1.10.0-3.
- Added `scip_solve` for dense and CSC linear/MIP models.
- Added persistent model-building API for continuous, binary, and integer
  variables.
- Added linear, quadratic, SOS1, SOS2, and indicator constraints.
- Added objective-sense and native SCIP parameter setters.
- Added common `scip_control` fields and arbitrary typed extra parameters.
- Added best-solution, solution-pool, status, gap, node, LP-iteration, and
  solve-time queries.
- Added a plain-C ABI shim with explicit `iso_c_binding` interfaces.
- Retained SCIP 10.0.2 and SoPlex 8.0.2 as vendored source backends rather
  than substituting a different optimizer.
- Added POSIX and MinGW/PowerShell vendor-build helpers.
- Added five strict regression tests and two examples.
