# Changelog

## 0.1.1 - 2026-08-06

- Route objective, gradient, proposal, and monitor callbacks through module-level
  trampoline procedures with explicit interfaces. This fixes
  `-Werror=implicit-interface` failures seen with some GNU Fortran/FPM builds.
- No numerical algorithm changes.

## 0.1.0

- Initial modern Fortran computational port of roptim 0.1.7.
- Added Nelder-Mead, BFGS, three nonlinear-CG update rules, L-BFGS-B 3.0,
  and simulated annealing.
- Added analytic and finite-difference gradient support.
- Added numerical Hessians, scaling, maximization, bounds, user data,
  monitor cancellation, deterministic SANN seeds, and custom proposals.
- Added FPM metadata, tests, examples, provenance, and license files.
