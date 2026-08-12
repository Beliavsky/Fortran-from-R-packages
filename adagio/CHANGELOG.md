# Changelog

## 0.1.0

- Initial modern Fortran/FPM translation of adagio 0.9.2 computational code.
- Added Nelder-Mead, bounded Nelder-Mead, Hooke-Jeeves, simpleEA, simpleDE,
  pureCMAES, and numerical gradients.
- Added subset sum, 0/1 knapsack, approximate bin packing, Hamiltonian search,
  maximum-sum, maximum-empty-rectangle, counting/occurrence and transfinite
  utilities.
- Added test functions and max-quadratic problem objects.
- Added `history_buffer` as the Fortran analogue of `Historize`.
- Integrated the supplied lpSolve-fortran v0.1.0 as a vendored FPM path
  dependency for assignment, change-making, set-cover, and multiple-knapsack
  MILP formulations.
- Added true objective-call counters alongside source-compatible simpleEA/DE
  accounting where useful.
- Fixed two translation-only short-circuit hazards in Hooke-Jeeves and bin
  packing.
- Added six regression executables and two examples.
