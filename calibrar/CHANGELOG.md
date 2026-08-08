# Changelog

## 0.1.0

- Initial modern Fortran/FPM translation of calibrar 0.9.0 computational core.
- Added explicit callback interfaces designed to be portable under `-Werror=implicit-interface`.
- Translated AHR-ES, phased calibration, numerical gradients, fitness functions, random/statistical helpers, splines, objective aggregation, and stopping rules.
- Added standalone optimizer compatibility paths while clearly separating external R-package algorithms from calibrar-owned code.
- Added eight strict regression programs and two examples.
