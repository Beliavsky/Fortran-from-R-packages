# Changelog

## 0.1.0

- Initial modern Fortran/FPM translation of DEoptim 2.2-8 computational code.
- Implemented all six exposed differential-evolution strategies.
- Implemented JADE/current-to-p-best controls `p` and `c`.
- Implemented bounds, initialization, mapped populations and duplicate removal,
  convergence controls, histories, stored populations, and evaluation counts.
- Added portable standalone uniform/normal/Cauchy random-number support.
- Added tests covering all six strategies, mapped/discrete optimization,
  storage, initial-population/VTR stopping, and DEoptim's unsupported `bs` flag.
- Added Rosenbrock and integer-map examples.
- Preserved GPL-2.0-or-later licensing and original package provenance.
