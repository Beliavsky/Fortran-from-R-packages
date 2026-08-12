# Changelog

## 0.1.0

- Initial modern Fortran/FPM translation of ppso 0.9-99994 computational code.
- Added asynchronous and synchronous PSO.
- Added DDS with all four particle-exchange modes.
- Added native random and Latin-hypercube initialization.
- Added initial-estimate and pending-estimate handling.
- Added native optimizer-state checkpoint save/load.
- Added package benchmark functions.
- Added deterministic standalone RNG.
- Added compatibility switch and dual call counts for the serial DDS pre-search
  omission present in the upstream R source.
- Omitted Rmpi transport, plotting, interactive UI, and R object infrastructure.
