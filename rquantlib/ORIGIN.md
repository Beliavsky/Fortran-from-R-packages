# Upstream origin

- Package: RQuantLib
- Upstream version: 0.4.28
- Package date: 2026-06-10
- Attached archive: `RQuantLib-master.zip`
- Attached archive SHA-256: `00a25495b33558ce67217b1aaf76ef50bd3c28607b85528b7bd20cf9e8c7ddf2`
- Upstream declared license: `GPL (>= 2)`
- Upstream system dependencies: QuantLib and Boost

RQuantLib is an R/Rcpp interface layer around QuantLib. Most pricing algorithms invoked by its C++ wrappers are implemented in the external QuantLib library, not in the attached RQuantLib repository. This translation is therefore an independent Fortran implementation of the feasible numerical surface rather than a line-by-line translation of QuantLib.

No QuantLib or Boost source code is included or linked. Their upstream license texts are retained in `licenses/` because the attached package distributed those notices and used those libraries.
