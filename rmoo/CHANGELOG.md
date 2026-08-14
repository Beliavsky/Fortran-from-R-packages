# Changelog

## 0.1.0

- Initial modern Fortran/FPM computational translation of rmoo 0.3.2.
- Added NSGA-I, NSGA-II, NSGA-III and R-NSGA-II high-level drivers.
- Added real, binary, discrete integer and permutation representations.
- Added Pareto ranking, crowding, sharing, reference-direction generation,
  NSGA-III association/niching and R-NSGA-II preference truncation.
- Added GD/GD+/IGD/IGD+ metrics and reference-point utilities.
- Reused the supplied GA Fortran RNG/operator layer.
- Added five permanent strict test programs and a ZDT1 example.
- Omitted plotting, parallel R infrastructure and S4/UI code.
