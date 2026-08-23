# Changelog

## 0.1.0

- Initial modern Fortran/FPM computational port of Rnanoflann 0.0.3.
- Added `nn()` with standard and radius searches and R-compatible matrix orientation.
- Added all 21 distance/dissimilarity formulas implemented by the upstream C++ source.
- Added raw-input Hellinger preprocessing, Euclidean/Hellinger squared modes, and Minkowski `p` support.
- Added optional OpenMP query parallelism.
- Corrected the upstream radius-query pointer bug.
- Prevented radius result-buffer overflow by capping stored output at `k` and reporting total matches separately.
- Honored the radius `sorted` argument, which is ignored by the upstream C++ call path.
- Added strict numerical and API tests plus an example program.
