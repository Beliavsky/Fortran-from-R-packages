# Validation

The test program checks:

- quantile-correlation values against independently calculated reference
  values;
- agreement between WQC estimates and an independent implementation of the
  quantile-correlation formula applied to `waveslim` MRA details;
- result dimensions and finite values;
- ordered lower and upper confidence limits;
- exact repeatability when the same RNG seed is supplied;
- multi-series analysis and series names;
- rejection of invalid input dimensions.

The package was compiled and tested with GNU Fortran 14.2.0 using both a debug
configuration with runtime checking and an optimized `-O3` configuration.

FPM was not installed in the build environment. The `fpm.toml` files were
parsed as TOML, and the same ordered sources and bundled dependencies were
compiled directly with GNU Fortran.
