# Validation

The translated source is validated with GNU Fortran using:

```text
-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all -O0
```

Regression programs cover:

1. dense median Frisch-Newton regression;
2. inequality-constrained QR and multiple taus;
3. lasso and SCAD paths;
4. Hyndman-Fan quantiles and order-statistic selection;
5. local linear QR;
6. nonlinear QR callback path;
7. XY bootstrap;
8. large-n preprocessing, recursive least squares and combinations.

`fpm.toml` is parsed with Python's standard TOML parser. The release process
also performs a clean unzip and recompiles/tests only the archive contents.
