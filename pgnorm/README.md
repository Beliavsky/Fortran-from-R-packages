# pgnorm-fortran

Modern Fortran translation of the computational code in the R package `pgnorm` 2.0.1.

## API

- `dpgnorm`, `ppgnorm`
- `rpgnorm`
- `rpgnorm_nardonpianca`
- `rpgnorm_pgenpolar`
- `rpgnorm_pgenpolarrej`
- `rpgnorm_montypython`
- `rpgnorm_ziggurat`
- `rpgangular`, `rpgunif`
- `zigsetup`
- `pgnorm_sd`

The default parameterization is the one implemented by the R source: `scale=1` gives

`f(x) = C_p exp(-|x|^p/p)`.

Its standard deviation is returned by `pgnorm_sd(p, scale)`.

Build with `fpm test` or compile with a Fortran 2018 compiler.
