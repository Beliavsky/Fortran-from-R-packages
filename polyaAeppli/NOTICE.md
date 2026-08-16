# Notice

This project translates the computational code of the R package
`polyaAeppli` 2.0.2 by Conrad Burden into modern Fortran/FPM.

The core Johnson-Kotz-Kemp/Nuel probability recurrence and the package's
special upper-tail recurrence are preserved. Base-R random generation and
gamma-quantile bracketing are replaced by standalone Fortran equivalents.
