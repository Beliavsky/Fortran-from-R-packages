# Notice

This project translates the computational code of the R package `frbinom`
1.0.0 by Jeonghwa Lee and Daniel Gernander into modern Fortran/FPM.

The generalized-Bernoulli waiting-time and count-distribution recurrences are
translated directly. Base-R `dbinom`/`runif` and R matrix/list infrastructure
are replaced by native Fortran numerical kernels.
