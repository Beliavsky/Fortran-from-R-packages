# Notice

This project translates the computational code of `mcmc` 0.9-8 by Charles J.
Geyer and Leif T. Johnson into modern Fortran/FPM.

The Metropolis, tempering, initial-sequence, and overlapping-batch-means
algorithms are derived from the supplied R/C sources. R S3/list/environment
infrastructure is replaced by native Fortran callbacks and derived types.
