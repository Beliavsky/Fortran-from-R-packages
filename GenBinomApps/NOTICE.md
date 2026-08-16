# Notice

This project translates the computational code of the R package
`GenBinomApps` 1.2.1 by Horst Lewitschnig and David Lenzi into modern
Fortran/FPM.

The generalized-binomial convolution and Clopper-Pearson formulas are
preserved. Base-R beta, root-finding, and uniform-RNG dependencies are replaced
by standalone Fortran implementations. See `PORTING_NOTES.md`.
