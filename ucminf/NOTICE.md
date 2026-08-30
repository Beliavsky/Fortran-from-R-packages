# Upstream attribution and provenance

This project is a modern Fortran translation of the computational code in the
R package **ucminf 1.2.3**, licensed upstream as **GPL (>= 2)**.

Upstream authors listed in `DESCRIPTION`:

- Hans Bruun Nielsen — author of the UCMINF algorithm and original Fortran code.
- Stig Bousgaard Mortensen — R implementation.
- K Herve Dakpo — contributor and current upstream maintainer in the supplied snapshot.

The supplied upstream documentation also credits modifications by Douglas Bates
and Tomas Kalibera in the R/native interface history.

The original algorithm is described in H. B. Nielsen (2000), *UCMINF - An
Algorithm For Unconstrained, Nonlinear Optimization*, Department of Mathematical
Modelling, Technical University of Denmark.

For auditability, the relevant original fixed-form Fortran, R wrapper, C bridge,
README, and DESCRIPTION from the user's source archive are retained under
`upstream/`.
