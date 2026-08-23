# Notice

This project is a modern Fortran translation of computational routines from **rvMF 0.1.2**, by Seungwoo Kang and Hee-Seok Oh.

Upstream source is preserved under `upstream/rvMF-master/`. The translation is not affiliated with or endorsed by the upstream authors.

The upstream R wrapper uses `Rfast::matrnorm` to generate independent standard normal variates. The supplied Rfast Fortran translation was inspected for dependency parity; this project implements that elementary RNG operation locally and therefore does not copy Rfast source code or require Rfast at build time.

The core non-rejection algorithm follows the upstream `rvMF64.cpp` implementation and Kang & Oh (2024), including integer probability truncation and condensed base-64 lookup tables.
