# Notice

This project is a modern Fortran translation of the computational R-level API of
`nloptr` 2.2.1.9000. The original package is licensed under the GNU Lesser
General Public License, version 3 or later. The translated Fortran sources are
distributed under the same license.

`nloptr` is an interface to the separate NLopt library. This Fortran project does
not copy or compile NLopt's C/C++ solver implementation. It provides portable
Fortran implementations behind compatible high-level algorithm names. See
`PORTING.md` for the exact mapping and limitations.
