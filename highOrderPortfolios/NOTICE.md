# Notice

This is a modern Fortran translation of `highOrderPortfolios` 0.1.1 by Rui
Zhou, Xiwen Wang, and Daniel P. Palomar. The upstream package is licensed under
GPL version 3.

The skew-t estimation modules are adapted from the previously produced
`fitHeavyTail-fortran` translation, itself based on the GPL-3.0-only R package
`fitHeavyTail` by Xiwen Wang and Daniel P. Palomar.

The upstream repository includes copied internal routines from
PerformanceAnalytics, with copyright notices for Kris Boudt, Brian G.
Peterson, Peter Carl, and other contributors. Those original C sources and
notices are retained under `original/src/`; this Fortran implementation uses
direct central-moment contractions rather than compiling that C code.

This translation is a modified work. It is provided without warranty under the
terms in `LICENSE`.
