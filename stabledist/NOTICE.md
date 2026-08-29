# Notices and attribution

This project is a modern Fortran translation of the computational portions of
R package `stabledist` 0.7-2.

Upstream authors/maintainers credited by `DESCRIPTION`:

* Diethelm Wuertz — original code
* Martin Maechler — checks, tests, fixes, and numerical improvements
* Yohan Chalabi — namespace/administrative contributions

The upstream package states `License: GPL (>= 2)`. Translated stabledist code
is therefore distributed as GPL-2.0-or-later. Historical source headers that
refer to the GNU Library General Public License are retained verbatim in the
archived upstream files.

The supplied `r_mod.f90` is separately MIT-licensed and remains under those
terms.

Algorithm/reference attribution retained from the upstream documentation:

* Chambers, Mallows and Stuck (1976), "A Method for Simulating Stable Random
  Variables", JASA 71, 340-344.
* J. P. Nolan (1997), "Numerical calculation of stable densities and
  distribution functions", Stochastic Models 13(4), 759-774.
* J. P. Nolan (2020), *Univariate Stable Distributions - Models for Heavy
  Tailed Data*, Springer.
* Samorodnitsky and Taqqu (1994), *Stable Non-Gaussian Random Processes*.
* Weron and Weron (1999), stable-variable simulation reference cited upstream.
* Royuela-del-Val, Simmross-Wattenberg and Alberola-Lopez (2017), libstable.

Original metadata, R sources, tests, extra Levy formulas, and package
attribution are retained under `upstream/`.
