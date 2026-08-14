# Upstream provenance

This project is a modern-Fortran translation of the computational code in the supplied KrigInv source archive:

- Package: KrigInv
- Upstream version: 1.4.2
- Date: 2022-09-06
- License: GPL-3
- Authors: Clement Chevalier, Dario Azzimonti, David Ginsbourger, Victor Picheny; Yann Richet (contributor)

The original KrigInv `DESCRIPTION` and `NAMESPACE` are retained as `UPSTREAM_DESCRIPTION` and `UPSTREAM_NAMESPACE`.

Version 0.2.0 additionally incorporates the numerical subset required from the previously translated DiceKriging package:

- Package: DiceKriging
- Upstream version: 1.6.1
- Upstream license: GPL-2 | GPL-3
- License option used in this combined KrigInv project: GPL-3
- Included modules: kinds, linear algebra, covariance/scaling, bounded optimizer, and `km` fit/predict/update core

DiceKriging provenance, original metadata, and license texts are retained under `licenses/DiceKriging/`.

Other third-party computational material included in the self-contained port:

- GPL-3 Fortran `anMC` modules for conservative excursion sets;
- Sobol direction-number data translated from randtoolbox/its upstream low-discrepancy source, retaining BSD notices.

No R, Rcpp, C, or C++ runtime is required.
