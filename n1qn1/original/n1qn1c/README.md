<!-- badges: start -->
[![Codecov test coverage](https://codecov.io/gh/nlmixr2/n1qn1c/graph/badge.svg)](https://app.codecov.io/gh/nlmixr2/n1qn1c)
[![R-CMD-check](https://github.com/nlmixr2/n1qn1c/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/nlmixr2/n1qn1c/actions/workflows/R-CMD-check.yaml)
[![CRAN version](http://www.r-pkg.org/badges/version/n1qn1)](https://cran.r-project.org/package=n1qn1)
<!-- badges: end -->

*Ported By:* Matthew Fidler, Wenping Wang

*Algorithm Authors:* C. Lemarechal, Stephen L. Campbell, Jean-Philippe
Chancelier, Ramine Nikoukhah

R port of the Scilab n1qn1 module.  This package provides n1qn1, or
Quasi-Newton BFGS "qn" without constraints.  This takes more memory
than traditional L-BFGS.  This routine is useful since it
allows prespecification of a Hessian; If the Hessian is near the truth
in optimization it can speed up the optimization problem.
