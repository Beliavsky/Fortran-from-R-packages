# R Interface for PIQP <img src="man/figures/logo.png" width="100" align="right" />

[![R-CMD-check](https://github.com/PREDICT-EPFL/piqp-r/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/PREDICT-EPFL/piqp-r/actions/workflows/R-CMD-check.yaml)
[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/piqp)](https://cran.r-project.org/package=piqp)
[![](https://cranlogs.r-pkg.org/badges/piqp)](https://CRAN.R-project.org/package=piqp)

R interface for [PIQP](https://predict-epfl.github.io/piqp/) a Proximal Interior Point Quadratic Programming solver.

This package/interface is partly based on [osqp-r](https://github.com/osqp/osqp-r).

Stable versions can be installed from CRAN as usual. Development
versions from this repo can be installed via:

```
## Install remotes packages if not available
if (! "remotes" %in% installed.packages()[, 1] ) {
	install.packages("remotes", repository = "https://cran.r-project.org")
}
remotes::install_github("PREDICT-EPFL/piqp-r")
```
