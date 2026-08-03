\name{multiAssetOptions-package}
\alias{multiAssetOptions-package}
\alias{multiAssetOptions}
\docType{package}
\title{
Finite Difference Method for Multi-Asset Option Valuation
}
\description{
  This package implements an efficient finite difference method for valuing multi-asset options in the Black-Scholes world. The model assumes static volatilities and correlations. The implementation allows users to vary the option setup (number of underlying assets, call vs. put, European vs. American, etc.) as well as some of the features of the numerical method (mesh spacing, timestepping scheme, etc.).
}
\details{
\tabular{ll}{
Package: \tab multiAssetOptions\cr
Type: \tab Package\cr
Version: \tab 0.1-2\cr
Date: \tab 2021-04-20\cr
License: \tab GPL-2 | GPL-3\cr
}
  This package contains the function \code{multiAssetOption}, which implements a generalized version of the finite difference method for option valuation. Several sub-functions are included in this package to improve code clarity and presentation.
}
\author{
  Michael Eichenberger and Carlo Rosa

  Maintainer: Michael Eichenberger <mike.eichenberger@gmail.com>
}
\keyword{ package }
