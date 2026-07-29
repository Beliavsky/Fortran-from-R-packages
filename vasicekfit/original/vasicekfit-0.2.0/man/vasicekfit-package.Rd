\name{vasicekfit-package}
\alias{vasicekfit-package}
\docType{package}
\title{Extended Vasicek Credit Loss Model with Macroeconomic Factors}
\description{
Fits the extended Vasicek single-factor credit loss model where the
probability of default depends on macroeconomic covariates. Maximum likelihood
estimates of all parameters, including asset value correlation, are obtained
via closed-form probit-transformed OLS regression. Also provides density,
distribution, quantile, and random generation functions for the Vasicek loss
distribution.
}
\author{
Dmitriy Mayorov
}
\references{
Vasicek, O. A. (2002). The distribution of loan portfolio value.
\emph{Risk}, 15(12), 160--162.

Yang, B. H. (2014). Estimating Long-Run PD, Asset Correlation, and
Portfolio Level PD by Vasicek Models. MPRA Paper No. 57244.
}
\seealso{
\code{\link{vasicekfit}}, \code{\link{predict.vasicekfit}},
\code{\link{dvasicek}}
}
