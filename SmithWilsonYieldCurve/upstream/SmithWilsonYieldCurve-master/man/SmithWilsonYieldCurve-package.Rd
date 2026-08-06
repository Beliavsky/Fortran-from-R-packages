\docType{package}
\name{SmithWilsonYieldCurve-package}
\alias{SmithWilsonYieldCurve}
\alias{SmithWilsonYieldCurve-package}
\title{Fit yield curves using the Smith-Wilson method}
\description{
  A package to fit yield curves using the Smith-Wilson
  method
}
\details{
  The main function exposed in this package is
  fFitSmithWilsonYieldCurve, which takes market data in the
  form of a vector of cashflow times, a matrix of cashflows
  and a vector of market prices. It returns an object of
  class "SmithWilsonYieldCurve".

  A convenience function
  fFitSmithWilsonYieldCurveToInstruments takes a dataframe
  containing market instrument data as type, tenor,
  frequency and rate. It extracts the required vectors and
  matrices and then calls fFitSmithWilsonYieldCurve.

  Objects of class SmithWilsonYieldCurve are a list, the
  first element of which is a function P(t), which returns
  the zero coupon bond price of the fitted curve at time t.

}
\examples{
dfInstruments <- data.frame(c("SWAP", "SWAP"), c(1,10), c(1,1), c(0.025, 0.05))
	colnames( dfInstruments ) <- c( "Type", "Tenor", "Frequency", "Rate" )
	Curve <- fFitSmithWilsonYieldCurveToInstruments( dfInstruments, 0.04, 0.1 )
	plot( Curve )
}
\author{
  Phil Joubert
  \email{phil.joubert@gmail.com}
}
\references{
  Smith A. and Wilson, T. - "Fitting Yield curves with long Term Constraints" (2001)
}

