#' A package to fit yield curves using the Smith-Wilson method
#' 
#' The main function exposed in this package is fFitSmithWilsonYieldCurve, which takes market data in the form of a vector of cashflow times, a matrix of cashflows and a vector of market prices. It returns an object of class "SmithWilsonYieldCurve". 
#' 
#' A convenience function fFitSmithWilsonYieldCurveToInstruments takes a dataframe containing market instrument data as type, tenor, frequency, rate and price. It extracts the required vectors and matrices and then calls fFitSmithWilsonYieldCurve.
#' 
#' Objects of class SmithWilsonYieldCurve are a list, the first element of which is a function P(t), which returns the zero coupon bond price of the fitted curve at time t.
#' 
#' For details including mathematics check the EIOPA document in references.
#' 
#' @name SmithWilsonYieldCurve-package
#' @aliases SmithWilsonYieldCurve
#' @docType package
#' @title Fit yield curves using the Smith-Wilson method
#' @author Phil Joubert \email{phil.joubert@@gmail.com}
#' @references The original CEIOPS paper referenced by this package has been removed from the EIOPA website, but googling "ceiops paper extrapolation risk free rates smith wilson" will get you to it.
#' @examples
#' 	dfInstruments <- data.frame(c("SWAP", "SWAP"), c(1,10), c(1,1), c(0.025, 0.05))
#' 	colnames( dfInstruments ) <- c( "Type", "Tenor", "Frequency", "Rate" )
#' 	Curve <- fFitSmithWilsonYieldCurveToInstruments( dfInstruments, 0.04, 0.1 )
#' 	plot( Curve )
NA