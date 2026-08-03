#'@title Computes the Cramer-von Mises statistic Sn for the test of randomness
#'
#'@description This function computes the Cramer-von Mises statistic for a tests of randomness for p consecutives values X(1), ..., X(p). Useful for traditional bootstrapping.
#'
#'@param x        Time series
#'@param p        Number of consecutive observations
#
#'
#'
#'@return \item{stats}{Cramer-von Mises statistics Sn}
#'
#'@references B.R Nasri (2022). Tests of serial dependence for arbitrary distributions
#'

#'#'@examples
#' X <- SimAR1Poisson(c(5,0.2),100)
#' out <- Sn_serial(X,5,3)

#'@keywords internal
#'
#'@export
#'
#'


Sn_serial = function(x,p){

  n = length(x)

  out0 = .C("Sn_serial0",
            as.double(x),
            as.integer(n),
            as.integer(p),
            Sn    = double(1),
            PACKAGE = "MixedIndTests"
  )
}
