#'@title Computes the Kendall's taus and Spearman's rho for tests of randomness
#'
#'@description This function Computes the Kendall's taus and Spearman's rho for tests of randomness.
#'
#'@param x          Time series
#'@param lag        Number of lags.
#
#'
#'@return \item{tau}{Kendall's taus for lags 1:lag}
#'@return \item{rho}{Spearman's rhos for lags 1:lag}
#'
#'@references B.R Nasri (2021). Tests of serial dependence for arbitrary distributions
#'
#'@examples
#' X <- SimAR1Poisson(c(5,0.2),100)
#' out <- stat_dep_ser(X,10)

#'
#'@keywords internal
#'
#'@export
#'
#'


stat_dep_ser =function(x,lag){
  n = length(x)

  out0 = .C("estdep_serial",
            as.double(x),
            as.integer(n),
            as.integer(lag),
            tau = double(1),
            rho = double(1),
            s2 = double(1),
            PACKAGE = "MixedIndTests"

  )
}
