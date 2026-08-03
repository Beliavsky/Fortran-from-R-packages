#'@title Computes the Moebius Cramer-von Mises statistics for the test of randomness
#'
#'@description This function computes Moebius Cramer-von Mises statistics for a tests of randomness for observations X(1), ..., X(p).
#'
#'@param x            Time series.
#'@param p            Number of consecutive observations for the test.
#'@param trunc.level  Only subsets of cardinality <= trunc.level (default=2) are considered for the Moebius statistics.
#
#'
#'
#'@return \item{stats}{Cramer-von Mises Moebius statistics}
#'@return \item{cardA}{Cardinality of subsets}
#'@return \item{M}{Matrix for multitpliers bootstrap for stats}
#'@return \item{Asets}{Vector of (0,1) for Moebius subsets}
#'@return \item{Sn}{Cramer-von Mises Sn statistic}
#'@return \item{J}{Matrix for multipliers bootstrap for Sn}
#'
#'@references B.R Nasri (2022). Tests of serial dependence for arbitrary distributions
#'

#'#'@examples
#' X <- SimAR1Poisson(c(5,0.2),100)
#' out <- Sn_Aserial(X,5,3)

#'@keywords internal
#'
#'@export
#'
#'

Sn_Aserial = function(x,p,trunc.level){
  v = c(1:(trunc.level-1))
  cA = sum(choose((p-1),v))
  n = length(x)

  out0 = .C("stats_serial",
            as.double(x),
            as.integer(n),
            as.integer(p),
            as.integer(trunc.level),
            stats = double(cA),
            cardA = double(cA),
            M = double(n*n*cA),
            Asets = double(p*cA),
            Sn    = double(1),
            J = double(n*n),
            PACKAGE = "MixedIndTests"
  )
}
