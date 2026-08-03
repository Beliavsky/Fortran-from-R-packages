#'@title Computes the Moebius Cramer-von Mises statistics for the test of randomness
#'
#'@description This function he Moebius Cramer-von Mises statistics for a tests of randomness for random vectors Y(1), ...l, X(p).
#'
#'@param Y            Time series.
#'@param p            Number of consecutive observations for the test.
#'@param trunc.level  Only subsets of cardinality <= trunc.level (default=2) are considered for the Moebius statistics.
#
#'
#'
#'@return \item{stats}{Cramer-von Mises Moebius statistics}
#'@return \item{cardA}{Cardinality of subsets}
#'@return \item{M}{Matrix for multitpliers bootstrap for stats}
#'@return \item{Asets}{vector of (0,1) for Moebius subsets}
#'@return \item{Sn}{Cramer-von Mises Sn statistic}
#'@return \item{J}{Matrix for multipliers bootstrap for Sn}
#'
#'@references B.R Nasri (2022). Tests of serial dependence for arbitrary distributions
#'

#'#'@examples
#' Y <- data(Y)
#' out <- Sn_Aserial(Y,5,2)

#'@keywords internal
#'
#'@export
#'
#'

Sn_AserialVec = function(Y,p,trunc.level=2){
  v = c(1:(trunc.level-1))
  cA = sum(choose((p-1),v))
  dd = dim(Y)
  n = dd[1]
  d = dd[2]

  out0 = .C("stats_serialVectors",
            as.double(Y),
            as.integer(n),
            as.integer(d),
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
