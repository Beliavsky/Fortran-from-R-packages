#'@title Computes the Moebius Cramer-von Mises statistics for the test of independence between random variables
#'
#'@description This function computes Moebius Cramer-von Mises statistics for a tests of independence for matrix of data x.
#'
#'@param x         Matrix of observations.
#'@param trunc.level  Only subsets of cardinality <= trunc.level (default=2) are considered for the Moebius statistics.
#
#'
#'@return \item{stats}{Cramer-von Mises Moebius statistics}
#'@return \item{cardA}{Cardinality of subsets}
#'@return \item{M}{Matrix for multitpliers bootstrap for stats}
#'@return \item{Asets}{Vector of (0,1) for Moebius subsets}
#'@return \item{Sn}{Cramer-von Mises Sn statistic}
#'@return \item{J}{Matrix for multipliers bootstrap for Sn}
#'
#'@references Genest, Neslehova, Remillard & Murphy (2019). Testing for independence in arbitrary distributions
#'

#'#'@examples
#' X <- rnorm(100,5)
#' out <- Sn_A(X,3)

#'@keywords internal
#'
#'@export
#'
#'



Sn_A = function(x,trunc.level){
  v = c(2:(trunc.level))
  dim0=dim(x)
  n = dim0[1]
  d = dim0[2]
  cA = sum(choose(d,v))

  out0 = .C("stats_nonserial",
            as.double(x),
            as.integer(n),
            as.integer(d),
            as.integer(trunc.level),
            stats = double(cA),
            cardA = double(cA),
            M = double(n*n*cA),
            Asets = double(d*cA),
            Sn    = double(1),
            J = double(n*n),
            PACKAGE = "MixedIndTests"
  )
}
