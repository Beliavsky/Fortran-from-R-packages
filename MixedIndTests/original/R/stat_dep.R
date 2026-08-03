#'@title Computes the Kendall's taus and Spearman's rho for tests of randomness
#'
#'@description This function Computes the Kendall's taus and Spearman's rho for tests of independence used in EstDep
#'@param x     Vector of length n
#'@param y     Vector of length n
#'
#'
#'@return \item{tau}{Kendall's taus for lags 1:lag}
#'@return \item{rho}{Spearman's rhos for lags 1:lag}
#'@return \item{s2}{Estimated variance of Spearman's rho}
#'
#'@references Genest, Neslehova, & Remillard (2017). Asymptotic behavior of the empirical multilinear copula process under broad conditions
#'
#'@examples
#' x <- SimAR1Poisson(c(5,0.2),100)
#' y <- SimAR1Poisson(c(5,0.2),100)
#' out <- stat_dep(x,y)

#'
#'@keywords internal
#'
#'@export
#'
#'


stat_dep =function(x,y){
  n = length(x)

  out0 = .C("estdep",
            as.double(x),
            as.double(y),
            as.integer(n),
            tau = double(1),
            rho = double(1),
            s2 = double(1),
            PACKAGE = "MixedIndTests"

  )
}
