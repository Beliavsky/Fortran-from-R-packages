##   file partial.mk.test.R part of package trend
##   Copyright (C) 2015-2018 Thorsten Pohlert
##
##    This program is free software: you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation, either version 3 of the License, or
##    (at your option) any later version.
##
##    This program is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with this program.  If not, see <http://www.gnu.org/licenses/>.
##
#' @title Partial Mann-Kendall Trend Test
#' @description Performs a partial Mann-Kendall Trend Test
#'
#' @param x a "vector" or "ts" object that contains the variable,
#' which is tested for trend (i.e. correlated with time)
#' @param y a "vector" or "ts" object that contains the variable,
#' which effect on "x" is partialled out
#' @param alternative character, the alternative method; defaults to
#' "two.sided"
#' 
#' @details
#' According to Libiseller and Grimvall (2002), the test statistic
#' for \code{x} with its covariate \code{y} is
#'
#' \deqn{
#' z = \frac{S_x - r_{xy} S_y}{\left[ \left( 1 - r_{xy}^2 \right)
#' n \left(n - 1 \right) \left(2 n + 5 \right) / 18 \right]^{0.5}}}{%
#' z = (S[x] - r * S[y]) / \sqrt(1 - r^2) *
#' n * (n - 1) * (2 * n + 5) / 18)}
#' 
#' where the correlation \eqn{r} is calculated as:
#'
#' \deqn{
#' r_{xy} = \frac{\sigma_{xy}}
#' {n \left(n - 1\right) \left(2 n + 5 \right) / 18}
#' }{%
#' r = \sigma /
#' (n (n - 1) (2n + 5) / 18)}
#'
#' The conditional covariance between \eqn{x} and \eqn{y} is
#'
#' \deqn{
#' \sigma_{xy} = \frac{1}{3}
#' \left[K + 4 \sum_{j=1}^n R_{jx} R_{jy} -
#' n \left(n + 1 \right) \left(n + 1 \right) \right]}
#'
#' with
#'
#' \deqn{
#' K = \sum_{1 \le i < j \le n} \mathrm{sgn} \left\{ \left( x_j - x_i \right)
#' \left( y_j - y_i \right) \right\}}
#'
#' and
#'
#' \deqn{
#' R_{jx} = \left\{ n + 1 + \sum_{i=1}^n
#' \mathrm{sgn} \left( x_j - x_i \right) \right\} / 2}
#'
#' @return
#' A list with class "htest"
#' 
#'  \item{method}{
#'    a character string indicating the chosen test
#'  }
#'  \item{data.name}{
#'    a character string giving the name(s) of the data
#'  }
#'  \item{statistic}{
#'    	the value of the test statistic
#'  }
#'  \item{estimate}{
#'    the Mann-Kendall score S, the variance varS and the correlation
#'    between x and y
#'  }
#'  \item{alternative}{
#'    a character string describing the alternative hypothesis
#'  }
#'  \item{p.value}{
#'    the p-value of the test
#'  }
#' \item{null.value}{
#' the null hypothesis
#' }
#'
#' @references
#' Libiseller, C. and Grimvall, A., (2002).
#' Performance of partial Mann-Kendall tests for trend detection in the
#' presence of covariates. Environmetrics 13, 71--84,
#' \doi{10.1002/env.507}.
#'
#' @note 
#'  Current Version is for complete observations only.
#'  The test statistic is not corrected for ties.
#'
#' @seealso
#'  \code{\link{partial.cor.trend.test}},
#'
#' @examples
#' data(maxau)
#' s <- maxau[,"s"]; Q <- maxau[,"Q"]
#' partial.mk.test(s,Q)
#'
#' @keywords ts nonparametric
#' @importFrom stats pnorm
#' @export
partial.mk.test <- function(x, y,
                            alternative = c("two.sided", "greater", "less"))
{
    ## Tests
    na.fail(x)
    na.fail(y)
    nx <- length(x)
    ny <- length(y)
    if (nx != ny) {
        stop("Error: objects x and y must be of same length")
    }
    alternative <- match.arg(alternative)

    ## See Libiseller and Grimvall (2002, p. 73--76)
    n <- nx

    Sx <- .mkScore(x)
    Sy <- .mkScore(y)

    Rx <- .R(x)
    Ry <- .R(y)
    K <- .K(x, y)

    sigma <- (K + 4 * sum(Rx * Ry) - n * (n + 1)^2) / 3
    rho <- sigma / (n * (n - 1) * (2 * n + 5 ) / 18)
    S <- Sx - rho * Sy
    
    varS <- ( 1 - rho^2) * n * (n - 1) * (2 * n + 5) / 18

    ##  test statistic
    z <- S / sqrt(varS)
    
    pval <- switch(alternative,
                   "two.sided" = {
                       2 * min(0.5 , pnorm(abs(z), lower.tail=FALSE))
                   },

                   "greater" = {
                       pnorm(z, lower.tail=FALSE)
                   },
                   
                   "less" = {
                       pnorm(z, lower.tail=TRUE)
                   }
                   )

    DNAMEX <- deparse(substitute(x))
    DNAMEY <- deparse(substitute(y))
    DNAME <- paste("t AND ",DNAMEX," . ",DNAMEY)
    
    ## prepare output
    res <- list(statistic= c(z = z),
                null.value = c(S = 0),
                alternative = alternative,
                estimates = c(S = S, varS = varS, cor = rho),
                p.value = pval,
                method = "Partial Mann-Kendall Trend Test",
                data.name = DNAME)
    class(res) <- "htest"
    return(res)
}
