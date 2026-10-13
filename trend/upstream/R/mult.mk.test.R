##    file mult.mk.test.R part of package trend
##
##    Copyright (C) 2015-2018  Thorsten Pohlert
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
##    This function computes the Multivariate (Multiside) Mann Kendall test.
##  
##
#' @rdname mult.mk.test
#' @title Multivariate (Multisite) Mann-Kendall Test
#' @description
#' Performs a Multivariate (Multisite) Mann-Kendall test.
#' 
#' @param x a time series object of class "ts" 
#' @param alternative the alternative hypothesis, defaults to \code{two.sided}
#'
#' @details
#'  The Mann-Kendall scores are first computed for each variate (side)
#' seperately.
#'
#' \deqn{
#' S = \sum_{k = 1}^{n-1} \sum_{j = k + 1}^n
#' \mathrm{sgn}\left(x_j - x_k\right)}{%
#' S = \sum_{k = 1}^{n-1} \sum_{j = k + 1}^n
#' sgn(x[j] - x[k])}
#' 
#' with \eqn{\mathrm{sgn}}{sgn} the signum function (see \code{\link{sign}}).
#'
#'  The variance - covariance matrix is computed according to
#'  Libiseller and Grimvall (2002).
#'
#' \deqn{
#' \Gamma_{xy} = \frac{1}{3}
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
#'  Finally, the corrected z-statistics for the entire series
#'  is calculated as follows, whereas a continuity correction is employed
#'  for \eqn{n \le 10}{n <= 10}:
#'
#'  \deqn{
#'  z = \frac{\sum_{i=1}^d S_i}{\sqrt{\sum_{j=1}^d \sum_{i=1}^d \Gamma_{ij}}}
#'  }{%
#'  z = sum(S) / sum(\Gamma)
#'  }
#'
#'  where
#' 
#'  \eqn{z} denotes the quantile of the normal distribution
#'  \eqn{S} is the vector of Mann-Kendall scores
#'  for each variate (site) \eqn{1 \le i \le d}{1 <= i <= d} and
#'  \eqn{\Gamma} denotes symmetric variance - covariance matrix.
#'
#' @return
#' An object with class "htest"
#' 
#' \item{data.name}{character string that denotes the input data}
#' \item{p.value}{the p-value for the entire series}
#' \item{statistic}{the z quantile of the standard normal distribution
#' for the entire series}
#' \item{null.value}{the null hypothesis}
#' \item{estimates}{the estimates S and varS for the entire series}
#' \item{alternative}{the alternative hypothesis}
#' \item{method}{character string that denotes the test}
#' \item{cov}{the variance - covariance matrix}
#'
#' @note
#' Ties are not corrected. Current Version is for complete observations only.
#'
#' @references
#' Hipel, K.W. and McLeod, A.I. (1994),
#' \emph{Time Series Modelling of Water Resources and Environmental Systems}. 
#' New York: Elsevier Science.
#'
#' Lettenmeier, D.P. (1988), Multivariate nonparametric tests for trend
#' in water quality. \emph{Water Resources Bulletin} 24, 505--512.
#' 
#' Libiseller, C. and Grimvall, A. (2002),
#' Performance of partial Mann-Kendall tests for trend detection in the
#' presence of covariates. \emph{Environmetrics} 13, 71--84,
#' \doi{10.1002/env.507}.
#'
#' @seealso
#' \code{\link{cor}},
#' \code{\link{cor.test}},
#' \code{\link{mk.test}},
#' \code{\link{smk.test}}
#' 
#' @examples
#' data(hcb)
#' mult.mk.test(hcb)
#'
#' @keywords ts nonparametric multivariate
#' @aliases "Multiside Mann-Kendall Test"
#' @importFrom stats na.fail pnorm
#' @export
mult.mk.test <- function(x, alternative = c("two.sided","greater", "less"))
{
    if(!is.ts(x)){
	stop("'x' must be a time series object")   
    }
    n <- length(x[,1])
    if(n < 3){
	stop("'x' must have at least 3 elements")
    }
    na.fail(x)
    alternative <- match.arg(alternative)
##   
    DNAME <- deparse(substitute(x))

    ## Convert to data frame
    x <- as.data.frame(x)
    ## number of columns
    d <- length(x[1,])
    ## number of rows
    n <- length(x[,1])

    ## Vector of scores
    S <- numeric(d)
    for (i in 1:d){
        S[i] <- .mkScore(x[,i])
    }

    ## Variance / Covariance matrix
    Gamma <- matrix(NA, ncol=d, nrow=d)
    for (i in 2:d)
    {
        for (j in 1:(i-1))
        {
            K <- .K(x[,i], x[,j])
            Ri <- .R(x[,i])
            Rj <- .R(x[,j])
            Gamma[i, j] <- (K + 4 * sum(Ri * Rj) -
                            n * (n + 1) * (n + 1)) / 3
            Gamma[j, i] <- Gamma[i, j]
        }
    }

    ## diagonal
    for (i in 1:d)
    {
        K <- .K(x[,i], x[,i])
        Ri <- .R(x[,i])
        Rj <- .R(x[,i])
        Gamma[i, i] <- (K + 4 * sum(Ri * Rj) -
                            n * (n + 1) * (n + 1)) / 3
    }
    
    ## one vector
    ##one <- cbind(rep(1, d))

    ## Total MK scores and variance
    ##Sg <- t(one) %*% S
    ##varSg <- (t(one) %*% Gamma %*% one)

    ## simpler
    Sg <- sum(S)
    varSg <- sum(Gamma)

    ## Continuity correction
    if (n <= 10)
    {
        sg <- sign(Sg)

        z <- sg * (abs(Sg) - 1) / sqrt(varSg)

        warning("'z' was corrected for continuity")
    } else {
        
        z <- Sg / sqrt(varSg)
    }

    ## p-value
    pval <- switch(alternative,

                   "two.sided" = 2 * min(0.5, pnorm(abs(z), lower.tail=FALSE)),

                   "greater" = pnorm(z, lower.tail=FALSE),

                   "less" = pnorm(z, lower.tail = TRUE)
                   )

    ## finalise
    ans <- list(statistic = c(z = z),
                null.value = c(S = 0),
                estimates = c(S = Sg, varS = varSg),
                p.value = pval,
                alternative = alternative,
                data.name = DNAME,
                cov = Gamma,
                method = "Multivariate Mann-Kendall Trend Test")
    class(ans) <- "htest"
    return(ans)
}
