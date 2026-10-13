##   file mk.test.R part of package trend
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
#' @title Mann-Kendall Trend Test
#'
#' @description Performs the Mann-Kendall Trend Test
#'
#' @param x a vector of class "numeric" or a time series object of class "ts"
#' @param alternative the alternative hypothesis, defaults to \code{two.sided}
#' @param continuity logical, indicates whether a continuity correction
#' should be applied, defaults to \code{TRUE}.
#'
#' @details
#' The null hypothesis is that the data come from a population with
#' independent realizations and are identically distributed.
#' For the two sided test, the alternative hypothesis
#' is that the data follow a monotonic trend. The Mann-Kendall test
#' statistic is calculated according to:
#'
#' \deqn{
#' S = \sum_{k = 1}^{n-1} \sum_{j = k + 1}^n
#' \mathrm{sgn}\left(x_j - x_k\right)}{%
#' S = \sum_{k = 1}^{n-1} \sum_{j = k + 1}^n
#' sgn(x[j] - x[k])}
#' 
#' with \eqn{\mathrm{sgn}}{sgn} the signum function (see \code{\link{sign}}).
#' 
#' The mean of \eqn{S} is \eqn{\mu = 0}. The variance including the
#' correction term for ties is
#'
#' \deqn{
#' \sigma^2 = \left\{n \left(n-1\right)\left(2n+5\right) -
#' \sum_{j=1}^p t_j\left(t_j - 1\right)\left(2t_j+5\right) \right\} / 18}
#'
#' where \eqn{p} is the number of the tied groups in the data set and
#' \eqn{t_j} is the number of data points in the \eqn{j}-th tied group.
#' The statistic \eqn{S} is approximately normally distributed, with
#'
#' \deqn{z = S / \sigma} 
#'
#' If \code{continuity = TRUE} then a continuity correction will be employed:
#'
#' \deqn{z = \mathrm{sgn}(S) ~ \left(|S| - 1\right) / \sigma}{%
#'   z = sgn(S) * (|S| -1) / \sigma}
#'
#' The statistic \eqn{S} is closely related to Kendall's \eqn{\tau}:
#' \deqn{\tau = S / D}
#'
#' where
#'
#' \deqn{
#' D = \left[\frac{1}{2}n\left(n-1\right)-
#'      \frac{1}{2}\sum_{j=1}^p t_j\left(t_j - 1\right)\right]^{1/2}
#'     \left[\frac{1}{2}n\left(n-1\right) \right]^{1/2}
#' }
#'
#' @return  A list with class "htest"
#' \item{data.name}{character string that denotes the input data}
#' \item{p.value}{the p-value}
#' \item{statistic}{the z quantile of the standard normal distribution}
#' \item{null.value}{the null hypothesis}
#' \item{estimates}{the estimates S, varS and tau}
#' \item{alternative}{the alternative hypothesis}
#' \item{method}{character string that denotes the test}
#' 
#' @references
#' Hipel, K.W. and McLeod, A.I. (1994),
#' \emph{Time Series Modelling of Water Resources and Environmental Systems}. 
#' New York: Elsevier Science.
#'
#' Libiseller, C. and Grimvall, A., (2002), Performance of partial
#' Mann-Kendall tests for trend detection in the presence of covariates.
#' \emph{Environmetrics} 13, 71--84, \doi{10.1002/env.507}.
#'
#' @note Current Version is for complete observations only.
#'
#' @seealso
#' \code{\link{cor.test}},
#' \code{\link[Kendall]{MannKendall}},
#' \code{\link{partial.mk.test}},
#' \code{\link{sens.slope}}
#'
#' @examples
#' data(Nile)
#' mk.test(Nile, continuity = TRUE)
#'
#' ##
#' n <- length(Nile)
#' cor.test(x=(1:n),y=Nile, meth="kendall", continuity = TRUE)
#'
#' @keywords ts
#' @keywords nonparametric
#' @keywords univar
#'
#' @importFrom stats pnorm
#' @export
mk.test <- function(x, alternative = c("two.sided","greater", "less"),
                            continuity = TRUE)
{
    if(!is.numeric(x)){
	stop("'x' must be a numeric vector")
    }
    n <- length(x)
    if(n < 3){
	stop("'x' must have at least 3 elements")
    }
    na.fail(x)
    alternative <- match.arg(alternative)

    S <- .mkScore(x)
    
    ## get ties
    t <- table(x)
    names(t) <- NULL
    varS <- .varmk(t, n)
    D <- .Dfn(t, n)
    tau <- S / D
    
    ## compute z
    if (continuity){
        sg <- sign(S)
        z <- sg * (abs(S) - 1) / sqrt(varS)
    } else {
        z <- S / sqrt(varS)
    }

    ## get the pvalue
    pval <- switch(alternative,
                   "two.sided" = 2 * min(0.5, pnorm(abs(z), lower.tail=FALSE)),
                   "greater" = pnorm(z, lower.tail=FALSE),
                   "less" = pnorm(z, lower.tail=TRUE)
                   )

    ## finalise
    DNAME <- deparse(substitute(x))
    ans <- list(data.name = DNAME, 
                p.value = pval, 
                statistic = c(z = z),
                null.value = c(S = 0),
                parameter = c(n = n),
                estimates = c(S = S, varS = varS, tau = tau),
                alternative = alternative,
                method = "Mann-Kendall trend test",
                pvalg = pval)
    class(ans) <- "htest"
    return(ans)
}
