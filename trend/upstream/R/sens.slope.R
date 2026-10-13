##    file sens.slope.R part of package trend
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
##    This function computes Sens's Slope.
##
#' @title Sen's slope
#' @description
#' Computes Sen's slope  for linear rate of change and corresponding
#' confidence intervalls
#'
#' @param x numeric vector or a time series object of class "ts"
#' @param conf.level numeric, the level of significance
#'
#' @details
#' This test computes both the slope (i.e. linear rate of change) and
#'  confidence levels according to Sen's method. First, a set of linear slopes is
#'  calculated as follows:
#'  \deqn{d_{k} = \frac{x_j - x_i}{j - i}}{%
#'    d(k) = (x(j) - x(i)) / (j - i)}
#'
#'  for \eqn{\left(1 \le i < j \le n \right)}{(1 <= i < j <= n)}, where d
#'  is the slope, x denotes the variable, n is the number of data, and i,
#'  j are indices.
#'
#'  Sen's slope is then calculated as the median from all slopes:
#'  \eqn{b_{Sen} = \textnormal{median}(d_k)}{b = Median(d(k))}.
#'
#'  This function also computes the upper and lower confidence limits for
#'  sens slope.
#'
#' @return
#' A list of class "htest".
#' 
#' \item{estimates}{numeric, Sen's slope}
#' \item{data.name}{character string that denotes the input data}
#' \item{p.value}{the p-value}
#' \item{statistic}{the z quantile of the standard normal distribution}
#' \item{null.value}{the null hypothesis}
#' \item{conf.int}{upper and lower confidence limit}
#' \item{alternative}{the alternative hypothesis}
#' \item{method}{character string that denotes the test}
#'
#' @references
#' Hipel, K.W. and McLeod, A.I. (1994),
#' \emph{Time Series Modelling of Water Resources and Environmental Systems}. 
#' New York: Elsevier Science.
#'
#'    Sen, P.K. (1968), Estimates of the regression coefficient based on
#'    Kendall's tau, \emph{Journal of the American Statistical Association} 63,
#'    1379--1389.
#'
#' @note Current Version is for complete observations only.
#'
#' @examples
#' data(maxau)
#' sens.slope(maxau[,"s"])
#' mk.test(maxau[,"s"])
#'
#' @keywords ts nonparametric univar
#' 
#' @importFrom stats na.fail median qnorm pnorm
#'
#' @export
sens.slope <-function(x, conf.level = 0.95)
{
    if(!is.numeric(x)){
        stop("'x' must be a numeric vector")
    }
    na.fail(x)
    n <- length(x)

    # get ties
    t <- table(x)
    names(t) <- NULL
    varS <- .varmk(t, n)
 

    k <- 0
    d <- rep(NA, n * (n-1)/2)
    for (i in 1:(n-1)) {
        for (j in (i+1):n){
            k <- k + 1
            d[k] <- (x[j] - x[i]) / ( j - i)
        }
    }
    b.sen <- median(d, na.rm=TRUE)
    
    C <- qnorm(1 - (1 - conf.level)/2) * sqrt(varS)
    rank.up <- round((k + C) / 2 + 1)
    rank.lo <- round((k - C) / 2)
    rank.d <- sort(d)
    lo <- rank.d[rank.lo]
    up <- rank.d[rank.up]

    S <- .mkScore(x)
    ## continuity correction
    sg <- sign(S)
    z <- sg * (abs(S) - 1) / sqrt(varS)
    
    pval <- 2 * min(0.5, pnorm(abs(z), lower.tail=FALSE))
            
    cint <- c(lo, up)
    attr(cint, "conf.level") <- conf.level

    ans <- list(estimates = c("Sen's slope" = b.sen),
                statistic = c(z = z),
                p.value = pval,
                null.value = c(z = 0),
                alternative = "two.sided",
                data.name = deparse(substitute(x)),
                method = "Sen's slope",
                parameter = c(n = n),
                conf.int = cint)
    class(ans) <- "htest"
    return(ans)
}


