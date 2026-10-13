##  file cs.test.R part of the package trend
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
#' @title Cox and Stuart Trend Test
#'
#' @description Performes the non-parametric Cox and Stuart trend test
#'
#' @param x a vector or a time series object of class "ts"
#'
#' @details
#' First, the series is devided by three. It is compared, whether the
#' data of the first third of the series are larger or smaller than the
#' data of the last third of the series.
#' The test statistic of the Cox-Stuart trend test for \eqn{n > 30}
#' is calculated as:
#'
#' \deqn{z = \frac{|S - \frac{n}{6}|}{\sqrt{\frac{n}{12}}}}{%
#' z = abs(S - n / 6) / sqrt(n / 12)}
#'
#'  where \eqn{S} denotes the maximum of the number of signs, i.e.
#' \eqn{+} or \eqn{-}, respectively. The \eqn{z}-statistic is
#' normally distributed. For \eqn{n \le 30} a continuity correction of
#' \eqn{-0.5} is included in the denominator.
#'
#' @return
#' An object of class "htest"
#'  \item{method}{a character string indicating the chosen test}
#'  \item{data.name}{a character string giving the name(s) of the data}
#' \item{statistic}{the Cox-Stuart z-value}
#' \item{alternative}{a character string describing the alternative hypothesis}
#' \item{p.value}{the p-value for the test}
#'
#' @note NA values are omitted. Many ties in the series will lead
#' to reject H0 in the present test.
#'
#' @references
#' L. Sachs (1997), \emph{Angewandte Statistik}. Berlin: Springer.
#'
#' C.-D. Schoenwiese (1992), \emph{Praktische Statistik}. Berlin:
#' Gebr. Borntraeger.
#'
#' D. R. Cox and A. Stuart (1955), Quick sign tests for trend in location
#' and dispersion. \emph{Biometrika} 42, 80-95.
#'
#' @seealso \code{\link{mk.test}}
#'
#' @examples
#' ## Example from Schoenwiese (1992, p. 114)
#' ## Number of frost days in April at Munich from 1957 to 1968
#' ## z = -0.5, Accept H0
#' frost <- ts(data=c(9,12,4,3,0,4,2,1,4,2,9,7), start=1957)
#' cs.test(frost)
#'
#' ## Example from Sachs (1997, p. 486-487)
#' ## z ~ 2.1, Reject H0 on a level of p = 0.0357
#' x <- c(5,6,2,3,5,6,4,3,7,8,9,7,5,3,4,7,3,5,6,7,8,9)
#' cs.test(x)
#'
#' cs.test(Nile)
#'
#' @keywords htest nonparametric
#' @importFrom stats pnorm na.omit
#' @export
cs.test <- function(x)
{
    dname <- deparse(substitute(x))
    stopifnot(is.numeric(x))
    x <- na.omit(x)
    n <- length(x)
    if(n < 3) {
        stop("sample size must be greater than 2")
    }
    l <- ceiling(n/3)
    u <- x[(n-l+1):n] - x[1:l]
    u <- u[u != 0]
    if (length(u) == 0){
        stop("entire sample contains identical values") 
    }
    sgn <- sign(u)
    S <- table(sgn)
    zstat <- function(S, n) {
        if (n > 30) {
            out <- abs(S - n / 6) / sqrt(n / 12)
            return(out)
        } else {
            out <- (abs(S - n / 6) - 0.5) / sqrt(n / 12)
            return(out)
        }
    }
    S <- max(S)
    z <- zstat(S, n)

    ## two sided
    p.value <- 2 * min(0.5, pnorm(abs(z), lower.tail = FALSE))

    out <- list(statistic = c(z = z),
                parameter = c(n = n),
                p.value = p.value,
                alternative = "monotonic trend",
                data.name = dname,
                method = "Cox and Stuart Trend test")
    class(out) <- "htest"
    out
}
