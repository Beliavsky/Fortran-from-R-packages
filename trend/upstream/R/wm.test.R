##  file wm.test.R part of the package trend
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
#' @title Wallis and Moore Phase-Frequency Test
#'
#' @description
#'  Performes the non-parametric Wallis and Moore phase-frequency
#'  test for testing the H0-hypothesis, whether the series
#'  comprises random data, against the HA-Hypothesis, that the series
#'  is significantly different from randomness (two-sided test).
#'
#' @param x a vector or a time series object of class "ts"
#'
#' @details
#' The test statistic of the phase-frequency test for \eqn{n > 30}
#' is calculated as:
#'
#' \deqn{z = \frac{| h - \frac{2 n - 7}{3}|}{\sqrt{\frac{16 n - 29}{90}}}}{%
#'   z = abs(h - (2 * n - 7) / 3) / sqrt((16 * n - 29)/ 90)}
#'
#' where \eqn{h} denotes the number of phases, whereas the first and the
#' last phase is not accounted. The \eqn{z}-statistic is
#' normally distributed. For \eqn{n \le 30} a continuity correction of
#' \eqn{-0.5} is included in the denominator.
#'
#' @return
#' An object of class "htest"
#'
#' \item{method}{
#'    a character string indicating the chosen test
#'  }
#'  \item{data.name}{
#'   a character string giving the name(s) of the data
#'  }
#'  \item{statistic}{
#'    the Wallis and Moore z-value
#'  }
#'  \item{alternative}{
#'    a character string describing the alternative hypothesis
#'  }
#'  \item{p.value}{
#'    the p-value for the test
#'  }
#'
#' @note
#'  NA values are omitted. Many ties in the series will lead to reject H0
#'  in the present test.
#'
#' @references
#'  L. Sachs (1997), \emph{Angewandte Statistik}. Berlin: Springer.
#'  
#'  C.-D. Schoenwiese (1992), \emph{Praktische Statistik}. Berlin:
#'  Gebr. Borntraeger.
#'
#'  W. A. Wallis and G. H. Moore (1941): A significance test for time series
#'  and other ordered observations. Tech. Rep. 1. National Bureau of
#'  Economic Research. New York.
#'
#' @seealso \code{\link{mk.test}}
#'
#' @examples
#' ## Example from Schoenwiese (1992, p. 113)
#' ## Number of frost days in April at Munich from 1957 to 1968
#' ## z = -0.124, Accept H0
#' frost <- ts(data=c(9,12,4,3,0,4,2,1,4,2,9,7), start=1957)
#' wm.test(frost)
#'
#' ## Example from Sachs (1997, p. 486)
#' ## z = 2.56, Reject H0 on a level of p < 0.05
#' x <- c(5,6,2,3,5,6,4,3,7,8,9,7,5,3,4,7,3,5,6,7,8,9)
#' wm.test(x)
#'
#' wm.test(nottem)
#'
#' @keywords htest nonparametric
#' @importFrom stats pnorm na.omit
#' @export
wm.test <- function(x)
{
    dname <- deparse(substitute(x))
    stopifnot(is.numeric(x))
    x <- na.omit(x)
    n <- length(x)
    if(n < 2) {
        stop("sample size must be greater than 1")
    }
    x.d <- diff(x)
    x.d <- x.d[x.d != 0]
    x.s <- sign(x.d)
    h <- 0
    nn <- length(x.s)
    if (nn > 2) {
        for (i in 2:nn) {
            if (x.s[i] != x.s[i-1]) h <- h + 1
        }
    }
    h <- h - 1
    h <- max(c(0,h))
    if (n > 30) {
        z <- abs(h - (2 * n - 7) / 3) / sqrt((16 * n - 29)/ 90)
    } else {
        z <- (abs(h - (2 * n - 7) / 3) - 0.5) /
            sqrt((16 * n - 29)/ 90)
    }

    p.value <- 2 * min(0.5, pnorm(abs(z), lower.tail = FALSE))
    
    out <- list(statistic = c(z = z),
                p.value = p.value,
                alternative = "The series is significantly different from randomness",
                data.name = dname,
                method = "Wallis and Moore Phase-Frequency test")
    class(out) <- "htest"
    out
}
