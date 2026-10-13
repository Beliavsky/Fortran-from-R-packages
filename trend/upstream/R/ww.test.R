##    file www.test.R part of package trend
##
##    Copyright (C) 2016-2018  Thorsten Pohlert
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
#' @title Wald-Wolfowitz Test for Independence and Stationarity
#'
#' @description
#' Performes the non-parametric Wald-Wolfowitz test for
#' independence and stationarity.
#'
#' @param x a vector or a time series object of class "ts"
#'
#' @details
#'  Let \eqn{x_1, x_2, ..., x_n} denote the sampled data, then the test
#'  statistic of the Wald-Wolfowitz test is calculated as:
#'
#'  \deqn{R = \sum_{i=1}^{n-1} x_i x_{i+1} + x_1 x_n}{%
#'    R = sum(x[1:(n-1)] * x[2:n]) + x[1] * x[n]}
#'
#'  The expected value of R is:
#'
#'  \deqn{E(R) = \frac{s_1^2 - s_2}{n - 1}}{%
#'    E(R) = (s1^2 - s2) / (n - 1)}
#'
#'  The expected variance is:
#'
#'  \deqn{V(R) = \frac{s_2^2 - s_4}{n - 1} - E(R)^2 + \frac{s_1^4 -
#'      4 s_1^2 s_2 + 4 s_1 s_3 + s_2^2 - 2 s_4}{(n - 1) (n - 2)}}{%
#'    V(R) = (s2^2 - s4) / (n - 1) - er^2 +
#'    (s1^4 - 4 * s1^2 * s2 + 4 * s1 * s3 + s2^2 - 2 * s4) /
#'    ((n - 1) * (n - 2))}
#'
#'  with:
#'	
#'  \deqn{s_t = \sum_{i=1}^{n} x_i^t, ~~ t = 1, 2, 3, 4}{%
#'    st = sum(x^t), ~~ t = 1, 2, 3, 4}
#'  
#'  For \eqn{n > 10} the test statistic is normally distributed, with:
#'
#'  \deqn{z = \frac{R - E(R)}{\sqrt{V(R)}}}{%
#'    z = (R - E(R)) / V(R)^0.5}
#'
#'  ww.test calculates p-values from the standard normal distribution for
#'  the two-sided case.
#'
#' @return An object of class "htest"
#' \item{method}{
#'    a character string indicating the chosen test
#'  }
#'  \item{data.name}{
#'    a character string giving the name(s) of the data
#'  }
#'  \item{statistic}{
#'    the Wald-Wolfowitz z-value
#'  }
#'  \item{alternative}{
#'    a character string describing the alternative hypothesis
#'  }
#'  \item{p.value}{
#'    the p-value for the test
#'  }
#'
#' @note NA values are omitted.
#'
#' @references
#'  R. K. Rai, A. Upadhyay, C. S. P. Ojha and L. M. Lye (2013),
#'  Statistical analysis of hydro-climatic variables. In:
#'  R. Y. Surampalli, T. C. Zhang, C. S. P. Ojha, B. R. Gurjar,
#'  R. D. Tyagi and C. M. Kao (ed. 2013), \emph{Climate change modelling,
#'    mitigation, and adaptation}. Reston, VA: ASCE. doi =
#'  10.1061/9780784412718.
#'  
#'  A. Wald and J. Wolfowitz (1943), An exact test for randomness in the
#'  non-parametric case based on serial correlation.
#'  \emph{Annual Mathematical Statistics} 14, 378--388.
#'
#'  WMO (2009), \emph{Guide to Hydrological Practices}. Volume II,
#'  Management of Water Resources and Application of Hydrological
#'  Practices, WMO-No. 168.
#'
#' @examples
#' ww.test(nottem)
#' ww.test(Nile)
#'
#' set.seed(200)
#' x <- rnorm(100)
#' ww.test(x)
#'
#' @keywords htest nonparametric
#' @importFrom stats pnorm na.omit
#' @export
ww.test <- function(x)
{
    DNAME <- deparse(substitute(x))
    stopifnot(is.numeric(x))
    x <- na.omit(x)

    n <- length(x)
    if (n < 2) {
        stop("sample size must be greater than 1")
    }
    
# The test statistic r
    r <- sum(x[1:(n-1)] * x[2:n]) + x[1] * x[n]

    s1 <- sum(x^1)
    s2 <- sum(x^2)
    s3 <- sum(x^3)
    s4 <- sum(x^4)

# expected value of r
    er <- (s1^2 - s2) / (n - 1)

# expected variance of r
    vr <- (s2^2 - s4) / (n - 1) - er^2 +
          (s1^4 - 4 * s1^2 * s2 + 4 * s1 * s3 + s2^2 - 2 * s4) /
          ((n - 1) * (n - 2))

# Standardised quantile (approaches normal distribution)	  
    z <- (r - er) / sqrt(vr)
    
# p-value, two-sided case
    pval <- 2 * min(0.5, pnorm(abs(z), lower.tail = FALSE))

    out <- list(statistic = c(z = z),
                p.value= pval,
                alternative = "The series is significantly different from \n independence and stationarity",
                method = "Wald-Wolfowitz test for independence and stationarity",
                data.name = DNAME,
                parameter = c(n = n))
			
    class(out) <- "htest"
    return(out)
}
