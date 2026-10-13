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
##    This function computes the partial correlation trend test.
##
#' @title Partial Correlation Trend Test
#' @description
#' Performs a partial correlation trend test with either Pearson's or
#' Spearman's correlation coefficients (\eqn{r(tx.z)}).
#'
#' @param x a "vector" or "ts" object that contains the variable,
#' which is tested for trend (i.e. correlated with time)
#' @param z a "vector" or "ts" object that contains the co-variate, which
#' will be partialled out
#' @param method a character string indicating which correlation coefficient
#'    is to be computed. One of "pearson" (default) or "spearman",
#'    can be abbreviated.
#' 
#' @details
#'  This function performs a partial correlation trend test using either
#'  the "pearson" correlation coefficient, or the "spearman" rank
#'  correlation coefficient (Hipel and McLoed (1994), p. 882).
#'  The partial correlation coefficient
#'  for the response variable "x" with time "t",
#'  when the effect of the explanatory variable "z" is partialled out,
#'  is defined as:
#'  \deqn{
#'  r_{tx.z} = \frac{r_{tx} - r_{tz}~r_{xz}}
#'  {\sqrt{1 - r_{tz}^2} ~ \sqrt{1-r_{xz}^2}}
#'  }{%
#'  r_{tx.z} = (r_{tx} - r_{tz}~r_{xz}) /
#'  (sqrt(1 - r_{tz}^2) sqrt(1-r_{xz}^2))}
#'
#' The H0: \eqn{r_{tx.z} = 0}{r(tx.z) = 0} (i.e. no trend for "x", when
#' effect of "z" is partialled out) is tested against the
#' alternate Hypothesis, that there is a trend for "x", when the effect of
#' "z" is partialled out.
#'
#'  The partial correlation coefficient is tested for significance with
#'  the student t distribution on \eqn{df = n - 2} degree of freedom.
#'
#' @return An object of class "htest"
#' \item{method}{
#'    a character string indicating the chosen test
#'  }
#'  \item{data.name}{
#'    a character string giving the name(s) of the data
#'  }
#'  \item{statistic}{
#'    	the value of the test statistic
#'  }
#'  \item{estimate}{
#'    the partial correlation coefficient \eqn{r(tx.z)}
#'  }
#'  \item{parameter}{
#'    the degrees of freedom of the test statistic in the case
#' that it follows a t distribution
#'  }
#'  \item{alternative}{
#'    a character string describing the alternative hypothesis
#'  }
#'  \item{p.value}{
#'    the p-value of the test
#'  }
#' \item{null.value}{The value of the null hypothesis}
#'
#' @references
#' Hipel, K.W. and McLeod, A.I. (1994),
#' \emph{Time Series Modelling of Water Resources and Environmental Systems}. 
#' New York: Elsevier Science.
#'
#' Bahrenberg, G., Giese, E. and Nipper, J., (1992): Statistische Methoden
#' in der Geographie, Band 2 Multivariate Statistik, Teubner, Stuttgart.
#'
#' @note Current Version is for complete observations only.
#'
#' @seealso
#'  \code{\link{cor}},
#'  \code{\link{cor.test}},
#'  \code{\link[psych]{partial.r}},
#'  \code{\link{partial.mk.test}},
#'
#' @examples
#'data(maxau)
#'a <- tsp(maxau) ; tt <- a[1]:a[2]
#'s <- maxau[,"s"] ; Q <- maxau[,"Q"]
#'maxau.df <- data.frame(Year = tt, s =s, Q = Q)
#'plot(maxau.df)
#'
#'partial.cor.trend.test(s,Q, method="pearson")
#'partial.cor.trend.test(s,Q, method="spearman")
#'
#' @keywords ts nonparametric multivariate
#' @importFrom stats cor
#' @importFrom stats pt
#' @export
partial.cor.trend.test <- function(x, z, method = c("pearson", "spearman"))
{
    method <- match.arg(method)
    na.fail(x)
    na.fail(z)
    n1 <- length(x)
    n2 <- length(z)
    if (n1 != n2) {
        stop("Error: time-series x and y must be of same length")
    }
    n <- n1
    dat <- data.frame(t=(1:n),x, z)
    rho <- cor(dat, method=method)
    a <- rho[1,2] - rho[2,3] * rho[1,3]
    b <- sqrt((1 - rho[2,3]^2)) * sqrt((1 - rho[1,3]^2))
    if (a == 0 & b == 0){
        rho.xt.z <- 0
    }
    else {
        rho.xt.z <- a / b
    }
    T <- (sqrt(n -2) * rho.xt.z) / sqrt((1 - rho.xt.z^2))
    pvalue <- 2 * (1 - pt(abs(T), df=(n-2)))

    res <- list(statistic=NULL,
                parameter = NULL,
                estimate = NULL,
                p.value =NULL,
                statistic =NULL,
                null.value = c("rho" = 0),
                alternative = "two.sided",
                method = NULL,
                data.name = NULL,
                cor)

    DNAMEX <- deparse(substitute(x))
    DNAMEY <- deparse(substitute(z))
    partrval <- paste("r(t",DNAMEX,".", DNAMEY,")", sep="")
    res$data.name <- paste("t AND ",DNAMEX," . ",DNAMEY)
    names(rho.xt.z) <- partrval
    res$estimate <- rho.xt.z
    names(T) <- "t"
    res$statistic <- T
    df <- n -2
    names(df) <- "df"
    res$parameter <- df
    names(pvalue) <- "p-value"
    res$p.value <- pvalue
    res$cor <- rho
    dimnames(res$cor)[[1]] <- c("t", DNAMEX, DNAMEY)
    dimnames(res$cor)[[2]] <- c("t", DNAMEX, DNAMEY)

    if(method == "pearson"){
        res$method <- "Pearson's Partial Correlation Trend Test"
    }
    else {
        res$method <- "Spearman's Partial Correlation Trend Test"
    }
    class(res) <- "htest"
    return(res)
}
