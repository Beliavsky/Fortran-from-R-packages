##
##    Copyright (C) 2015-2018 Thorsten Pohlert
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
##    This function computes Pettitt's test.
#     Changes
#     2017-07-01
#     - Rd ported to roxygen2
#     - Function was re-written
#
#' @title Pettitt's Test for Change-Point Detection
#' 
#' @description
#' Performes a non-parametric test after Pettitt in order to test for a
#' shift in the central tendency of a time series. The H0-hypothesis,
#' no change, is tested against the HA-Hypothesis, change.
#' 
#' @param x a vector of class "numeric" or a time series object of class "ts"
#'
#' @details
#' In this function, the test is implemented as given by Verstraeten et. al.
#' (2006), where the ranks \eqn{r_1, \ldots, r_n}{r[1], ..., r[n]} of
#' the \eqn{X_i, \ldots, X_n}{X[i], ..., X[n]} are used for the statistic:
#'
#' \deqn{U_k = 2 \sum_{i=1}^k r_i - k \left(n + 1\right) \qquad k = 1, \ldots, n}{%
#' U[k] = 2 * \sum r_i - k (n + 1)   k = 1, ..., n}
#'
#' The test statistic is the maximum of the absolute value of the vector:
#' 
#' \deqn{\hat{U} = \max |U_k|}{U* = max |P[k]|}.
#'
#' The probable change-point \eqn{K} is located where \eqn{\hat{U}}{U*} has its maximum.
#' The approximate probability for a two-sided test is calculated
#' according to 
#'
#' \deqn{p = 2 \exp^{-6K^2 / (T^3 + T^2)}}{%
#'   p =2 exp[-6K^2 / (T^3 + T^2)].}
#' 
#' @return A list with class "htest" and "cptest"
#'
#' @note
#'  The current function is for complete observations only.
#'  The approximate probability is good for \eqn{p \le 0.5}.
#'
#' @references
#' CHR (ed., 2010), Das Abflussregime des Rheins und seiner Nebenfluesse im
#' 20. Jahrhundert, Report no I-22 of the CHR, p. 172.
#'
#' Pettitt, A. N. (1979), A non-parametric approach to the change point
#' problem. \emph{Journal of the Royal Statistical Society Series C}, Applied
#' Statistics 28, 126-135.
#'
#' G. Verstraeten, J. Poesen, G. Demaree, C. Salles (2006),
#' Long-term (105 years) variability in rain erosivity as derived from 10-min
#' rainfall depth data for Ukkel (Brussels, Belgium):
#' Implications for assessing soil erosion rates.
#' \emph{Journal of Geophysical Research} 111, D22109.
#' @seealso
#' \code{\link[strucchange]{efp}}
#' \code{\link[strucchange]{sctest.efp}}
#' 
#' @examples
#' data(maxau) ; plot(maxau[,"s"])
#' s.res <- pettitt.test(maxau[,"s"])
#' n <- s.res$nobs
#' i <- s.res$estimate
#' s.1 <- mean(maxau[1:i,"s"])
#'s.2 <- mean(maxau[(i+1):n,"s"])
#' s <- ts(c(rep(s.1,i), rep(s.2,(n-i))))
#' tsp(s) <- tsp(maxau[,"s"])
#' lines(s, lty=2)
#' print(s.res)
#'
#'
#' data(PagesData) ; pettitt.test(PagesData)
#'  
#' @keywords ts nonparametric univar htest
#' @export
pettitt.test <- function(x){    
    na.fail(x)
    n <- length(x)
    DNAME <- deparse(substitute(x))
    k <- 1:n
    r <- rank(x)
    Uk <- sapply(k, function(x) 2 * sum(r[1:x]) - x*(n+1))
    Uka <- abs(Uk)
    U <- max(Uka)
    K <- k[Uka == U]
    ## ensure not to exceed 1
    pval <- min(1, 2.0 * exp(( -6.0 * U^2) / (n^3 + n^2)))

    if (is.ts(x)){
        fr <- frequency(x)
        st <- start(x)
        ed <- end(x)
        Uk <- ts(Uk, start=st, end = ed, frequency= fr)
    }

    attr(Uk, 'nm') <- "Uk"
    
    names(K) <- "probable change point at time K"
    retval <- list(nobs = n, 
                   statistic = c("U*" = U),
                   estimate= K,
                   p.value =  pval,
                   data.name= DNAME,
                   alternative="two.sided",
                   data = Uk,
                   method = "Pettitt's test for single change-point detection")
    class(retval) <- c("htest", "cptest")
    return(retval)
}
