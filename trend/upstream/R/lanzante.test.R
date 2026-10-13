##
##    Copyright (C) 2018 Thorsten Pohlert
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
##    This function computes Lanzante's test.
#     2017-11-13
#     - initial writing
#
#' @title Lanzante's Test for Change Point Detection
#' 
#' @description
#' Performes a non-parametric test after Lanzante in order to test for a
#' shift in the central tendency of a time series. The null hypothesis,
#' no shift, is tested against the alternative, shift.
#' 
#' @param x a vector of class "numeric" or a time series object of class "ts"
#' @param method the test method. Defaults to \code{"wilcox.test"}. 
#' 
#' @details
#' Let \eqn{X} denote a continuous random variable, then the following model
#' with a single shift (change-point) can be proposed:
#'
#' \deqn{
#'   x_i = \left\{
#'       \begin{array}{lcl}
#'        \theta + \epsilon_i, & \qquad & i = 1, \ldots, m \\
#'        \theta + \Delta + \epsilon_i & \qquad & i = m + 1, \ldots, n \\
#'       \end{array} \right.}{%
#'   x[i] = \theta + \epsilon[i] for i = 1, ..., m and x[i] = \theta + \Delta
#'   + \epsilon_i for i = m + 1, ..., n}
#'
#' with \eqn{\theta(\epsilon) = 0}. The null hypothesis, H:\eqn{\Delta = 0}
#' is tested against the alternative A:\eqn{\Delta \ne 0}{\delta != 0}.
#'
#' First, the data are transformed into increasing ranks
#' and for each time-step the adjusted rank sum is computed:
#'
#' \deqn{U_k = 2 \sum_{i=1}^k r_i - k \left(n + 1\right) \qquad k = 1, \ldots, n}{%
#' U[k] = 2 * \sum r_i - k (n + 1)   k = 1, ..., n}
#'
#' The probable change point is located at the absolute maximum
#' of the statistic:
#' 
#' \deqn{m = k(\max |U_k|)}{m = k(max |P[k]|)}.
#'
#' For \code{method = "wilcox.test"} the Wilcoxon-Mann-Whitney two-sample
#' test is performed, using \eqn{m} to split the series. Otherwise,
#' the robust rank-order distributional test (\code{\link{rrod.test}} is
#' performed.
#'  
#' @return A list with class "htest" and "cptest".
#' 
#' @references
#' Lanzante, J. R. (1996), Resistant, robust and non-parametric
#' techniques for the analysis of climate data: Theory and examples,
#' including applications to historical radiosonde station data,
#' \emph{Int. J. Clim.}, 16, 1197--1226.
#' 
#' @seealso
#' \code{\link{pettitt.test}}
#' 
#' @examples
#' data(maxau) ; plot(maxau[,"s"])
#' s.res <- lanzante.test(maxau[,"s"])
#' n <- s.res$nobs
#' i <- s.res$estimate
#' s.1 <- mean(maxau[1:i,"s"])
#' s.2 <- mean(maxau[(i+1):n,"s"])
#' s <- ts(c(rep(s.1,i), rep(s.2,(n-i))))
#' tsp(s) <- tsp(maxau[,"s"])
#' lines(s, lty=2)
#' print(s.res)
#'
#'
#' data(PagesData) ; lanzante.test(PagesData)
#' @keywords ts nonparametric univar htest
#' @importFrom stats wilcox.test
#' @export
lanzante.test <-function(x, method = c("wilcox.test", "rrod.test"))
{
    na.fail(x)
    n <- length(x)
    method <- match.arg(method)
    DNAME <- deparse(substitute(x))
    k <- 1:n
    r <- rank(x)
    Uk <- sapply(k, function(x) 2 * sum(r[1:x]) - x*(n+1))
    Uka <- abs(Uk)
    U <- max(Uka)
    K <- k[Uka == U]

    if (is.ts(x)){
        fr <- frequency(x)
        st <- start(x)
        ed <- end(x)
        Uk <- ts(Uk, start=st, end = ed, frequency= fr)
    }
    attr(Uk, 'nm') <- "Uk"
    n.x <- K
    n.y <- n - n.x
    if(n.y < 1L)
        stop("not enough 'y' observations")

    METHOD <- "Lanzante's procedure for single change-point detection"
    
    if (method == "wilcox.test") {
        out <- do.call("wilcox.test", list(x = x[1:K],
                                           y = x[(K+1):n],
                                           alternative = "two.sided"))
        METHOD <- c(METHOD, " with Wilcoxon-Mann-Whitney Test") 
     } else {
         out <- do.call("rrod.test", list(x = x[1:K],
                                          y = x[(K+1):n],
                                          alternative = "two.sided"))
         METHOD <- c(METHOD, " with Robust Rank-Order Distributional Test")
     }
    
    names(K) <- "probable change point at time K"
    retval <- list(nobs = n, 
                   statistic = out$statistic,
                   estimate= K,
                   p.value =  out$p.value,
                   data.name= DNAME,
                   alternative="two.sided",
                   data = Uk,
                   method = METHOD)
    class(retval) <- c("htest", "cptest")
    return(retval)
}
