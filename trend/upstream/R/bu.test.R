##    file bu.test.R part of package trend
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
##   
#     Changes
#     2017-07-01
#     - initial writing
#
#' @title Buishand U Test for Change-Point Detection
#' 
#' @description
#' Performes the Buishand U test for change-point detection
#' of a normal variate.
#' 
#' @param x a vector of class "numeric" or a time series object of class "ts"
#' @param m numeric, number of Monte-Carlo replicates, defaults to 20000
#' 
#' @details
#' Let \eqn{X} denote a normal random variate, then the following model
#' with a single shift (change-point) can be proposed:
#'
#' \deqn{
#'   x_i = \left\{
#'       \begin{array}{lcl}
#'        \mu + \epsilon_i, & \qquad & i = 1, \ldots, m \\
#'        \mu + \Delta + \epsilon_i & \qquad & i = m + 1, \ldots, n \\
#'       \end{array} \right.}{%
#'   x[i] = \mu + \epsilon[i] for i = 1, ..., m and x[i] = \mu + \delta
#'   + \epsilon_i for i = m + 1, ..., n}
#'
#' with \eqn{\epsilon \approx N(0,\sigma)}. The null hypothesis \eqn{\Delta = 0}
#' is tested against the alternative \eqn{\Delta \ne 0}{\delta != 0}.
#' 
#' In the Buishand U test, the rescaled adjusted partial sums
#' are calculated as
#'
#' \deqn{S_k = \sum_{i=1}^k \left(x_i - \bar{x}\right) \qquad (1 \le i \le n)}{%
#' S[k] = \sum (x[i] - xmean)   (1, <= i <= n)}
#'
#' The sample standard deviation is
#' \deqn{
#'    D_x = \sqrt{n^{-1} \sum_{i=1}^n \left(x_i - \bar{x}\right)}}{%
#'    Dx = sqrt(n^(-1) \sum(x - \mu))}
#' 
#' The test statistic is calculated as:
#' \deqn{U = \left[n \left(n + 1 \right) \right]^{-1} \sum_{k=1}^{n-1}
#' \left(S_k / D_x \right)^2 }{%
#' U = 1 / [n * (n + 1)] * \sum_{k=1}^{n-1} (S[k] - Dx)^2}. 
#'
#' The \code{p.value} is estimated with a Monte Carlo simulation
#' using \code{m} replicates.
#' 
#' Critical values based on \eqn{m = 19999} Monte Carlo simulations
#' are tabulated for \eqn{U} by Buishand (1982, 1984).
#' 
#' @return A list with class "htest" and "cptest"
#' \item{data.name}{character string that denotes the input data}
#' \item{p.value}{the p-value}
#' \item{statistic}{the test statistic}
#' \item{null.value}{the null hypothesis}
#' \item{estimates}{the time of the probable change point}
#' \item{alternative}{the alternative hypothesis}
#' \item{method}{character string that denotes the test}
#' \item{data}{numeric vector of Sk for plotting}
#' 
#' @note
#'  The current function is for complete observations only.
#'
#' @references
#' T. A. Buishand (1982), Some Methods for Testing the Homogeneity
#' of Rainfall Records, \emph{Journal of Hydrology} 58, 11--27.
#'
#' T. A. Buishand (1984), Tests for Detecting a Shift in the Mean of
#' Hydrological Time Series, \emph{Journal of Hydrology} 73, 51--69.
#' 
#' @seealso
#' \code{\link[strucchange]{efp}}
#' \code{\link[strucchange]{sctest.efp}}
#' 
#' @examples
#' data(Nile)
#' (out <- bu.test(Nile))
#' plot(out)
#'
#' data(PagesData)
#' bu.test(PagesData)
#'  
#' @keywords ts univar htest
#' @useDynLib 'trend', .registration = TRUE
#' @importFrom stats na.fail is.ts start end frequency ts
#' @export
bu.test <- function(x, m = 20000){    
    na.fail(x)
    DNAME <- deparse(substitute(x))
    
    xmean <- mean(x)
    n <- length(x)
    k <- 1:n
    Sk <- sapply(k, function(i) sum(x[1:i] - xmean))
    sigma <- sd(x)
    U <- 1 / (n * ( n + 1)) * sum((Sk[1:(n-1)] / sigma)^2)
    Ska <- abs(Sk)
    S <- max(Ska)
    K <- k[Ska == S]

    ## standardised value
    Skk <- (Sk / sigma)
    if (is.ts(x)){
        fr <- frequency(x)
        st <- start(x)
        ed <- end(x)
        Skk <- ts(Sk, start=st, end = ed, frequency= fr)
    }
        
        
    PVAL <- .Fortran("mcbu",
                     stat = as.double(U),
                     n = as.integer(n),
                     m = as.integer (m),
                     pval = double(1))$pval
    attr(Skk, 'nm') <- "Sk**"
    names(K) <- "probable change point at time K"
    retval <- list(statistic = c(U =U),
                   parameter = c(n = n),
                   null.value = c(delta = 0),
                   estimate= K,
                   p.value =  PVAL,
                   data.name= DNAME,
                   alternative="two.sided",
                   data = Skk,
                   method = "Buishand U test")
    class(retval) <- c("htest", "cptest")
    return(retval)
}
