##    file snh.test.R part of package trend
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
#     2017-07-06
#     - initial writing
#
#' @title Standard Normal Homogeinity Test (SNHT) for Change-Point Detection
#' 
#' @description
#' Performes the Standard Normal Homogeinity Test (SNHT)
#' for change-point detection of a normal variate.
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
#' The test statistic for the SNHT test is calculated as follows:
#'
#' \deqn{T_k = k z_1^2 + \left(n - k\right) z_2^2 \qquad (1 \le k < n)}{%
#' T[k] = k z[1]^2 + (n - k) z[2]^2    (1 <= k < n)}
#'
#' where
#' 
#' \deqn{
#' \begin{array}{l l}
#' z_1 = \frac{1}{k} \sum_{i=1}^k \frac{x_i - \bar{x}}{\sigma} &
#' z_2 = \frac{1}{n-k} \sum_{i=k+1}^n \frac{x_i - \bar{x}}{\sigma}. \\
#' \end{array}}{%
#' z[1] = 1 / k * \sum((x[1:k] - xmean) / \sigma) and
#' z[2] = (n - k) * (1 / (n - k) * \sum((x[(k+1):n] - xmean) / \sigma)).
#' }
#' 
#' The critical value is:
#' \deqn{T = \max T_k.}{T = \max(T[k]).} 
#'
#' The \code{p.value} is estimated with a Monte Carlo simulation
#' using \code{m} replicates.
#' 
#' Critical values based on \eqn{m = 1,000,000} Monte Carlo simulations
#' are tabulated for \eqn{T} by Khaliq and Ouarda (2007).
#' 
#' @return A list with class "htest" and "cptest"
#' \item{data.name}{character string that denotes the input data}
#' \item{p.value}{the p-value}
#' \item{statistic}{the test statistic}
#' \item{null.value}{the null hypothesis}
#' \item{estimates}{the time of the probable change point}
#' \item{alternative}{the alternative hypothesis}
#' \item{method}{character string that denotes the test}
#' \item{data}{numeric vector of Tk for plotting}
#' 
#' @note
#'  The current function is for complete observations only.
#'
#' @references
#' H. Alexandersson (1986), A homogeneity test applied to precipitation data,
#' \emph{Journal of Climatology} 6, 661--675.
#' 
#' M. N. Khaliq, T. B. M. J. Ouarda (2007), On the critical values of the
#' standard normal homogeneity test (SNHT),
#' \emph{International Journal of Climatology} 27, 681--687.
#' 
#' G. Verstraeten, J. Poesen, G. Demaree, C. Salles (2006),
#' Long-term (105 years) variability in rain erosivity as derived from 10-min
#' rainfall depth data for Ukkel (Brussels, Belgium):
#' Implications for assessing soil erosion rates.
#' \emph{Journal of Geophysical Research} 111, D22109.
#' 
#' @seealso
#' \code{\link[strucchange]{efp}}
#' \code{\link[strucchange]{sctest.efp}}
#' 
#' @examples
#' data(Nile)
#' (out <- snh.test(Nile))
#' plot(out)
#'
#' data(PagesData) ; snh.test(PagesData)
#'  
#' @keywords ts univar htest
#' @importFrom stats sd
#' @useDynLib 'trend', .registration = TRUE
#' @export
snh.test <- function(x, m = 20000){
    if(!is.numeric(x)){
        stop("'x' must be a numeric vector")
    }
    na.fail(x)
    DNAME <- deparse(substitute(x))
    
    xmean <- mean(x)
    sigma <- sd(x)
    n <- length(x)
    k <- 1:(n-1)


    Tk <- sapply(k,
                 function(k)
                 {
                     k * (1 / k * sum((x[1:k] - xmean) / sigma))^2 +
                         (n - k) * (1 / (n - k) * sum((x[(k+1):n] - xmean)
                                                      / sigma))^2
                 }
                 )
           
    T <- max(Tk)
    K <- k[Tk == T]

        # standardised value
    
    if (is.ts(x)){
        fr <- frequency(x)
        st <- start(x)
        ed <- end(x)
        Tk <- ts(Tk, start=st, end = ed, frequency= fr)
    }
    

    PVAL <- .Fortran("mcsnht",
                     stat = as.double(T),
                     n = as.integer(n),
                     m = as.integer(m),
                     pval = double(1))$pval

    names(K) <- "probable change point at time K"
    attr(Tk, 'nm') <- "Tk"
    retval <- list(statistic = c(T = T),
                   parameter = c(n = n),
                   null.value = c(delta = 0),
                   estimate= K,
                   p.value =  PVAL,
                   data.name= DNAME,
                   alternative="two.sided",
                   data = Tk,
                   method = "Standard Normal Homogeneity Test (SNHT)")
    class(retval) <- c("htest", "cptest")
    return(retval)
}
