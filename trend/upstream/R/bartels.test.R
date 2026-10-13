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
##    This function computes Bartels's test
##    (Rank version of von Neumann's Ratio Test).
#     Changes
#     2017-07-12
#     - Initial writing
#
#' @title Bartels Test for Randomness
#' 
#' @description
#' Performes a rank version of von Neumann's ratio test as proposed by
#' Bartels. The null hypothesis of randomness is tested against the
#' alternative hypothesis
#' 
#' @param x a vector of class "numeric" or a time series object of class "ts"
#'
#' @details
#' In this function, the test is implemented as given by Bartels (1982),
#' where the ranks \eqn{r_1, \ldots, r_n}{r[1], ..., r[n]} of
#' the \eqn{X_i, \ldots, X_n}{X[i], ..., X[n]} are used for the statistic:
#'
#' \deqn{
#' T = \frac{\sum_{i=1}^n (r_i - r_{i+1})^2}{\sum_{i=1}^n (r_i - \bar{r})^2}
#' }{%
#' T = \sum_{i=1}^n (r[i] - r[i+1])^2 / \sum(r[i] - meanr)^2 }
#'
#' As proposed by Bartels (1982), the \eqn{p}-value is calculated
#' for sample sizes in the range of \eqn{(10 \le n < 100)}{(10 <= n < 100)} with
#' the non-standard beta distribution for the range \eqn{0 \le x \le 4}{0 <= x <= 4}
#' with parameters:
#'
#' \deqn{
#'  a = b = \frac{5 n \left( n + 1\right) \left(n - 1\right)^2}
#'          {2 \left(n - 2\right)  \left(5n^2 - 2n - 9\right)} - \frac{1}{2}
#' }{%
#'  a = b = 5 * n * ( n + 1) * (n - 1)^2 /
#'            (2 * ( n - 2) * (5 * n^2 - 2 * n - 9)) - 1/2
#' }
#'
#' For sample sizes \eqn{n \ge 100}{n >= 100} a normal approximation with
#' \eqn{N(2, 20/(5n + 7))} is used for \eqn{p}-value calculation.
#' 
#' @return A list with class "htest"
#'
#' \item{data.name}{character string that denotes the input data}
#' \item{p.value}{the p-value}
#' \item{statistic}{the test statistic}
#' \item{alternative}{the alternative hypothesis}
#' \item{method}{character string that denotes the test}
#' 
#' @note
#'  The current function is for complete observations only.
#'  
#'
#' @references
#' R. Bartels (1982), The Rank Version of von Neumann's Ratio Test
#' for Randomness, \emph{Journal of the American Statistical Association}
#' 77, 40--46.
#' 
#' @seealso
#' \code{\link{ww.test}},
#' \code{\link{wm.test}}
#' 
#' @examples
#' # Example from Schoenwiese (1992, p. 113)
#' ## Number of frost days in April at Munich from 1957 to 1968
#' ## 
#' frost <- ts(data=c(9,12,4,3,0,4,2,1,4,2,9,7), start=1957)
#' bartels.test(frost)
#'
#' ## Example from Sachs (1997, p. 486)
#' x <- c(5,6,2,3,5,6,4,3,7,8,9,7,5,3,4,7,3,5,6,7,8,9)
#' bartels.test(x)
#'
#' ## Example from Bartels (1982, p. 43)
#' x <- c(4, 7, 16, 14, 12, 3, 9, 13, 15, 10, 6, 5, 8, 2, 1, 11, 18, 17)
#' bartels.test(x)
#'  
#' @keywords ts nonparametric univar htest
#' @concept von-Neumann
#' @importFrom extraDistr pnsbeta
#' @importFrom stats pnorm
#' @export
bartels.test <- function(x){    
    na.fail(x)
    n <- length(x)
    if (n < 10){
        stop("'x' must have at least 10 elements")
    }
    DNAME <- deparse(substitute(x))
 
    r <- rank(x)
    RVN <- function(r){
        rmean <- mean(r)
        n <- length(r)
        NM <- sum((r - rmean)^2)
        s <- 0
        for (i in 1:(n-1)){
            s <- s + (r[i] - r[i+1])^2
        }
        return(s / NM)
    }

    T <- RVN(r)
    if(n >= 10 & n < 100){
        ## beta distribution
        a <- 5 * n * ( n + 1) * (n - 1)^2 /
            (2 * ( n - 2) * (5 * n^2 - 2 * n - 9)) - 1/2
        b <- a
        pval <- pnsbeta(T, a, b, min = 0, max = 4)
        
    } else {
        ## normal approximation
        pval <- pnorm(T, 2, sqrt(20/(5 * n + 7)))
    }
    
    retval <- list(statistic = c("RVN" = T),
                   p.value =  pval,
                   data.name= DNAME,
                   alternative="The series is significantly different from randomness",
                   method = "Bartels's test for randomness")
    class(retval) <- c("htest")
    return(retval)
}
