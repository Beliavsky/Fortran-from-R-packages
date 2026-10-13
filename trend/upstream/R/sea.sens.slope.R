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
##    This function computes seasonal Sens's Slope.
##
#' @title Seasonal Sen's Slope
#' @description Computes seasonal Sen's slope  for linear rate of change
#'
#' @param x a time series object of class "ts"
#'
#' @details
#' Acccording to Hirsch et al. (1982) the seasonal Sen's slope
#' is calculated as follows:
#'
#' \deqn{
#' d_{ijk} = \frac{x_{ij} - x_{ik}}{j - k}
#' }{%
#' d[ijk] = (x[ij] - x[ik]) / (j - k)
#' }
#'
#' for each \eqn{(x_{ij}, x_{ik})}{(x[ij], x[ik])}
#' pair \eqn{i = 1, \ldots, m},
#' where \eqn{(1 \leq k < j \leq n_i} and \eqn{n_i}{n[i]} is
#' the number of known values in the \eqn{i-th} season.
#' The seasonal slope estimator is the median of
#' the \eqn{d_{ijk}}{d[ijk]} values.
#'
#' @return numeric, Seasonal Sen's slope.
#'
#' @references
#' Hipel, K.W. and McLeod, A.I. (1994),
#' \emph{Time Series Modelling of Water Resources and Environmental Systems}. 
#' New York: Elsevier Science.
#'
#' Hirsch, R., J. Slack, R. Smith (1982), T
#' echniques of Trend Analysis for Monthly Water Quality Data.
#' \emph{Water Resources Research} 18, 107-121.
#'
#' Sen, P.K. (1968), Estimates of the regression coefficient based on
#' Kendall's tau, \emph{Journal of the American Statistical Association} 63,
#' 1379--1389.
#'
#' @note Current Version is for complete observations only.
#'
#' @seealso \code{\link{smk.test}},
#'
#' @examples
#' sea.sens.slope(nottem)
#'
#' @keywords ts nonparametric univar 
#'
#' @importFrom stats na.fail is.ts frequency median
#' 
#' @export
sea.sens.slope <- function(x)
{
        na.fail(x)
        if (!is.ts(x)){
            stop("error: input must be ts object")
        }
        fr <- frequency(x)
        if (fr < 2){
            stop("error: minimum of two seasons required.")
        }

 
        ## diff function
        .d <- function(x)
        {
           n <- length(x)
           i <- 0
           d <- rep(NA, n * (n-1)/2)
           for (k in 1:(n-1)) {
               for (j in (k+1):n){
                   i <- i + 1
                   d[i] <- (x[j] - x[k]) / ( j - k)
               }
           }
           return(d)
        }
        
        ## loop over seasons
        p <- length(x) / fr
        d <- matrix(NA, ncol=fr, nrow = p * (p - 1) / 2)
        
        for (i in 1:fr)
        {
            dat <- x[cycle(x) == i]
            d[,i] <- .d(dat)
        }

        slp <- median(d)
        return(slp)
    }
