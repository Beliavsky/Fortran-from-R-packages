##    file csmk.test.R part of package trend
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
##    This function computes the Correlated Seasonal Mann Kendall test.
##    It calls mult.mk.test()
##
#' @title Correlated Seasonal Mann-Kendall Test
#' @description
#' Performs a Seasonal Mann-Kendall test under the presence of correlated
#'  seasons.
#' 
#' @param x a time series object with class \code{ts} comprising >= 2 seasons;
#' \code{NA} values are not allowed
#' @param alternative the alternative hypothesis, defaults to \code{two.sided}
#'
#' @details{
#' 
#'  The Mann-Kendall scores are first computed for each season seperately.
#'  The variance - covariance matrix is computed according to Libiseller and Grimvall (2002).
#'  Finally the corrected Z-statistics for the entire series
#'  is calculated as follows, whereas a continuity correction is employed
#'  for \eqn{n \le 10}{n <= 10}:
#'
#'  \deqn{
#'    z = \frac{\mathbf{1}^T \mathbf{S}}
#'    {\sqrt{\mathbf{1}^T \mathbf{\Gamma}~\mathbf{1}}}
#'    }{%
#'    z = 1^T S / sqrt(1^T \Gamma 1)
#'    }
#'
#'  where
#' 
#'  \eqn{z} denotes the quantile of the normal distribution, 1 indicates a vector
#'  with all elements equal to one, \eqn{\mathbf{S}} is the vector of Mann-Kendall scores
#'  for each season and \eqn{\mathbf{\Gamma}}} denotes the variance - covariance matrix.
#'
#' @return
#' An object with class "htest"
#' \item{data.name}{character string that denotes the input data}
#' \item{p.value}{the p-value for the entire series}
#' \item{statistic}{the z quantile of the standard normal distribution
#' for the entire series}
#' \item{null.value}{the null hypothesis}
#' \item{estimates}{the estimates S and varS for the entire series}
#' \item{alternative}{the alternative hypothesis}
#' \item{method}{character string that denotes the test}
#' \item{cov}{the variance - covariance matrix}
#'
#' @note
#' Ties are not corrected. Current Version is for complete observations only.
#'
#' @inherit mk.test references 
#
#' @seealso
#' \code{\link{cor}},
#' \code{\link{cor.test}},
#' \code{\link{mk.test}},
#' \code{\link{smk.test}}
#' @examples
#' csmk.test(nottem)
#' @keywords ts nonparametric multivariate
#' @export
csmk.test <- function(x, alternative = c("two.sided", "greater", "less"))
{
    if(!is.ts(x)){
        stop("'x' must be objects of class 'ts'")
    }
    p <- frequency(x)
    if (p < 2){
        stop("'x' must have at least 2 seasons")
    }
 
    na.fail(x)
    alternative = match.arg(alternative)
    n <- length(x)

    dat <- matrix(NA, ncol = p, nrow = n / p)

    for (i in 1:p)
    {
        dat[,i] <- x[cycle(x) == i]
    }

    inp <- list(
        x = ts(dat),
        alternative = alternative)

    out <- do.call("mult.mk.test", inp)

    DNAME <- deparse(substitute(x))
    METHOD <- "Correlated Seasonal Mann-Kendall Test"
    out$method <- METHOD
    out$data.name <- DNAME

    return(out)
}
