##   file smk.test.R part of package trend
##   Copyright (C) 2015-2018 Thorsten Pohlert
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
#' @title Seasonal Mann-Kendall Trend Test
#' @concept Hirsch-Slack-Test
#' @aliases Hirsch-Slack-Test
#' 
#' @description
#' Performs a Seasonal Mann-Kendall Trend Test (Hirsch-Slack Test)
#'
#' @param x a time series object with class \code{ts} comprising >= 2 seasons;
#' \code{NA} values are not allowed
#' @param alternative the alternative hypothesis, defaults to \code{two.sided}
#' @param continuity logical, indicates, whether a continuity correction
#' should be done; defaults to \code{TRUE}
#'
#' @details
#' The Mann-Kendall statistic for the $g$-th season is calculated as:
#'
#' \deqn{
#' S_g = \sum_{i = 1}^{n-1} \sum_{j = i + 1}^n
#' \mathrm{sgn}\left(x_{jg} - x_{ig}\right), \qquad (1 \le g \le m)}
#'
#' with \eqn{\mathrm{sgn}}{sgn} the signum function (see \code{\link{sign}}).
#' 
#' The mean of \eqn{S_g} is \eqn{\mu_g = 0}. The variance including the
#' correction term for ties is
#'
#' \deqn{
#' \sigma_g^2 = \left\{n \left(n-1\right)\left(2n+5\right) -
#' \sum_{j=1}^p t_{jg}\left(t_{jg} - 1\right)\left(2t_{jg}+5\right) \right\} / 18
#' ~~ (1 \le g \le m)}
#'
#' The seasonal Mann-Kendall statistic for the entire series is calculated
#' according to
#'
#' \deqn{
#' \begin{array}{ll}
#' \hat{S} = \sum_{g = 1}^m S_g &
#' \hat{\sigma}_g^2 =  \sum_{g = 1}^m \sigma_g^2
#' \end{array}}
#'
#' The statistic \eqn{S_g} is approximately normally distributed, with
#'
#' \deqn{z_g = S_g / \sigma_g} 
#'
#' If \code{continuity = TRUE} then a continuity correction will be employed:
#'
#' \deqn{z = \mathrm{sgn}(S_g) ~ \left(|S_g| - 1\right) / \sigma_g}
#'
#' @return
#' An object with class "htest" and "smktest"
#' \item{data.name}{character string that denotes the input data}
#' \item{p.value}{the p-value for the entire series}
#' \item{statistic}{the z quantile of the standard normal distribution
#' for the entire series}
#' \item{null.value}{the null hypothesis}
#' \item{estimates}{the estimates S and varS for the entire series}
#' \item{alternative}{the alternative hypothesis}
#' \item{method}{character string that denotes the test}
#' \item{Sg}{numeric vector that contains S scores for each season}
#' \item{varSg}{numeric vector that contains varS for each season}
#' \item{pvalg}{numeric vector that contains p-values for each season}
#' \item{Zg}{numeric vector that contains z-quantiles for each season}
#' 
#' @references
#' Hipel, K.W. and McLeod, A.I. (1994),
#' \emph{Time Series Modelling of Water Resources and Environmental Systems}. 
#' New York: Elsevier Science.
#'
#' Libiseller, C. and Grimvall, A. (2002), Performance of partial
#' Mann-Kendall tests for trend detection in the presence of covariates.
#' \emph{Environmetrics} 13, 71--84, \doi{10.1002/env.507}.
#'
#' R. Hirsch, J. Slack, R. Smith (1982), Techniques of Trend Analysis for
#' Monthly Water Quality Data, \emph{Water Resources Research} 18, 107--121.
#'
#' @examples
#' res <- smk.test(nottem)
#' ## print method
#' res
#' ## summary method
#' summary(res)
#'
#' @keywords ts
#' @keywords nonparametric
#' @keywords univar
#' @importFrom stats pnorm na.fail cycle
#' @export
smk.test <- function(x, alternative = c("two.sided", "greater", "less"),
                     continuity = TRUE)
{
    if(!is.ts(x)){
        stop("'x' must be an object of class 'ts'")
    }
    p <- frequency(x)
    if (p < 2){
        stop("'x' must have at least 2 seasons")
    }
    na.fail(x)
    alternative = match.arg(alternative)
    
    Sg <- numeric(p)
    varSg <- numeric(p)
    taug <- numeric(p)
    
    for (i in 1:p){
        dat <- x[cycle(x) == i]
        n <- length(dat)
        t <- table(dat)
        names(t) <- NULL
        Sg[i] <- .mkScore(dat)
        varSg[i] <- .varmk(t, n)
        taug[i] <- Sg[i] / .Dfn(t, n)
    }

    ## total statistics
    Sp <- sum(Sg)
    varSp <- sum(varSg)

    ## combine
    Stotal <- c(Sg, Sp)
    varStotal <- c(varSg, varSp)
    
    if(continuity){
        sg <- sign(Stotal)
        ztotal <- sg * (abs(Stotal) - 1) / sqrt(varStotal)
    } else {
        ztotal <- Stotal / sqrt(varStotal)
    }
    
    ## get the pvalue
    pval <- sapply(ztotal,
                   function(z) {
                       switch(alternative,
                              "two.sided" = {
                                  2 * min(0.5, pnorm(abs(z),
                                                     lower.tail=FALSE))
                              },

                              "greater" = {
                                  pnorm(z, lower.tail=FALSE)
                              },
                              
                              "less" = {
                                  pnorm(z, lower.tail=TRUE)
                              }
                              )
                   }
                   )
    
    ## finalise
    DNAME <- deparse(substitute(x))
    ans <- list(data.name = DNAME, 
                p.value = pval[p+1], 
                statistic = c(z = ztotal[p+1]),
                null.value = c(S = 0),
                estimates = c(S = Stotal[p+1], varS = varStotal[p+1]),
                alternative = alternative,
                method = "Seasonal Mann-Kendall trend test (Hirsch-Slack test)",
                Sg = Stotal[1:p],
                varSg = varStotal[1:p],
                pvalg = pval[1:p],
                Zg = ztotal[1:p],
                taug = taug)
    
    class(ans) <- c("htest", "smktest")
    return(ans)
}
