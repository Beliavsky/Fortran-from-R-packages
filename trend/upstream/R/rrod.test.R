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
##     2017-11-13
##     - initial writing
##

#' @name rrod.test
#' @title Robust Rank-Order Distributional Test
#' @description
#' Performs Fligner-Pollicello robust rank-order distributional test
#' for location.
#'
#' @details
#' The non-parametric RROD two-sample test can be used to test
#' for differences in location, whereas it does not assume variance
#' homogeneity.
#'
#' Let \eqn{X} and \eqn{Y} denote two samples with sizes \eqn{n_x}{nx} and \eqn{n_y}{ny}
#' of a continuous variable.First, the combined sample is transformed
#' into ranks in increasing order.
#' Let \eqn{S_{xi}}{Sx} and \eqn{S_{yj}}{Sy} denote the counts of \eqn{Y} \eqn{(X)}
#' values having a lower rank than \eqn{x_i} \eqn{(y_j)}. The mean counts are:
#' 
#' \deqn{\bar{S}_x = \sum_{i=1}^{n_x} S_{xi} / n_x}{%%
#' Sx = sum(Sxi) / nx}
#'
#' \deqn{\bar{S}_y = \sum_{j=1}^{n_y} S_{yj} / n_y}{%%
#' Sy = sum(Syi) / ny}
#'
#' The variances are:
#' \deqn{s^2_{Sx} = \sum_{i=1}^{n_x} \left( S_{xi} - \bar{S}_x \right)^2}{%%
#' sSxsq = sum((Sxi - Sx)^2)}
#'
#' \deqn{s^2_{Sy} = \sum_{j=1}^{n_y} \left( S_{yj} - \bar{S}_y \right)^2}{%%
#' sSysq = sum((Syj - Sy)^2)}
#'
#' The test statistic is:
#' \deqn{
#' z = \frac{1}{2}~
#' \frac{n_x \bar{S}_x - n_y \bar{S}_y}
#' {\left( \bar{S}_x \bar{S}_y + s^2_{Sx} + s^2_{Sy} \right)^{1/2}}}{%%
#' z = 1/2 * (nx * Sx - ny * Sy) / (Sx * Sy + sSxsq + sSysq)^0.5}
#'
#' The two samples have significantly different location parameters,
#' if \eqn{|z| > z_{1-\alpha/2}}{|z| > z(1-alpha/2)}.
#' The function calculates the \eqn{p}-values of the null hypothesis
#' for the selected alternative than can be \code{"two.sided"}, \code{"greater"}
#' or \code{"less"}.
#' 
#' @return
#' A list with class \code{"htest"}.
#' 
#' @seealso
#' \code{\link[stats]{wilcox.test}}
#' 
#' @keywords htest nonparametric
#'
#' @references
#' Fligner, M. A., Pollicello, G. E. III. (1981), Robust Rank Procedures for 
#' the Behrens-Fisher Problem, \emph{Journal of the 
#' American Statistical Association}, 76, 162--168.
#' 
#' Lanzante, J. R. (1996), Resistant, robust and non-parametric
#' techniques for the analysis of climate data: Theory and examples,
#' including applications to historical radiosonde station data,
#' \emph{Int. J. Clim.}, 16, 1197--1226.
#' 
#' Siegel, S. and Castellan, N. (1988),
#' \emph{Nonparametric Statistics For The Behavioural Sciences},
#' New York: McCraw-Hill.
#' 
#' @examples
#' ## Two-sample test.
#' ## Hollander & Wolfe (1973), 69f.
#' ## Permeability constants of the human chorioamnion (a placental
#' ##  membrane) at term (x) and between 12 to 26 weeks gestational
#' ##  age (y).  The alternative of interest is greater permeability
#' ##  of the human chorioamnion for the term pregnancy.
#' x <- c(0.80, 0.83, 1.89, 1.04, 1.45, 1.38, 1.91, 1.64, 0.73, 1.46)
#' y <- c(1.15, 0.88, 0.90, 0.74, 1.21)
#' rrod.test(x, y, alternative = "g")
#'
#' ## Formula interface.
#' boxplot(Ozone ~ Month, data = airquality)
#' rrod.test(Ozone ~ Month, data = airquality,
#'             subset = Month %in% c(5, 8)) 
#' @export
rrod.test <- function(x, ...) UseMethod("rrod.test")

#' @rdname rrod.test
#' @aliases rrod.test.default
#' @method rrod.test default
#' @param x a vector of data values.
#' @param y an optional numeric vector of data values.
#' @param alternative the alternative hypothesis. Defaults to \code{"two.sided"}.
#' @importFrom stats pnorm
#' @export
rrod.test.default <-
function(x, y, alternative = c("two.sided", "less", "greater"), ...)
{
    alternative <- match.arg(alternative)

    if(!is.numeric(x)) stop("'x' must be numeric")
   
    if(!is.numeric(y)) stop("'y' must be numeric")
    DNAME <- paste(deparse(substitute(x)), "and",
                   deparse(substitute(y)))
  
    x <- x[is.finite(x)]
    y <- y[is.finite(y)]
        
    ## sample sizes
    n.x <- length(x)
    n.y <- length(y)
    n <- n.x + n.y
    R <- rank(c(x, y))
    rx <- R[1:n.x]
    ry <- R[(n.x + 1):n]
    NXi <- sapply(rx, function(r) sum(ifelse(r > ry, 1, 0)))
    NYi  <- sapply(ry, function(r) sum(ifelse(r > rx, 1, 0)))

    NX <- mean(NXi)
    NY <- mean(NYi)
    ssqnx <- sum((NXi - NX)^2)
    ssqny <- sum((NYi - NY)^2)

    Z <- 0.5 * ( n.x * NX - n.y * NY) / sqrt(NX * NY + ssqnx + ssqny)
    PVAL <- switch(alternative,
                   two.sided = 2 * min(0.5, pnorm(abs(Z), lower.tail=FALSE)),
                   greater = pnorm(Z, lower.tail=FALSE),
                   less = pnorm(Z)
                   )
    STATISTIC <- c(z = Z)
    METHOD <- "Robust Rank-Order Distributional Test"
    RVAL <- list(statistic = STATISTIC,
                 parameter = NULL,
                 p.value = as.numeric(PVAL),
                 alternative = alternative,
                 method = METHOD,
                 null.value = c("location shift" = 0),
                 data.name = DNAME)
    class(RVAL) <- "htest"
    RVAL
}

#' @rdname rrod.test
#' @aliases rrod.test.formula
#' @method rrod.test formula
#' @param formula a formula of the form \code{response ~ group} where
#'    \code{response} gives the data values and \code{group} a vector or
#'    factor of the corresponding groups.
#' @param data an optional matrix or data frame (or similar: see
#'  \code{\link{model.frame}}) containing the variables in the
#'  formula \code{formula}.  By default the variables are taken from
#'  \code{environment(formula)}.
#' @param subset an optional vector specifying a 
#'  subset of observations to be used.
#' @param na.action a function which indicates what should happen when
#'    the data contain \code{NA}s.  Defaults to \code{getOption("na.action")}.
#' @param \dots further arguments to be passed to or from methods.
#' @importFrom stats terms
#' @importFrom stats setNames
#' @export
rrod.test.formula <-
function(formula, data, subset, na.action, ...)
{
    if(missing(formula)
       || (length(formula) != 3L)
       || (length(attr(terms(formula[-2L]), "term.labels")) != 1L))
        stop("'formula' missing or incorrect")
    m <- match.call(expand.dots = FALSE)
    if(is.matrix(eval(m$data, parent.frame())))
        m$data <- as.data.frame(data)
    ## need stats:: for non-standard evaluation
    m[[1L]] <- quote(stats::model.frame)
    m$... <- NULL
    mf <- eval(m, parent.frame())
    DNAME <- paste(names(mf), collapse = " by ")
    names(mf) <- NULL
    response <- attr(attr(mf, "terms"), "response")
    g <- factor(mf[[-response]])
    if(nlevels(g) != 2L)
        stop("grouping factor must have exactly 2 levels")
    DATA <- setNames(split(mf[[response]], g),
		     c("x", "y"))
    y <- do.call("rrod.test", c(DATA, list(...)))
    y$data.name <- DNAME
    y
}
