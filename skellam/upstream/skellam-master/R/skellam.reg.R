#' Skellam Regression
#'
#' Fits a regression model assuming a Skellam distribution for the response variable.
#'
#' @param y A vector of integers (positive or negative)
#' @param x A matrix, vector or data.frame of covariates
#'
#' @details
#' The function uses an exponential link function to ensure positive values for both
#' rate parameters (\eqn{\lambda_1} and \eqn{\lambda_2}). Optimization is performed
#' using \code{\link[stats]{nlm}}.
#'
#' @return
#' A list with components:
#' \describe{
#'   \item{loglik}{Maximized log-likelihood value}
#'   \item{param1}{Matrix for \eqn{\lambda_1} parameters:
#'     \itemize{
#'       \item Column 1: Estimated coefficients
#'       \item Column 2: Standard errors
#'       \item Column 3: t-values (coef/se)
#'       \item Column 4: p-values (Wald test)
#'     }
#'   }
#'   \item{param2}{Matrix for \eqn{\lambda_2} parameters (same structure as param1)}
#' }
#'
#' @references
#' \itemize{
#'   \item Skellam, J. G. (1946) The frequency distribution of the difference between
#'   two Poisson variates belonging to different populations. \emph{Journal of the
#'   Royal Statistical Society, Series A} \bold{109}(3), 296.
#'
#'   \item Strackee, J.; van der Gon, J. J. D. (1962) The frequency distribution of
#'   the difference between two Poisson variates. \emph{Statistica Neerlandica}
#'   \bold{16}(1), 17-23.
#'
#'   \item Karlis D. and Ntzoufras I. (2009) \emph{Analysis of sports data using bivariate
#'   Poisson models}. IMA Conference Presentation.
#'   \url{http://www2.stat-athens.aueb.gr/~jbn/papers/files/20_Karlis_Ntzoufras_2009_IMA_presentation_handouts_v01.pdf}
#' }
#'
#' @examples
#' set.seed(0)
#' x <- rnorm(100)
#' y1 <- rpois(100, exp(1 + 1 * x))
#' y2 <- rpois(100, exp(-1 + 1 * x))
#' y <- y2 - y1
#' skellam.reg(y, x)
#'
#' @export
#' @keywords models regression
#' @author Michail Tsagris
#' @export
skellam.reg <- function(y, x) {

  n <- length(y)
  x <- stats::model.matrix( ~., data.frame(x) )
  p <- dim(x)[2]

  skelreg <- function(pa) {
    b1 <- pa[1:p]   ;   b2 <- pa[ -c(1:p) ]
    a1 <- x %*% b1      ;     a2 <- x %*% b2
    lam1 <- exp(a1)     ;     lam2 <- exp(a2)
    a <- 2 * sqrt(lam1 * lam2)
    sum(lam1 + lam2)  + 0.5 * sum(y * (a1 - a2) ) - sum( log( besselI(a, y) ) )
  }

  options(warn = -1)
  mod <- stats::nlm(skelreg, stats::rnorm(2 * p), iterlim = 5000 )
  mod <- stats::nlm(skelreg, mod$estimate, iterlim = 5000 )
  mod <- stats::optim(mod$estimate, skelreg, hessian = TRUE, control = list(maxit = 5000) )
  b1 <- mod$par[1:p]    ;    b2 <- mod$par[ -c(1:p) ]
  s <- diag( solve(mod$hessian) )
  s1 <- sqrt(s[1:p])    ;    s2 <- sqrt(s[ -c(1:p) ])
  param1 <- cbind(b1, s1, b1 / s1, stats::pchisq( (b1 / s1)^2, 1, lower.tail = FALSE) )
  param2 <- cbind(b2, s2, b2 / s2, stats::pchisq( (b2 / s2)^2, 1, lower.tail = FALSE) )
  rownames(param1) <- rownames(param2) <- colnames(x)
  colnames(param1) <- colnames(param2) <- c("Estimate", "Std. Error", "Wald value", "p-value")

  list(loglik = -mod$value, param1 = param1, param2 = param2)
}

