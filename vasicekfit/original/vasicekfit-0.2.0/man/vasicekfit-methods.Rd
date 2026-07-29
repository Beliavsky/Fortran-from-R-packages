\name{vasicekfit-methods}
\alias{print.vasicekfit}
\alias{summary.vasicekfit}
\alias{coef.vasicekfit}
\alias{fitted.vasicekfit}
\alias{residuals.vasicekfit}
\alias{vcov.vasicekfit}
\alias{confint.vasicekfit}
\title{Methods for vasicekfit Objects}
\description{
Print, summary, coefficient extraction, covariance, confidence interval,
and accessor methods for \code{"vasicekfit"} objects.
}
\usage{
\method{print}{vasicekfit}(x, \dots)
\method{summary}{vasicekfit}(object, type = c("iid", "HAC"), \dots)
\method{coef}{vasicekfit}(object, \dots)
\method{fitted}{vasicekfit}(object, \dots)
\method{residuals}{vasicekfit}(object, \dots)
\method{vcov}{vasicekfit}(object, type = c("iid", "HAC"), \dots)
\method{confint}{vasicekfit}(object, parm, level = 0.95,
                              type = c("iid", "HAC"), \dots)
}
\arguments{
  \item{x, object}{a \code{\link{vasicekfit}} object.}
  \item{parm}{a specification of which parameters to compute confidence
    intervals for, either a vector of names or indices. If missing, all
    parameters are considered.}
  \item{level}{the confidence level required.}
  \item{type}{covariance type. \code{"iid"} (default) uses the closed-form
    block-diagonal Gaussian covariance. \code{"HAC"} uses a kernel long-run
    covariance of the joint influence functions for
    \eqn{(\hat\beta, \hat\sigma^2)} and requires the \pkg{sandwich} package.}
  \item{\dots}{when \code{type = "HAC"}, extra arguments are forwarded to
    \code{\link[sandwich]{lrvar}}, e.g. \code{prewhite}, \code{kernel},
    \code{bw}, or \code{adjust}. (\code{lrvar}'s own \code{type} argument
    cannot be set this way because the name collides; if you need
    Newey-West-style bandwidth selection, pass \code{kernel = "Bartlett"}
    and an explicit \code{bw}.) Ignored when \code{type = "iid"}.}
}
\details{
\code{coef} returns the recovered Vasicek parameters (p, rho, and kappas)
in the original parameter space.

\code{fitted} and \code{residuals} return values in probit (Y) space,
consistent with the underlying OLS regression.

\code{vcov} returns the delta-method covariance matrix for the
original-space parameters \eqn{(p, \rho, \kappa_1, \ldots, \kappa_m)},
obtained by propagating the joint covariance of \eqn{(\hat\beta, \hat\sigma^2)}
through the closed-form recovery transform.

With \code{type = "iid"} the \eqn{(\hat\beta, \hat\sigma^2)} block is
treated as block-diagonal by the independence of sample mean and variance
under normality. The variance of \eqn{\hat\sigma^2} uses
\eqn{2 \sigma^4 (N - m - 1)/N^2} for the MLE and \eqn{2 \sigma^4 / (N - m - 1)}
when \code{bias_correct = TRUE}.

With \code{type = "HAC"} the joint covariance of
\eqn{(\hat\beta, \hat\sigma^2)} is estimated as the long-run variance of
the stacked influence functions
\eqn{\psi_{\beta,i} = (X'X/N)^{-1} x_i e_i} and
\eqn{\psi_{\sigma^2,i} = e_i^2 - \hat\sigma^2}, computed via
\code{\link[sandwich]{lrvar}}. This is appropriate when the data exhibit
serial correlation (e.g. quarterly macro covariates), at the cost of some
finite-sample bias. Defaults follow \code{lrvar} (Andrews quadratic-spectral
kernel, prewhitened); pass \code{prewhite}, \code{kernel}, \code{bw}, or
\code{adjust} via \code{\dots} to override.

The HAC standard errors do not change with \code{bias_correct}: the
long-run variance is estimated empirically from the influence functions
rather than from the parametric Gaussian formula, so the
\eqn{N/(N-m-1)} scaling does not apply. The \code{bias_correct} flag
still affects the point estimates of \eqn{p}, \eqn{\rho}, and
\eqn{\kappa} (via \eqn{\hat\sigma^2}) regardless of \code{type}.

\code{confint} returns Wald intervals based on \code{vcov}.
}
\value{
\code{print} and \code{summary} return \code{object} invisibly.

\code{coef} returns a named numeric vector.

\code{fitted} and \code{residuals} return numeric vectors.

\code{vcov} returns a \eqn{(m+2) \times (m+2)} covariance matrix.

\code{confint} returns a two-column matrix of lower and upper bounds.
}
\seealso{
\code{\link{vasicekfit}}, \code{\link{predict.vasicekfit}}
}
