\name{vasicekfit}
\alias{vasicekfit}
\title{Fit an Extended Vasicek Credit Loss Model}
\description{
Fits the extended Vasicek single-factor credit loss model where the
probability of default depends on macroeconomic covariates. All parameters,
including asset value correlation, are estimated via closed-form
probit-transformed OLS.
}
\usage{
vasicekfit(formula, data, bias_correct = FALSE, portfolio_size = NULL)
}
\arguments{
  \item{formula}{a formula of the form \code{response ~ predictors} where the
    response is an observed default/loss rate in (0, 1) (or [0, 1] when
    \code{portfolio_size} is supplied). The model must include an intercept.}
  \item{data}{a data frame containing the variables in \code{formula}.}
  \item{bias_correct}{logical; if \code{TRUE}, apply the small-sample bias
    correction to the variance estimate by multiplying by \eqn{N/(N - m - 1)},
    where \eqn{N} is number of observations and \eqn{m} is number of predictors.}
  \item{portfolio_size}{optional positive integer. If supplied, a
    finite-portfolio variance correction is applied to the response before
    fitting (see Yang, 2014, section 4.3). The correction shrinks the
    observed rates towards their mean, mapping rates of exactly 0 or 1 into
    (0, 1), so such observations are permitted when \code{portfolio_size}
    is given.}
}
\value{
An object of class \code{"vasicekfit"}, which is a list containing:
\describe{
  \item{p}{estimated probability of default (baseline PD).}
  \item{rho}{estimated asset value correlation.}
  \item{kappa}{named numeric vector of macro-factor sensitivities.}
  \item{sigma2}{MLE variance estimate in probit space.}
  \item{lm_fit}{the underlying \code{\link{lm}} object.}
  \item{fitted.values}{fitted values in probit (Y) space.}
  \item{residuals}{residuals in probit (Y) space.}
  \item{formula}{the model formula.}
  \item{call}{the matched call.}
  \item{terms}{the \code{\link{terms}} object.}
  \item{model}{the model frame.}
  \item{bias_correct}{logical flag used.}
  \item{portfolio_size}{portfolio size used, or \code{NULL}.}
}
}
\references{
Vasicek, O. A. (2002). The distribution of loan portfolio value.
\emph{Risk}, 15(12), 160--162.

Yang, B. H. (2014). Estimating Long-Run PD, Asset Correlation, and
Portfolio Level PD by Vasicek Models. MPRA Paper No. 57244.
}
\examples{
set.seed(42)
n <- 100
u1 <- rnorm(n)
u2 <- rnorm(n)
p_true <- 0.03; rho_true <- 0.02
kappa_true <- c(0.13, -0.07)
z <- rnorm(n)
x <- pnorm((qnorm(p_true) + kappa_true[1] * u1 + kappa_true[2] * u2 +
  sqrt(rho_true) * z) / sqrt(1 - rho_true))
d <- data.frame(default_rate = x, unemp = u1, hpi = u2)
fit <- vasicekfit(default_rate ~ unemp + hpi, data = d)
fit
}
\seealso{
\code{\link{predict.vasicekfit}}, \code{\link{dvasicek}}
}
