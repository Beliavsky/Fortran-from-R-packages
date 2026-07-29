\name{predict.vasicekfit}
\alias{predict.vasicekfit}
\title{Predict Method for vasicekfit Objects}
\description{
Obtain predictions from a fitted Vasicek model, optionally at specified
confidence levels.
}
\usage{
\method{predict}{vasicekfit}(object, newdata = NULL,
        type = c("link", "response"), alpha = NULL, \dots)
}
\arguments{
  \item{object}{a \code{\link{vasicekfit}} object.}
  \item{newdata}{an optional data frame of new predictor values. If omitted,
    the training data are used.}
  \item{type}{character; \code{"link"} (default) returns predictions in probit
    space, \code{"response"} returns the conditional mean loss rate (the
    effective PD) on the (0, 1) scale.}
  \item{alpha}{optional numeric vector of confidence levels in (0, 1). When
    supplied, predictions are conditional quantiles of the loss distribution
    at each confidence level.}
  \item{\dots}{additional arguments (currently unused).}
}
\details{
With \code{type = "response"} and \code{alpha = NULL}, the prediction is
the conditional mean of the loss distribution given the covariates (the
effective PD), \eqn{\Phi(\Phi^{-1}(p) + \sum_j \kappa_j u_j)}, matching
the convention of \code{\link{predict.glm}}. Note that this is not
\code{pnorm} of the \code{"link"} prediction: the model has an error term
in probit space and the probit transform is nonlinear, so the
back-transformed linear predictor \eqn{\Phi(\hat\eta)} is the conditional
\emph{median} of the loss distribution (available as \code{alpha = 0.5}),
which for \eqn{p < 0.5} is smaller than the mean.
}
\value{
If \code{alpha} is \code{NULL}, a numeric vector of predictions.
If \code{alpha} is supplied, a matrix with rows corresponding to observations
and columns to confidence levels. When \code{alpha} is a scalar, the result
is a vector.
}
\examples{
set.seed(42)
n <- 100
u <- rnorm(n)
x <- pnorm((qnorm(0.03) + 0.1 * u + sqrt(0.02) * rnorm(n)) / sqrt(1 - 0.02))
d <- data.frame(default_rate = x, macro = u)
fit <- vasicekfit(default_rate ~ macro, data = d)

# Conditional mean PD (expected loss rate given the covariates)
predict(fit, type = "response")

# Conditional median loss rate
predict(fit, alpha = 0.5, type = "response")

# 99th percentile loss rate under stress
predict(fit, newdata = data.frame(macro = 2), alpha = 0.99,
        type = "response")

# Multiple confidence levels
predict(fit, newdata = data.frame(macro = c(0, 1, 2)),
        alpha = c(0.95, 0.99, 0.999), type = "response")
}
\seealso{
\code{\link{vasicekfit}}
}
