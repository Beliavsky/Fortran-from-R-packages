\name{confint.onls}
\alias{confint.onls}

\title{Confidence intervals for 'onls' model parameters}

\description{
Computes confidence intervals for one or more parameters of an \code{\link{onls}} model. As in \code{MASS:::confint.nls}, these are based on profile likelihoods, using \code{onls:::profile.onls} and \code{onls:::confint.profile.onls}.
}

\usage{
\method{confint}{onls}(object, parm, level = 0.95, ...)
}

\arguments{
  \item{object}{an object returned from \code{\link{onls}}.}
  \item{parm}{a specification of which parameters are to be given confidence intervals, either a vector of numbers or a vector of names. If missing, all parameters are considered.}
  \item{level}{the confidence level required.}
  \item{...}{additional argument(s) for methods.}       
}

\details{
Profiling the likelihood uses the following strategy:\cr
If \eqn{\theta} is the parameter to be profiled and \eqn{\delta} the vector of remaining parameters,\cr
1) compute the log-likelihood of the model \eqn{\mathcal{L}(\theta^{*}, \delta^{*})} using the converged parameters,\cr
2) compute a lower bound \eqn{\theta^{*} - 0.6 \cdot \sigma(\theta^{*})} for the lower confidence limit,\cr
3) define a grid of values ranging from \eqn{\theta^{'}} to \eqn{\theta^{*}} (e.g., 100 equidistant points),\cr
4) for each grid value \eqn{\theta_i}, compute the profile log-likelihood value \eqn{\mathcal{L}_1(\theta_i)} by maximizing \eqn{\mathcal{L}(\theta_i, \delta)} over \eqn{\delta}-values by fixing \eqn{\theta} at \eqn{\theta_i},\cr 
5) find the confidence \code{level} by interpolation of the profile traces obtained from 4).
}

\value{
A matrix (or vector) with columns giving lower and upper confidence limits for each parameter. These will be labelled as (1 - \code{level})/2 and 1 - (1 - \code{level})/2 in \% (by default 2.5\% and 97.5\%).
}
     
\author{
Andrej-Nikolai Spiess, taken and modified from the \code{\link{nls}} functions.
}
        
\examples{
\donttest{
DNase1 <- subset(DNase, Run == 1)
DNase1$density <- sapply(DNase1$density, function(x) rnorm(1, x, 0.1 * x))
mod1 <- onls(density ~ Asym/(1 + exp((xmid - log(conc))/scal)), 
             data = DNase1, start = list(Asym = 3, xmid = 0, scal = 1))
confint(mod1)
}
}

\keyword{models}
\keyword{nonlinear}
