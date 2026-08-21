\name{km.mrl}
\alias{km.mrl}
\title{
Mean Residual Life using Kaplan-Meier estimate
}
\usage{
km.mrl(times, cens)
}
\description{
  This function computes the mean residual life for censored data
  using the Kaplan-Meier estimate of the survival function. If
  \eqn{S(t)} is the K-M estimate, the MRL for a censored observation
  is computed as \eqn{(\int_t^{\infty} S(u)du)/S(t)}. We take
  \eqn{S(t)=0} when \eqn{t} is greater than the largest observation,
  regardless of whether that observation was censored.

  When there are ties between censored and uncensored observations,
  for definiteness our ordering places the censored observations
  before uncensored.

  This function is used by \code{\link{locfit.censor}} to compute
  censored regression estimates.
}
\arguments{
\item{times}{
  Obsereved survival times.
}
\item{cens}{
  Logical variable indicating censoring. The coding is \code{1}
  or \code{TRUE} for censored; \code{0} or \code{FALSE} for uncensored.
}
}
\value{
A vector of the estimated mean residual life. For uncensored observations,
the corresponding estimate is 0.
}
\examples{
# censored regression using the Kaplan-Meier estimate.
data(heart, package="locfit")
fit <- locfit.censor(log10(surv+0.5)~age, cens=cens, data=heart, km=TRUE)
plotbyfactor(heart$age, 0.5+heart$surv, heart$cens, ylim=c(0.5,16000), log="y")
lines(fit, tr=function(x)10^x)
}
\references{
Buckley, J. and James, I. (1979). Linear Regression with censored data.
  Biometrika 66, 429-436.

Loader, C. (1999). Local Regression and Likelihood. Springer, NY (Section 7.2).
}
\seealso{
\code{\link{locfit.censor}}
}
%\keyword{locfit}
\keyword{smooth}
% Converted by Sd2Rd version 0.2-a5.
