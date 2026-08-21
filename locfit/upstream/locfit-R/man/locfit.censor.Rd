\name{locfit.censor}
\alias{locfit.censor}

\title{
Censored Local Regression
}
\usage{
locfit.censor(x, y, cens, ..., iter=3, km=FALSE)
}
\description{
  \code{locfit.censor} produces local regression estimates for censored
  data. The basic idea is to use an EM style algorithm, where one
  alternates between estimating the regression and the true values
  of censored observations.

  \code{locfit.censor} is designed as a front end
  to \code{\link{locfit.raw}} with data vectors, or as an intemediary
  between \code{\link{locfit}} and \code{\link{locfit.raw}} with a
  model formula. If you can stand the syntax, the second calling
  sequence above will be slightly more efficient than the third.
}
\arguments{
\item{x}{
  Either a \code{\link{locfit}} model formula or a numeric vector
  of the predictor variable.
}
\item{y}{
  If \code{x} is numeric, \code{y} gives the response variable.
}
\item{cens}{
  Logical variable indicating censoring. The coding is \code{1}
  or \code{TRUE} for censored; \code{0} or \code{FALSE} for uncensored.
}
\item{...}{
  Other arguments to \code{\link{locfit.raw}}
}
\item{iter}{Number of EM iterations to perform}
\item{km}{
  If \code{km=TRUE}, the estimation of censored observations uses
  the Kaplan-Meier estimate, leading to a local version of the
  Buckley-James estimate. If \code{km=F}, the estimation is based
  on a normal model (Schmee and Hahn). Beware of claims that B-J
  is nonparametric; it makes stronger assumptions on the upper tail
  of survival distributions than most authors care to admit.
}
}
\value{
\code{locfit} object.
}
\seealso{
\code{\link{km.mrl}},
\code{\link{locfit}},
\code{\link{locfit.raw}}
}
\examples{
data(heart, package="locfit")
fit <- locfit.censor(log10(surv+0.5) ~ age, cens=cens, data=heart)
## Can also be written as:
\dontrun{fit <- locfit(log10(surv + 0.5) ~ age, cens=cens, data=heart, lfproc=locfit.censor)}
with(heart, plotbyfactor(age, 0.5 + surv, cens, ylim=c(0.5, 16000), log="y"))
lines(fit, tr=function(x) 10^x)
}
\references{
Buckley, J. and James, I. (1979). Linear Regression with censored data.
  Biometrika 66, 429-436.

Loader, C. (1999). Local Regression and Likelihood. Springer, NY (Section 7.2).

Schmee, J. and Hahn, G. J. (1979). A simple method for linear regression
  analysis with censored data (with discussion). Technometrics 21, 417-434.
}
\keyword{smooth}
% Converted by Sd2Rd version 0.2-a5.
