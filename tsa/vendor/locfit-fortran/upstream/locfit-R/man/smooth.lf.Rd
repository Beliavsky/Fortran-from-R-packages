\name{smooth.lf}
\alias{smooth.lf}
\title{
Local Regression, Likelihood and Density Estimation.
}
\usage{
smooth.lf(x, y, xev=x, direct=FALSE, ...)
}
\description{
  \code{smooth.lf} is a simple interface to the Locfit library.
  The input consists of a predictor vector (or matrix) and response.
  The output is a list with vectors of fitting points and fitted values.
  Most \code{\link{locfit.raw}} options are valid.
}
\arguments{
\item{x}{
Vector (or matrix) of the independent variable(s).
}
\item{y}{
Response variable. If omitted, \code{x} is treated as the response and
the predictor variable is \code{1:n}.
}
\item{xev}{
Fitting Points. Default is the data vector \code{x}.
}
\item{direct}{
Logical variable. If \code{T}, local regression is performed directly
at each fitting point. If \code{F}, the standard Locfit method combining
fitting and interpolation is used.
}
\item{...}{
  Other arguments to \code{\link{locfit.raw}()}.
}
}
\value{
A list with components \code{x} (fitting points) and \code{y} (fitted values).
Also has a \code{call} component, so \code{update()} will work.
}
\examples{
# using smooth.lf() to fit a local likelihood model.
data(morths)
fit <- smooth.lf(morths$age, morths$deaths, weights=morths$n,
                 family="binomial")
plot(fit,type="l")

# update with the direct fit
fit1 <- update(fit, direct=TRUE)
lines(fit1,col=2)
print(max(abs(fit$y-fit1$y)))
}
\seealso{
  \code{\link{locfit}()},
  \code{\link{locfit.raw}()},
  \code{\link{density.lf}()}.
}
\keyword{smooth}
