\name{ConfRadius}
\alias{ConfRadius}
\title{Helper to compute confidence ellipsoids}
\description{
  Computes the quantile of the Mahalanobis distance needed to draw
  confidence ellipsoids.
}
\usage{
ConfRadius(model,prob=1-alpha,alpha)
          }
\arguments{
  \item{model}{A multivariate linear model}
  \item{prob}{The confidence probability}
  \item{alpha}{The alpha error allowed, i.e. the complement of the
    confidence probability}
}
\value{
  a scalar
}
\details{
  Calculates the radius to be used in confidence ellipses for the
  parameters based on the Hottelings \eqn{T^2}{T^2} distribution. 
}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{lm}}, \code{\link{mvar}}, \code{\link{AIC}} 
}
\examples{
data(SimulatedAmounts)
model <- lm(ilr(sa.groups)~sa.groups.area)
cf = coef(model)
plot(ilrInv(cf, x=sa.groups))
for(i in 1:nrow(cf)){
  vr = vcovAcomp(model)[,,i,i]
  vr = ilrvar2clr(vr)
  ellipses(ilrInv(cf[i,]), vr, r=ConfRadius(model, alpha=0.05) )
 }
}
\keyword{multivariate}
