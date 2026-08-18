\name{rSquared}
\alias{rSquared}
\title{Calculate R squared value}

\description{
   Calculate R squared value.
}

\usage{rSquared( y, resid )}

\arguments{
   \item{ y }{vector of endogenous variables}
   \item{ resid }{vector of residuals}
}

\author{Arne Henningsen}

\examples{
   data( "Electricity", package = "Ecdat" )
   reg <- lm( cost ~ q + pl + pk + pf, Electricity )
   rSquared( Electricity$cost, reg$residuals )
   summary( reg )$r.squared  # returns the same value
}

\keyword{univar}
\keyword{multivariate}
\keyword{array}

