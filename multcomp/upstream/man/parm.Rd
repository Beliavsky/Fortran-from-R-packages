\name{parm}
\alias{parm}
\title{ Model Parameters }
\description{
  Directly specify estimated model parameters and their covariance matrix.
}
\usage{
parm(coef, vcov, df = 0)
}
\arguments{
  \item{coef}{ estimated coefficients.}
  \item{vcov}{ estimated covariance matrix of the coefficients.}
  \item{df}{ an optional specification of the degrees of freedom to be used in
             subsequent computations. }
}
\details{
  
  When only estimated model parameters and the corresponding
  covariance matrix is available for simultaneous inference
  using \code{\link{glht}} (for example, when only the results
  but not the original data are available or, even worse, when the model
  has been fitted outside R), function \code{parm} sets up an
  object \code{\link{glht}} is able to compute on (mainly 
  by offering \code{coef} and \code{vcov} methods).

  Note that the linear function in \code{\link{glht}} can't 
  be specified via \code{\link{mcp}} since the model terms
  are missing.

}
\value{
  An object of class \code{parm} with elements
  \item{coef}{model parameters}
  \item{vcov}{covariance matrix of model parameters}
  \item{df}{degrees of freedom}
}
\examples{

## example from
## Bretz, Hothorn, and Westfall (2002). 
## On multiple comparisons in R. R News, 2(3):14-17.

beta <- c(V1 = 14.8, V2 = 12.6667, V3 = 7.3333, V4 = 13.1333)
Sigma <- 6.7099 * (diag(1 / c(20, 3, 3, 15)))
confint(glht(model = parm(beta, Sigma, 37),
             linfct = c("V2 - V1 >= 0", 
                        "V3 - V1 >= 0", 
                        "V4 - V1 >= 0")), 
        level = 0.9)

}
\keyword{misc}
