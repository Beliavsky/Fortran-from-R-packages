\name{CovClassic}
\alias{CovClassic}
\alias{Cov}
\title{Classical Estimates of Multivariate Location and Scatter}
\concept{High breakdown point}
\description{
    Computes the classical estimates of multivariate location and scatter.
    Returns an S4 class \code{CovClassic} with the estimated \code{center},
    \code{cov}, Mahalanobis distances and weights based on these distances.
}
\usage{
    CovClassic(x, unbiased=TRUE)
    Cov(x, unbiased=TRUE)
}
\arguments{
  \item{x}{a matrix or data frame. As usual, rows are observations and
        columns are variables. }
  \item{unbiased}{ whether to return the unbiased estimate of
    the covariance matrix. Default is  \code{unbiased = TRUE}}
}
%\details{}
\value{
  An object of class \code{"CovClassic"}.
}
\references{
  Todorov V & Filzmoser P (2009), An Object Oriented Framework for Robust Multivariate Analysis. 
  \emph{Journal of Statistical Software}, \bold{32}(3), 1--47.
  \doi{10.18637/jss.v032.i03}.
}
%\note{}
\author{Valentin Todorov \email{valentin.todorov@chello.at}}
\seealso{
  \code{\link{Cov-class}}, \code{\link{CovClassic-class}}
}
\examples{
data(hbk)
hbk.x <- data.matrix(hbk[, 1:3])
cv <- CovClassic(hbk.x)
cv
summary(cv)
plot(cv)
}
\keyword{classes}
\keyword{robust}
\keyword{multivariate}
