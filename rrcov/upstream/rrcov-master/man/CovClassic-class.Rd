\name{CovClassic-class}
\docType{class}
\alias{CovClassic-class}
\title{Class "CovClassic" - classical estimates of multivariate
    location and scatter }
\description{ The class \code{CovClassic} represents an estimate of the
    multivariate location and scatter of a data set. The objects of class \code{CovClassic}
    contain the classical estimates.
}
\section{Objects from the Class}{
    Objects can be created by calls of the form \code{new("CovClassic", ...)},
    but the usual way of creating \code{CovClassic} objects is a call to the function
    \code{CovClassic} which serves as a constructor.
}
\section{Slots}{
    \describe{
    \item{\code{call}:}{Object of class \code{"language"} }
    \item{\code{cov}:}{covariance matrix}
    \item{\code{center}:}{location}
    \item{\code{n.obs}:}{number of observations used for the computation of the estimates}
    \item{\code{mah}:}{mahalanobis distances}
    \item{\code{method}:}{a character string describing the method used to
        compute the estimate: "Classic"}
    \item{\code{singularity}:}{a list with singularity information for the
        ocvariance matrix (or \code{NULL} of not singular)}
    \item{\code{X}:}{data}
  }
}
\section{Methods}{
  \describe{
    \item{getCenter}{\code{signature(obj = "CovClassic")}: location vector }
    \item{getCov}{\code{signature(obj = "CovClassic")}: covariance matrix }
    \item{getCorr}{\code{signature(obj = "CovClassic")}: correlation matrix }
    \item{getData}{\code{signature(obj = "CovClassic")}: data frame }
    \item{getDistance}{\code{signature(obj = "CovClassic")}: distances }
    \item{getEvals}{\code{signature(obj = "CovClassic")}: Computes and returns
        the eigenvalues of the covariance matrix }
    \item{plot}{\code{signature(x = "CovClassic")}: plot the object }
    \item{show}{\code{signature(object = "CovClassic")}: display the object }
    \item{summary}{\code{signature(object = "CovClassic")}: calculate summary information }
     }
}
\references{
  Todorov V & Filzmoser P (2009), An Object Oriented Framework for Robust Multivariate Analysis. 
  \emph{Journal of Statistical Software}, \bold{32}(3), 1--47.
  \doi{10.18637/jss.v032.i03}.
}
%\note{}
\author{Valentin Todorov \email{valentin.todorov@chello.at}}
%\seealso{}
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
