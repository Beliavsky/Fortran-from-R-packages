\name{CovOgk-class}
\Rdversion{1.1}
\docType{class}
\alias{CovOgk-class}

\title{ OGK Estimates of Multivariate Location and Scatter }
\description{
    This class, derived from the virtual class \code{"CovRobust"} accomodates
    OGK Estimates of multivariate location and scatter computed by the
    algorithm proposed by Marona and Zamar (2002).
}
\section{Objects from the Class}{
    Objects can be created by calls of the form \code{new("CovOgk", ...)},
    but the usual way of creating \code{CovOgk} objects is a call to the function
    \code{CovOgk} which serves as a constructor.
}
\section{Slots}{
    \describe{
        \item{\code{raw.cov}:}{Object of class \code{"matrix"} the raw
            (not reweighted) estimate of covariance matrix }
        \item{\code{raw.center}:}{Object of class \code{"vector"} - the raw
            (not reweighted) estimate of the location vector}
        \item{\code{raw.mah}:}{Object of class \code{"Uvector"} - mahalanobis
            distances of the observations based on the raw estimate of the
            location and scatter }
        \item{\code{raw.wt}:}{Object of class \code{"Uvector"} - weights of
            the observations based on the raw estimate of the location and scatter }
        \item{\code{iter}, \code{crit}, \code{wt}:}{from the
            \code{"\linkS4class{CovRobust}"} class.
        }
        \item{\code{call}, \code{cov}, \code{center},
            \code{n.obs}, \code{mah}, \code{method},
            \code{singularity}, \code{X}:}{from the \code{"\linkS4class{Cov}"} class.
        }
    }
}
\section{Extends}{
    Class \code{"\linkS4class{CovRobust}"}, directly.
    Class \code{"\linkS4class{Cov}"}, by class \code{"\linkS4class{CovRobust}"}.
}
\section{Methods}{
No methods defined with class "CovOgk" in the signature.
}
\references{
  Todorov V & Filzmoser P (2009), An Object Oriented Framework for Robust Multivariate Analysis. 
  \emph{Journal of Statistical Software}, \bold{32}(3), 1--47.
  \doi{10.18637/jss.v032.i03}.
}
\author{ Valentin Todorov \email{valentin.todorov@chello.at}
}
\seealso{ \code{\link{CovMcd-class}}, \code{\link{CovMest-class}} }
\examples{
showClass("CovOgk")
}
\keyword{classes}
\keyword{robust}
\keyword{multivariate}
