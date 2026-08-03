\name{CovMcd-class}
\Rdversion{1.1}
\docType{class}
\alias{CovMcd-class}

\title{ MCD Estimates of Multivariate Location and Scatter }
\description{
    This class, derived from the virtual class \code{"CovRobust"} accomodates
    MCD Estimates of multivariate location and scatter computed by the
    \sQuote{Fast MCD} algorithm.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("CovMcd", ...)},
but the usual way of creating \code{CovMcd} objects is a call to the function
\code{CovMcd} which serves as a constructor.
}
\section{Slots}{
    \describe{
        \item{\code{alpha}:}{Object of class \code{"numeric"} - the size of the
            subsets over which the determinant is minimized (the default is (n+p+1)/2) }
        \item{\code{quan}:}{Object of class \code{"numeric"} - the number of
            observations on which the MCD is based. If \code{quan} equals
            \code{n.obs}, the MCD is the classical covariance matrix. }
        \item{\code{best}:}{Object of class \code{"Uvector"} - the best subset
            found and used for computing the raw estimates. The size of \code{best}
            is equal to \code{quan} }
        \item{\code{raw.cov}:}{Object of class \code{"matrix"} the raw
            (not reweighted) estimate of location }
        \item{\code{raw.center}:}{Object of class \code{"vector"} - the raw
            (not reweighted) estimate of scatter}
        \item{\code{raw.mah}:}{Object of class \code{"Uvector"} - mahalanobis
            distances of the observations based on the raw estimate of the
            location and scatter }
        \item{\code{raw.wt}:}{Object of class \code{"Uvector"} - weights of
            the observations based on the raw estimate of the location and scatter }
        \item{\code{raw.cnp2}:}{Object of class \code{"numeric"} - a vector of length
            two containing the consistency correction factor and the finite sample
            correction factor of the raw estimate of the covariance matrix }
        \item{\code{cnp2}:}{Object of class \code{"numeric"} - a vector of length two
            containing the consistency correction factor and the finite sample
            correction factor of the final estimate of the covariance matrix.
        }
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
No methods defined with class \code{"CovMcd"} in the signature.
}
\references{
  Todorov V & Filzmoser P (2009), An Object Oriented Framework for Robust Multivariate Analysis. 
  \emph{Journal of Statistical Software}, \bold{32}(3), 1--47.
  \doi{10.18637/jss.v032.i03}.
}
\author{ Valentin Todorov \email{valentin.todorov@chello.at} }
\seealso{
   \code{\link{CovMcd}}, \code{\link{Cov-class}}, \code{\link{CovRobust-class}}
}
\examples{
showClass("CovMcd")
}
\keyword{classes}
\keyword{robust}
\keyword{multivariate}
