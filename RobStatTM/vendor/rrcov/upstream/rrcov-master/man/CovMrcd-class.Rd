\name{CovMrcd-class}
\Rdversion{1.1}
\docType{class}
\alias{CovMrcd-class}

\title{ MRCD Estimates of Multivariate Location and Scatter }
\description{
    This class, derived from the virtual class \code{"CovRobust"} accomodates
    MRCD Estimates of multivariate location and scatter computed by a variant of the 
    \sQuote{Fast MCD} algorithm.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("CovMrcd", ...)},
but the usual way of creating \code{CovMrcd} objects is a call to the function
\code{CovMrcd} which serves as a constructor.
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
        \item{\code{cnp2}:}{Object of class \code{"numeric"} - containing the consistency 
            correction factor of the estimate of the covariance matrix.
        }
        \item{\code{icov}:}{The inverse of the covariance matrix.
        }
        \item{\code{rho}:}{The estimated regularization parameter.
        }
        \item{\code{target}:}{The estimated target matrix.
        }
        \item{\code{crit}:}{from the
            \code{"\linkS4class{CovRobust}"} class.
        }
        \item{\code{call}, \code{cov}, \code{center},
            \code{n.obs}, \code{mah}, \code{method},
            \code{X}:}{from the \code{"\linkS4class{Cov}"} class.
        }
    }
}
\section{Extends}{
    Class \code{"\linkS4class{CovRobust}"}, directly.
    Class \code{"\linkS4class{Cov}"}, by class \code{"\linkS4class{CovRobust}"}.
}
\section{Methods}{
    No methods defined with class \code{"CovMrcd"} in the signature.
}
\references{
  Todorov V & Filzmoser P (2009), An Object Oriented Framework for Robust Multivariate Analysis. 
  \emph{Journal of Statistical Software}, \bold{32}(3), 1--47.
  \doi{10.18637/jss.v032.i03}.
}
\author{ Valentin Todorov \email{valentin.todorov@chello.at} }
\seealso{
   \code{\link{CovMrcd}}, \code{\link{Cov-class}}, \code{\link{CovRobust-class}}, \code{\link{CovMcd-class}}
}
\examples{
showClass("CovMrcd")
}
\keyword{classes}
\keyword{robust}
\keyword{multivariate}
