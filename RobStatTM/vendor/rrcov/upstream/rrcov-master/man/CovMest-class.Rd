\name{CovMest-class}
\Rdversion{1.1}
\docType{class}
\alias{CovMest-class}

\title{Constrained M-estimates of Multivariate Location and Scatter }
\description{
    This class, derived from the virtual class "CovRobust" accomodates
    constrained M-Estimates of multivariate location and scatter
    based on the translated biweight function (\sQuote{t-biweight}) using
    a High breakdown point initial estimate (Minimum Covariance Determinant -
    \sQuote{Fast MCD})}
\section{Objects from the Class}{
    Objects can be created by calls of the form \code{new("CovMest", ...)},
    but the usual way of creating \code{CovMest} objects is a call to the function
    \code{CovMest} which serves as a constructor.

}
\section{Slots}{
    \describe{
        \item{\code{vt}:}{Object of class \code{"vector"} - vector of weights (v)}
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
   No methods defined with class "CovMest" in the signature.
}
\references{
  Todorov V & Filzmoser P (2009), An Object Oriented Framework for Robust Multivariate Analysis. 
  \emph{Journal of Statistical Software}, \bold{32}(3), 1--47.
  \doi{10.18637/jss.v032.i03}.
}
\author{Valentin Todorov \email{valentin.todorov@chello.at}}
%\note{}
\seealso{
   \code{\link{CovMest}}, \code{\link{Cov-class}}, \code{\link{CovRobust-class}}
}
\examples{
showClass("CovMest")
}
\keyword{classes}
\keyword{robust}
\keyword{multivariate}
