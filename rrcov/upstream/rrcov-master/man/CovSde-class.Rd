\name{CovSde-class}
\Rdversion{1.1}
\docType{class}
\alias{CovSde-class}

\title{
Stahel-Donoho Estimates of Multivariate Location and Scatter
}
\description{ This class, derived from the virtual class \code{"CovRobust"}
    accomodates Stahel-Donoho estimates of multivariate location and scatter.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("CovSde", ...)},
but the usual way of creating \code{CovSde} objects is a call to the function
\code{CovSde} which serves as a constructor.
}
\section{Slots}{
    \describe{
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
No methods defined with class "CovSde" in the signature.
}
\references{
  Todorov V & Filzmoser P (2009), An Object Oriented Framework for Robust Multivariate Analysis. 
  \emph{Journal of Statistical Software}, \bold{32}(3), 1--47.
  \doi{10.18637/jss.v032.i03}.
}
\author{ Valentin Todorov \email{valentin.todorov@chello.at} }
\seealso{
   \code{\link{CovSde}}, \code{\link{Cov-class}}, \code{\link{CovRobust-class}}
}
\examples{
showClass("CovSde")
}
\keyword{classes}
\keyword{robust}
\keyword{multivariate}
