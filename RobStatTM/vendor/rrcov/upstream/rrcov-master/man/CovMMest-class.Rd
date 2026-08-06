\name{CovMMest-class}
\Rdversion{1.1}
\docType{class}
\alias{CovMMest-class}

\title{ MM Estimates of Multivariate Location and Scatter }
\description{ This class, derived from the virtual class \code{"CovRobust"}
    accomodates MM Estimates of multivariate location and scatter.}
\section{Objects from the Class}{
    Objects can be created by calls of the form \code{new("CovMMest", ...)},
    but the usual way of creating \code{CovSest} objects is a call to the function
    \code{CovMMest} which serves as a constructor.
}
\section{Slots}{
    \describe{
        \item{\code{det}, \code{flag}, \code{iter}, \code{crit}:}{from the
            \code{"\linkS4class{CovRobust}"} class.
        }
        \item{c1}{tuning parameter of the loss function for MM-estimation 
        (depend on control parameters \code{eff} and \code{eff.shape}). 
        Can be computed by the internal function 
        \code{.csolve.bw.MM(p, eff, eff.shape=TRUE)}.
        For the tuning parameters of the underlying S-estimate see the slot \code{sest} and
        \code{"\linkS4class{CovSest}"}.
        }
        \item{sest}{an \code{CovSest} object containing the initial S-estimate.
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
No methods defined with class "CovMMest" in the signature.
}
\references{
  Todorov V & Filzmoser P (2009), An Object Oriented Framework for Robust Multivariate Analysis. 
  \emph{Journal of Statistical Software}, \bold{32}(3), 1--47.
  \doi{10.18637/jss.v032.i03}.
}
\author{ Valentin Todorov \email{valentin.todorov@chello.at} }
\seealso{
   \code{\link{CovMMest}}, \code{\link{Cov-class}}, \code{\link{CovRobust-class}}
}
\examples{
showClass("CovMMest")
}
\keyword{classes}
\keyword{robust}
\keyword{multivariate}
