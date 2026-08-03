\name{CovSest-class}
\Rdversion{1.1}
\docType{class}
\alias{CovSest-class}

\title{ S Estimates of Multivariate Location and Scatter }
\description{ This class, derived from the virtual class \code{"CovRobust"}
    accomodates S Estimates of multivariate location and scatter computed
    by the  \sQuote{Fast S} or \sQuote{SURREAL} algorithm.}
\section{Objects from the Class}{
    Objects can be created by calls of the form \code{new("CovSest", ...)},
    but the usual way of creating \code{CovSest} objects is a call to the function
    \code{CovSest} which serves as a constructor.
}
\section{Slots}{
    \describe{
        \item{\code{iter}, \code{crit}, \code{wt}:}{from the
            \code{"\linkS4class{CovRobust}"} class.
        }
        \item{\code{iBest}, \code{nsteps}, \code{initHsets}:}{parameters for 
        deterministic S-estimator (the best initial subset, number of concentration 
        steps to convergence for each of the initial subsets, and the computed 
        initial subsets, respectively).
        }
        \item{\code{cc}, \code{kp}}{tuning parameters used in Tukey biweight 
        loss function, as determined by bdp. Can be computed by the internal function 
        \code{.csolve.bw.S(bdp, p)}.
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
No methods defined with class "CovSest" in the signature.
}
\references{
  Todorov V & Filzmoser P (2009), An Object Oriented Framework for Robust Multivariate Analysis. 
  \emph{Journal of Statistical Software}, \bold{32}(3), 1--47.
  \doi{10.18637/jss.v032.i03}.
}
\author{ Valentin Todorov \email{valentin.todorov@chello.at} }
\seealso{
   \code{\link{CovSest}}, \code{\link{Cov-class}}, \code{\link{CovRobust-class}}
}
\examples{
showClass("CovSest")
}
\keyword{classes}
\keyword{robust}
\keyword{multivariate}
