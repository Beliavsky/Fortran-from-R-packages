\name{CovControlMve-class}
\docType{class}
\alias{CovControlMve-class}
\alias{restimate,CovControlMve-method}

\title{Class 'CovControlMve' - contains control parameters for CovMve }
\description{
    This class extends the \code{CovControl} class
    and contains the control parameters for \code{"CovMve"}
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("CovControlMve", ...)}
or by calling the constructor-function \code{CovControlMve}.
}
\section{Slots}{
     \describe{
    \item{\code{alpha}:}{numeric parameter controlling the size of the subsets
        over which the determinant is minimized, i.e., \code{alpha*n}
        observations are used for computing the determinant.  Allowed values
        are between 0.5 and 1 and the default is 0.5.}
    \item{\code{nsamp}:}{ number of subsets used for initial estimates or \code{"best"}
        or \code{"exact"}.  Default is \code{nsamp = 500}.  For
        \code{nsamp="best"} exhaustive enumeration is done, as long as the
        number of trials does not exceed 5000. For \code{"exact"},
        exhaustive enumeration will be attempted however many samples are
        needed.  In this case a warning message will be displayed saying
        that the computation can take a very long time.}
    \item{\code{seed}:}{starting value for random generator. Default is \code{seed = NULL}}
    \item{\code{trace}, \code{tolSolve}:}{from the
        \code{"\linkS4class{CovControl}"} class.
    }
  }
}
\section{Extends}{
    Class \code{"\linkS4class{CovControl}"}, directly.
}
\section{Methods}{
  \describe{
    \item{restimate}{\code{signature(obj = "CovControlMve")}: the generic
    function \code{restimate} allowes the different methods for robust estimation to be
    used polymorphically -  this function will call \code{CovMve} passing it the control
    object and will return the obtained \code{CovRobust} object}
}}
\references{
  Todorov V & Filzmoser P (2009), An Object Oriented Framework for Robust Multivariate Analysis. 
  \emph{Journal of Statistical Software}, \bold{32}(3), 1--47.
  \doi{10.18637/jss.v032.i03}.
}
\author{Valentin Todorov \email{valentin.todorov@chello.at}}
\examples{
    ## the following two statements are equivalent
    ctrl1 <- new("CovControlMve", alpha=0.75)
    ctrl2 <- CovControlMve(alpha=0.75)

    data(hbk)
    CovMve(hbk, control=ctrl1)
}
\keyword{classes}
\keyword{robust}
\keyword{multivariate}
