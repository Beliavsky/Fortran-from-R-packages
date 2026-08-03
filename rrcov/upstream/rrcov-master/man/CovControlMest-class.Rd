\name{CovControlMest-class}
\docType{class}
\alias{CovControlMest-class}
\title{Class 'CovControlMest' - contains control parameters for "CovMest"}
\description{
    This class extends the \code{CovControl} class
    and contains the control parameters for CovMest}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("CovControlMest", ...)}
or by calling the constructor-function \code{CovControlMest}.
}

\section{Slots}{
 \describe{
    \item{\code{r}:}{a numeric value specifying the required
    breakdown point.  Allowed values are between
    \code{(n - p)/(2 * n)} and 1 and the default is 0.45}
    \item{\code{arp}:}{a numeric value specifying the asympthotic
    rejection point, i.e. the fraction of points receiving zero
    weight (see Rocke (1996)).  Default is \code{0.05} }
    \item{\code{eps}:}{a numeric value specifying the
    relative precision of the solution of the M-estimate.
    Defaults to \code{1e-3}}
    \item{\code{maxiter}:}{maximum number of iterations allowed
    in the computation of the M-estimate. Defaults to 120 }
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
    \item{restimate}{\code{signature(obj = "CovControlMest")}: the generic
    function \code{restimate} allowes the different methods for robust estimation to be
    used polymorphically -  this function will call \code{CovMest} passing it the control
    object and will return the obtained \code{CovRobust} object}
}

}
\references{
  Todorov V & Filzmoser P (2009), An Object Oriented Framework for Robust Multivariate Analysis. 
  \emph{Journal of Statistical Software}, \bold{32}(3), 1--47.
  \doi{10.18637/jss.v032.i03}.
}
\author{Valentin Todorov \email{valentin.todorov@chello.at}}
\examples{
    ## the following two statements are equivalent
    ctrl1 <- new("CovControlMest", r=0.4)
    ctrl2 <- CovControlMest(r=0.4)

    data(hbk)
    CovMest(hbk, control=ctrl1)

}
\keyword{classes}
\keyword{robust}
\keyword{multivariate}
