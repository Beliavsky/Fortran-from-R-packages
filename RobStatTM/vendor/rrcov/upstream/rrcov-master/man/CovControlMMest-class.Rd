\name{CovControlMMest-class}
\Rdversion{1.1}
\docType{class}
\alias{CovControlMMest-class}
\alias{restimate,CovControlMMest-method}
\title{Class 'CovControlMMest' - contains control parameters for "CovMMest"}
\description{
    This class extends the \code{CovControl} class
    and contains the control parameters for CovMMest
}
\section{Objects from the Class}{
    Objects can be created by calls of the form \code{new("CovControlMMest", ...)}
    or by calling the constructor-function \code{CovControlMMest}.
}
\section{Slots}{
     \describe{
        \item{bdp}{a numeric value specifying the required
            breakdown point.  Allowed values are between
            0.5 and 1 and the default is \code{bdp=0.5}.
        }
        \item{eff}{a numeric value specifying the required efficiency
            for the MM estimates.  Default is \code{eff=0.95}.
        }
        \item{sest}{an \code{CovControlSest} object containing control parameters for the initial S-estimate.
        }
        \item{maxiter}{maximum number of iterations allowed
            in the computation of the MM-estimate.
            Default is \code{maxiter=50}.
        }
        \item{\code{trace}, \code{tolSolve}:}{from the
            \code{"\linkS4class{CovControl}"} class. \code{tolSolve} is used as 
            a convergence tolerance for the MM-iteration.
        }
    }
}
\section{Extends}{
    Class \code{"\linkS4class{CovControl}"}, directly.
}
\section{Methods}{
\describe{
    \item{restimate}{\code{signature(obj = "CovControlMMest")}: the generic
    function \code{restimate} allowes the different methods for robust estimation to be
    used polymorphically -  this function will call \code{CovMMest} passing it the control
    object and will return the obtained \code{CovRobust} object}
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

    ## the following two statements are equivalent
    ctrl1 <- new("CovControlMMest", bdp=0.25)
    ctrl2 <- CovControlMMest(bdp=0.25)

    data(hbk)
    CovMMest(hbk, control=ctrl1)

}
\keyword{classes}
\keyword{robust}
\keyword{multivariate}
