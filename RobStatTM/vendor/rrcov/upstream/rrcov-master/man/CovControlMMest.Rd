\name{CovControlMMest}
\alias{CovControlMMest}
\title{ Constructor function for objects of class "CovControlMMest" }
\concept{High breakdown point}
\description{
    This function will create a control object \code{CovControlMMest}
    containing the control parameters for \code{CovMMest}
}
\usage{
    CovControlMMest(bdp = 0.5, eff=0.95, maxiter = 50, sest=CovControlSest(),
        trace = FALSE, tolSolve = 1e-7)
}
\arguments{
  \item{bdp}{ a numeric value specifying the required
    breakdown point.  Allowed values are between
    0.5 and 1 and the default is 0.5
  }
    \item{eff}{a numeric value specifying the required efficiency
        for the MM estimates.  Default is \code{eff=0.95}.
    }
    \item{sest}{an \code{CovControlSest} object containing control parameters for the initial S-estimate.
    }

  \item{maxiter}{maximum number of iterations allowed
    in the computation of the MM-estimate.
    Defaults to 150.
  }
  \item{trace}{whether to print intermediate results. Default is \code{trace = FALSE}.}
  \item{tolSolve}{numeric tolerance to be used as 
            a convergence tolerance for the MM-iteration.
  }
}
%\details{}
\value{
A \code{CovControlSest} object.
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
