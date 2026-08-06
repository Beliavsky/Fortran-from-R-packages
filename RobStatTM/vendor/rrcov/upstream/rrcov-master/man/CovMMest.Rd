\encoding{latin1}
\name{CovMMest}
\alias{CovMMest}
\title{ MM Estimates of Multivariate Location and Scatter }
\concept{High breakdown point}
\description{
    Computes MM-Estimates of multivariate location and scatter starting from an initial S-estimate
}
\usage{
    CovMMest(x, bdp = 0.5, eff = 0.95, eff.shape=TRUE, maxiter = 50, 
        trace = FALSE, tolSolve = 1e-7, control)
}
\arguments{
  \item{x}{ a matrix or data frame. }
  \item{bdp}{a numeric value specifying the required
    breakdown point.  Allowed values are between
    0.5 and 1 and the default is \code{bdp=0.5}.
  }
    \item{eff}{a numeric value specifying the required efficiency
        for the MM estimates.  Default is \code{eff=0.95}.
    }
    \item{eff.shape}{logical; if TRUE, eff is with regard to shape-efficiency, otherwise location-efficiency. Default is \code{eff.shape=FALSE}.
    }


  \item{maxiter}{maximum number of iterations allowed
    in the computation of the S-estimate (bisquare and Rocke type).
    Default is \code{maxiter=50}.
  }
  \item{trace}{whether to print intermediate results. Default is \code{trace = FALSE}.}
  \item{tolSolve}{numeric tolerance to be used as a convergence tolerance for the MM-iteration}
  \item{control}{ a control object (S4) of class \code{\link{CovControlMMest-class}}
    containing estimation options - same as these provided in the fucntion
    specification. If the control object is supplied, the parameters from it
    will be used. If parameters are passed also in the invocation statement, they will
    override the corresponding elements of the control object.}
}
\details{
    Computes MM-estimates of multivariate location and scatter starting from an initial S-estimate.
}
\value{
  An S4 object of class \code{\link{CovMMest-class}} which is a subclass of the
  virtual class \code{\link{CovRobust-class}}.
}
\references{
    Tatsuoka, K.S. and Tyler, D.E. (2000).  
    The uniqueness of S and M-functionals under non-elliptical distributions.
    \emph{Annals of Statistics} 28, 1219--1243

    M. Salibian-Barrera,  S. Van Aelstt and G. Willems (2006). 
    Principal components analysis based on multivariate MM-estimators with fast and robust bootstrap.
    \emph{Journal of the American Statistical Association} 101, 1198--1211.
  
    R. A. Maronna, D. Martin and V. Yohai (2006). 
    \emph{Robust Statistics: Theory and Methods}.
    Wiley, New York.

  Todorov V & Filzmoser P (2009), An Object Oriented Framework for Robust Multivariate Analysis. 
  \emph{Journal of Statistical Software}, \bold{32}(3), 1--47.
  \doi{10.18637/jss.v032.i03}.
}
\author{ Valentin Todorov \email{valentin.todorov@chello.at}}
%\note{}
%\seealso{}
\examples{

library(rrcov)
data(hbk)
hbk.x <- data.matrix(hbk[, 1:3])
CovMMest(hbk.x)

## the following four statements are equivalent
c0 <- CovMMest(hbk.x)
c1 <- CovMMest(hbk.x, bdp = 0.25)
c2 <- CovMMest(hbk.x, control = CovControlMMest(bdp = 0.25))
c3 <- CovMMest(hbk.x, control = new("CovControlMMest", bdp = 0.25))

## direct specification overrides control one:
c4 <- CovMMest(hbk.x, bdp = 0.40,
             control = CovControlMMest(bdp = 0.25))
c1
summary(c1)
plot(c1)

## Deterministic MM-estmates
CovMMest(hbk.x, control=CovControlMMest(sest=CovControlSest(method="sdet")))

}
\keyword{robust}
\keyword{multivariate}
