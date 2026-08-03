\name{CovControlSest}
\alias{CovControlSest}
\title{ Constructor function for objects of class "CovControlSest" }
\concept{High breakdown point}
\description{
    This function will create a control object \code{CovControlSest}
    containing the control parameters for \code{CovSest}
}
\usage{
    CovControlSest(bdp = 0.5, arp = 0.1, eps = 1e-5, maxiter = 120,
        nsamp = 500, seed = NULL, trace = FALSE, tolSolve = 1e-14, method= "sfast")
}
\arguments{
  \item{bdp}{ a numeric value specifying the required
    breakdown point.  Allowed values are between
    \code{(n - p)/(2 * n)} and 1 and the default is 0.45
  }
  \item{arp}{a numeric value specifying the asympthotic
    rejection point (for the Rocke type S estimates),
    i.e. the fraction of points receiving zero
    weight (see Rocke (1996)).  Default is \code{0.1}
  }
  \item{eps}{a numeric value specifying the
    relative precision of the solution of the S-estimate (bisquare and Rocke type).
    Defaults to \code{1e-5}.
  }
  \item{maxiter}{maximum number of iterations allowed
    in the computation of the S-estimate (bisquare and Rocke type).
    Defaults to 120.
  }
  \item{nsamp}{the number of random subsets considered. Default is \code{nsamp = 500}.}
  \item{seed}{starting value for random generator. Default is \code{seed = NULL}.}
  \item{trace}{whether to print intermediate results. Default is \code{trace = FALSE}.}
  \item{tolSolve}{numeric tolerance to be used for inversion
    (\code{\link{solve}}) of the covariance matrix in
    \code{\link{mahalanobis}}.
  }
  \item{method}{ Which algorithm to use: 'sfast'=FAST-S or 'surreal'=SURREAL }
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
\author{Valentin Todorov \email{valentin.todorov@chello.at}}
\examples{
    ## the following two statements are equivalent
    ctrl1 <- new("CovControlSest", bdp=0.4)
    ctrl2 <- CovControlSest(bdp=0.4)

    data(hbk)
    CovSest(hbk, control=ctrl1)
}
\keyword{classes}
\keyword{robust}
\keyword{multivariate}
