\name{CovControlMcd}
\alias{CovControlMcd}
%
\title{ Constructor function for objects of class "CovControlMcd"}
\description{
  This function will create a control object \code{CovControlMcd}
  containing the control parameters for \code{CovMcd}
}
\usage{
CovControlMcd(alpha = 0.5, nsamp = 500, scalefn=NULL, maxcsteps=200, 
seed = NULL, trace= FALSE, use.correction = TRUE)
}

\arguments{
  \item{alpha}{numeric parameter controlling the size of the subsets
    over which the determinant is minimized, i.e., \code{alpha*n}
    observations are used for computing the determinant.  Allowed values
    are between 0.5 and 1 and the default is 0.5.}
  \item{nsamp}{ number of subsets used for initial estimates or \code{"best"},
    \code{"exact"} or \code{"deterministic"}.  Default is \code{nsamp = 500}.  
    For \code{nsamp="best"} exhaustive enumeration is done, as long as the
    number of trials does not exceed 5000. For \code{"exact"},
    exhaustive enumeration will be attempted however many samples are
    needed.  In this case a warning message will be displayed saying
    that the computation can take a very long time.

    For \code{"deterministic"}, the \emph{deterministic} MCD is computed; as
    proposed by Hubert et al. (2012) it starts from the \eqn{h} most
    central observations of \emph{six} (deterministic) estimators.   
    }
  \item{scalefn}{\code{\link{function}} to compute a robust scale
    estimate or character string specifying a rule determining such a
    function, see \code{\link[robustbase]{rrcov.control}}.}
  \item{maxcsteps}{maximal number of concentration steps in the
    deterministic MCD; should not be reached.}
  \item{seed}{starting value for random generator. Default is \code{seed = NULL}}
  \item{trace}{whether to print intermediate results. Default is \code{trace = FALSE}}
  \item{use.correction}{ whether to use finite sample correction factors.
    Default is \code{use.correction=TRUE}}
}
%\details{}
\value{
A \code{CovControlMcd} object
}
\references{
  Todorov V & Filzmoser P (2009), An Object Oriented Framework for Robust Multivariate Analysis. 
  \emph{Journal of Statistical Software}, \bold{32}(3), 1--47.
  \doi{10.18637/jss.v032.i03}.
}
\author{Valentin Todorov \email{valentin.todorov@chello.at}}
\examples{
    ## the following two statements are equivalent
    ctrl1 <- new("CovControlMcd", alpha=0.75)
    ctrl2 <- CovControlMcd(alpha=0.75)

    data(hbk)
    CovMcd(hbk, control=ctrl1)
}
\keyword{classes}
\keyword{robust}
\keyword{multivariate}
