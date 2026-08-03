\name{CovControlMrcd}
\alias{CovControlMrcd}
%
\title{ Constructor function for objects of class "CovControlMrcd"}
\description{
  This function will create a control object \code{CovControlMrcd}
  containing the control parameters for \code{CovMrcd}
}
\usage{
CovControlMrcd(alpha = 0.5, h=NULL, maxcsteps=200, rho=NULL, 
    target=c("identity", "equicorrelation"), maxcond=50,
    trace= FALSE)
}

\arguments{
  \item{alpha}{numeric parameter controlling the size of the subsets
    over which the determinant is minimized, i.e., \code{alpha*n}
    observations are used for computing the determinant.  Allowed values
    are between 0.5 and 1 and the default is 0.5.}
  \item{h}{the size of the subset (can be between ceiling(n/2) and n). 
    Normally NULL and then it \code{h} will be calculated as 
    \code{h=ceiling(alpha*n)}. If \code{h} is provided, \code{alpha} 
    will be calculated as \code{alpha=h/n}.}
  \item{maxcsteps}{maximal number of concentration steps in the
    deterministic MCD; should not be reached.}
  \item{rho}{regularization parameter. Normally NULL and will be estimated from the data.}
  \item{target}{structure of the robust positive definite target matrix: 
  a) "identity": target matrix is diagonal matrix with robustly estimated 
  univariate scales on the diagonal or b) "equicorrelation": non-diagonal 
  target matrix that incorporates an equicorrelation structure 
  (see (17) in paper). Default is \code{target="identity"}}
  \item{maxcond}{maximum condition number allowed 
  (see step 3.4 in algorithm 1). Default is \code{maxcond=50}}
  \item{trace}{whether to print intermediate results. Default is \code{trace = FALSE}}
}
%\details{}
\value{
    A \code{CovControlMrcd} object
}
\references{
  Todorov V & Filzmoser P (2009), An Object Oriented Framework for Robust Multivariate Analysis. 
  \emph{Journal of Statistical Software}, \bold{32}(3), 1--47.
  \doi{10.18637/jss.v032.i03}.
}
\author{Valentin Todorov \email{valentin.todorov@chello.at}}
\examples{
    ## the following two statements are equivalent
    ctrl1 <- new("CovControlMrcd", alpha=0.75)
    ctrl2 <- CovControlMrcd(alpha=0.75)

    data(hbk)
    CovMrcd(hbk, control=ctrl1)
}
\keyword{classes}
\keyword{robust}
\keyword{multivariate}
