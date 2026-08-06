\name{CovControlSde}
\Rdversion{1.1}
\alias{CovControlSde}
\title{ Constructor function for objects of class "CovControlSde" }
\concept{High breakdown point}
\description{
    This function will create a control object \code{CovControlSde}
    containing the control parameters for \code{CovSde}
}
\usage{
CovControlSde(nsamp = 0, maxres = 0, tune = 0.95, eps = 0.5, prob = 0.99,
    seed = NULL, trace = FALSE, tolSolve = 1e-14)
}
\arguments{
  \item{nsamp}{a positive integer giving the number of resamples required;
    \code{nsamp} may not be reached if too many of the \code{p}-subsamples,
    chosen out of the observed vectors, are in a hyperplane.
    If \code{nsamp = 0} all possible subsamples are taken.
    If \code{nsamp} is omitted, it is calculated to provide a breakdown point
    of \code{eps} with probability \code{prob}.}
  \item{maxres}{a positive integer specifying the maximum number of
    resamples to be performed including those that are discarded due to linearly
    dependent subsamples.  If \code{maxres} is omitted it will be set to 2 times \code{nsamp}.}
  \item{tune}{a numeric value between 0 and 1 giving the fraction of the data to receive non-zero weight.
  Defaults to \code{0.95}. }
  \item{prob}{a numeric value between 0 and 1 specifying the probability of high breakdown point;
    used to compute \code{nsamp} when \code{nsamp} is omitted. Defaults to \code{0.99}. }
  \item{eps}{a numeric value between 0 and 0.5 specifying the breakdown point;
    used to compute \code{nsamp} when \code{nresamp} is omitted. Defaults to \code{0.5}. }
  \item{seed}{starting value for random generator. Default is \code{seed = NULL}.}
  \item{trace}{whether to print intermediate results. Default is \code{trace = FALSE}.  }
  \item{tolSolve}{numeric tolerance to be used for inversion
    (\code{\link{solve}}) of the covariance matrix in
    \code{\link{mahalanobis}}.}
}
\value{
A \code{CovControlSde} object.
}
\references{
  Todorov V & Filzmoser P (2009), An Object Oriented Framework for Robust Multivariate Analysis. 
  \emph{Journal of Statistical Software}, \bold{32}(3), 1--47.
  \doi{10.18637/jss.v032.i03}.
}
\author{Valentin Todorov \email{valentin.todorov@chello.at}}
\examples{
    ## the following two statements are equivalent
    ctrl1 <- new("CovControlSde", nsamp=2000)
    ctrl2 <- CovControlSde(nsamp=2000)

    data(hbk)
    CovSde(hbk, control=ctrl1)
}
\keyword{classes}
\keyword{robust}
\keyword{multivariate}
