\encoding{latin1}
\name{CovSde}
\Rdversion{1.1}
\alias{CovSde}
%
\title{Stahel-Donoho Estimates of Multivariate Location and Scatter}
\description{
    Compute a robust estimate of location and scale using the Stahel-Donoho projection based estimator
}
\usage{
CovSde(x, nsamp, maxres, tune = 0.95, eps = 0.5, prob = 0.99, 
seed = NULL, trace = FALSE, control)
}
\arguments{
  \item{x}{ a matrix or data frame.
}
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
    Defaults to \code{0.95}}
  \item{prob}{a numeric value between 0 and 1 specifying the probability of high breakdown point;
    used to compute \code{nsamp} when \code{nsamp} is omitted. Defaults to \code{0.99}.}
  \item{eps}{a numeric value between 0 and 0.5 specifying the breakdown point; used to compute
    \code{nsamp} when \code{nresamp} is omitted. Defaults to \code{0.5}. }
  \item{seed}{starting value for random generator. Default is \code{seed = NULL}.}
  \item{trace}{whether to print intermediate results. Default is \code{trace = FALSE}.  }
  \item{control}{ a control object (S4) of class \code{\link{CovControlSde-class}}
    containing estimation options - same as these provided in the fucntion
    specification. If the control object is supplied, the parameters from it
    will be used. If parameters are passed also in the invocation statement, they will
    override the corresponding elements of the control object.}
}

\details{
The projection based Stahel-Donoho estimator posses very good statistical properties, 
but it can be very slow if the number of variables is too large. It is recommended to use
this estimator if \code{n <= 1000} and \code{p<=10} or \code{n <= 5000} and \code{p<=5}.
The number of subsamples required is calculated to provide a breakdown point of 
\code{eps} with probability \code{prob} and can reach values larger than 
the larger integer value - in such case it is limited to \code{.Machine$integer.max}.
Of course you could provide \code{nsamp} in the call, i.e. \code{nsamp=1000} but 
this will not guarantee the required breakdown point of th eestimator.
For larger data sets it is better to use \code{\link{CovMcd}} or \code{\link{CovOgk}}. 
If you use \code{\link{CovRobust}}, the estimator will be selected automatically 
according on the size of the data set.}
\value{
  An S4 object of class \code{\link{CovSde-class}} which is a subclass of the
  virtual class \code{\link{CovRobust-class}}.
}

\references{
  R. A. Maronna and V.J. Yohai (1995) The Behavior of the Stahel-Donoho Robust Multivariate
  Estimator. \emph{Journal of the American Statistical Association} \bold{90} (429), 330--341.

  R. A. Maronna, D. Martin and V. Yohai (2006). \emph{Robust Statistics: Theory and Methods}.
  Wiley, New York.

  Todorov V & Filzmoser P (2009), An Object Oriented Framework for Robust Multivariate Analysis. 
  \emph{Journal of Statistical Software}, \bold{32}(3), 1--47.
  \doi{10.18637/jss.v032.i03}.
}

\author{
Valentin Todorov \email{valentin.todorov@chello.at} and
Kjell Konis \email{kjell.konis@epfl.ch}
}

\note{
The Fortran code for the Stahel-Donoho method was taken almost with no changes from
package \code{robust} which in turn has it from the \emph{Insightful Robust Library}
(thanks to by Kjell Konis).
}

\examples{
data(hbk)
hbk.x <- data.matrix(hbk[, 1:3])
CovSde(hbk.x)

## the following four statements are equivalent
c0 <- CovSde(hbk.x)
c1 <- CovSde(hbk.x, nsamp=2000)
c2 <- CovSde(hbk.x, control = CovControlSde(nsamp=2000))
c3 <- CovSde(hbk.x, control = new("CovControlSde", nsamp=2000))

## direct specification overrides control one:
c4 <- CovSde(hbk.x, nsamp=100,
             control = CovControlSde(nsamp=2000))
c1
summary(c1)
plot(c1)

## Use the function CovRobust() - if no estimation method is
##  specified, for small data sets CovSde() will be called
cr <- CovRobust(hbk.x)
cr

}
\keyword{robust}
\keyword{multivariate}
