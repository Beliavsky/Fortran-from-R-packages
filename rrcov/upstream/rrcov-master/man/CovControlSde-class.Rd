\name{CovControlSde-class}
\Rdversion{1.1}
\docType{class}
\alias{CovControlSde-class}
\alias{restimate,CovControlSde-method}

\title{Class 'CovControlSde' - contains control parameters for "CovSde"}
\description{
    This class extends the \code{CovControl} class
    and contains the control parameters for CovSde
}
\section{Objects from the Class}{
    Objects can be created by calls of the form \code{new("CovControlSde", ...)}
    or by calling the constructor-function \code{CovControlSde}.
}
\section{Slots}{
    \describe{
        \item{\code{nsamp}}{a positive integer giving the number of resamples required}
        \item{\code{maxres}}{a positive integer specifying the maximum number of
        resamples to be performed including those that are discarded due to linearly
        dependent subsamples.}
        \item{\code{tune}}{a numeric value between 0 and 1 giving the fraction of
        the data to receive non-zero weight. Default is \code{tune = 0.95}.}
        \item{\code{prob}}{a numeric value between 0 and 1 specifying
        the probability of high breakdown point; used to compute
        \code{nsamp} when \code{nsamp} is omitted. Default is \code{prob = 0.99}.}
        \item{\code{eps}}{a numeric value between 0 and 0.5 specifying the breakdown point;
        used to compute \code{nsamp} when \code{nresamp} is omitted.
        Default is \code{eps = 0.5}.}
        \item{seed}{starting value for random generator. Default is \code{seed = NULL}.}
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
    \item{restimate}{\code{signature(obj = "CovControlSde")}: ... }
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
    ctrl1 <- new("CovControlSde", nsamp=2000)
    ctrl2 <- CovControlSde(nsamp=2000)

    data(hbk)
    CovSde(hbk, control=ctrl1)
}
\keyword{classes}
\keyword{robust}
\keyword{multivariate}
