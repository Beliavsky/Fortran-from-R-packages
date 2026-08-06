\name{CovRobust-class}
\Rdversion{1.1}
\docType{class}
\alias{CovRobust-class}
\alias{isClassic,CovRobust-method}
\alias{getMeth,CovRobust-method}
\alias{show,CovRobust-method}
\alias{summary,CovRobust-method}
\alias{getRaw,CovRobust-method}

\title{Class "CovRobust" - virtual base class for robust
 estimates of multivariate location and scatter}
\description{  \code{CovRobust} is a virtual base class used for deriving the concrete classes
representing different robust estimates of multivariate location and scatter. Here are implemeted the
standard methods common for all robust estimates like \code{show}, \code{summary} and \code{plot}.
The derived classes can override these methods and can define new ones.}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Slots}{
    \describe{
        \item{\code{iter}:}{number of iterations used to compute the estimates}
        \item{\code{crit}:}{value of the criterion function}
        \item{\code{wt}:}{weights}
        \item{\code{call}, \code{cov}, \code{center},
            \code{n.obs}, \code{mah}, \code{method},
            \code{singularity}, \code{X}:}{from the \code{"\linkS4class{Cov}"} class.
        }
    }
}
\section{Extends}{
    Class \code{"\linkS4class{Cov}"}, directly.
}
\section{Methods}{
\describe{
\item{isClassic}{\code{signature(obj = "CovRobust")}: Will return FALSE, since this is a 'Robust' object }
\item{getMeth}{\code{signature(obj = "CovRobust")}: Return the name of the particular robust method used (as a character string) }
\item{show}{\code{signature(object = "CovRobust")}: display the object }
\item{plot}{\code{signature(x = "CovRobust")}: plot the object }
\item{getRaw}{\code{signature(obj = "CovRobust")}: Return the object with the reweighted estimates replaced by the raw ones (only relevant for CovMcd, CovMve and CovOgk) }
}
}
\references{
  Todorov V & Filzmoser P (2009), An Object Oriented Framework for Robust Multivariate Analysis. 
  \emph{Journal of Statistical Software}, \bold{32}(3), 1--47.
  \doi{10.18637/jss.v032.i03}.
}
\author{Valentin Todorov \email{valentin.todorov@chello.at}}
\seealso{
  \code{\link{Cov-class}}, \code{\link{CovMcd-class}}, \code{\link{CovMest-class}}, \code{\link{CovOgk-class}}
}
\examples{
     data(hbk)
     hbk.x <- data.matrix(hbk[, 1:3])
     cv <- CovMest(hbk.x)               # it is not possible to create an object of
                                        # class CovRobust, since it is a VIRTUAL class
     cv
     summary(cv)                        # summary method for class CovRobust
     plot(cv)                           # plot method for class CovRobust
}
\keyword{classes}
\keyword{robust}
\keyword{multivariate}
