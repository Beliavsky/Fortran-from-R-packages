\name{CovControlSest-class}
\Rdversion{1.1}
\docType{class}
\alias{CovControlSest-class}
\alias{restimate,CovControlSest-method}
\title{Class 'CovControlSest' - contains control parameters for "CovSest"}
\description{
    This class extends the \code{CovControl} class
    and contains the control parameters for CovSest
}
\section{Objects from the Class}{
    Objects can be created by calls of the form \code{new("CovControlSest", ...)}
    or by calling the constructor-function \code{CovControlSest}.
}
\section{Slots}{
     \describe{
        \item{bdp}{a numeric value specifying the required
            breakdown point.  Allowed values are between
            \code{(n - p)/(2 * n)} and 1 and the default is \code{bdp=0.45}.
        }
        \item{arp}{a numeric value specifying the asympthotic
            rejection point (for the Rocke type S estimates),
            i.e. the fraction of points receiving zero
            weight (see Rocke (1996)).  Default is \code{arp=0.1}.
        }
        \item{eps}{a numeric value specifying the
            relative precision of the solution of the S-estimate
            (bisquare and Rocke type). Default is to \code{eps=1e-5}.
        }
        \item{maxiter}{maximum number of iterations allowed
            in the computation of the S-estimate (bisquare and Rocke type).
            Default is \code{maxiter=120}.
        }
        \item{nsamp}{the number of random subsets considered.
            Default is \code{nsamp = 500}.}
        \item{seed}{starting value for random generator. Default is \code{seed = NULL}.}
        \item{method}{Which algorithm to use: 'sfast'=FAST-S, 'surreal'=Ruppert's SURREAL algorithm,
            'bisquare'=Bisquare S-estimation with HBDP start or 'rocke' for Rocke type S-estimates}
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
    \item{restimate}{\code{signature(obj = "CovControlSest")}: the generic
    function \code{restimate} allowes the different methods for robust estimation to be
    used polymorphically -  this function will call \code{CovSest} passing it the control
    object and will return the obtained \code{CovRobust} object}
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
    ctrl1 <- new("CovControlSest", bdp=0.4)
    ctrl2 <- CovControlSest(bdp=0.4)

    data(hbk)
    CovSest(hbk, control=ctrl1)

}
\keyword{classes}
\keyword{robust}
\keyword{multivariate}
