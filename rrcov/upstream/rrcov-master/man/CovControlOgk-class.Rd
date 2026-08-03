\name{CovControlOgk-class}
\Rdversion{1.1}
\docType{class}
\alias{CovControlOgk-class}
\alias{restimate,CovControlOgk-method}

\title{ Class 'CovControlOgk' - contains control parameters for CovOgk }
\description{
    This class extends the \code{CovControl} class
    and contains the control parameters for \code{"CovOgk"}
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("CovControlOgk", ...)}
or by calling the constructor-function \code{CovControlOgk}.
}
\section{Slots}{
    \describe{
        \item{niter}{ number of iterations, usually 1 or 2 since iterations
            beyond the second do not lead to improvement. }
        \item{beta}{ coverage parameter for the final reweighted estimate }
        \item{mrob}{ function for computing the robust univariate location
            and dispersion - defaults to the \code{tau scale} defined in
            Yohai and Zamar (1998) }
        \item{vrob}{ function for computing robust estimate
            of covariance between two random vectors - defaults the one proposed by
            Gnanadesikan and Kettenring (1972) }
        \item{smrob}{A string indicating the name of the function for computing
            the robust univariate location and dispersion - defaults to \code{scaleTau2} -
            the scale 'tau' function defined in Yohai and Zamar (1998)}
        \item{svrob}{A string indicating the name of the function for computing
            robust estimate of covariance between two random vectors -
            defaults to \code{gk}, the one proposed by Gnanadesikan and Kettenring (1972).
        }
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
    \item{restimate}{\code{signature(obj = "CovControlOgk")}: the generic
    function \code{restimate} allowes the different methods for robust estimation to be
    used polymorphically -  this function will call \code{CovOgk} passing it the control
    object and will return the obtained \code{CovRobust} object}
}}
\references{
  Todorov V & Filzmoser P (2009), An Object Oriented Framework for Robust Multivariate Analysis. 
  \emph{Journal of Statistical Software}, \bold{32}(3), 1--47.
  \doi{10.18637/jss.v032.i03}.
}
\author{Valentin Todorov \email{valentin.todorov@chello.at}}
\examples{
    ## the following two statements are equivalent
    ctrl1 <- new("CovControlOgk", beta=0.95)
    ctrl2 <- CovControlOgk(beta=0.95)

    data(hbk)
    CovOgk(hbk, control=ctrl1)
}
\keyword{classes}
\keyword{robust}
\keyword{multivariate}
