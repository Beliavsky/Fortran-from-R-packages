\name{CovControlOgk}
\alias{CovControlOgk}
%
\title{ Constructor function for objects of class "CovControlOgk" }
\description{
    This function will create a control object \code{CovControlOgk}
    containing the control parameters for \code{CovOgk}
}
\usage{
CovControlOgk(niter = 2, beta = 0.9, mrob = NULL, 
vrob = .vrobGK, smrob = "scaleTau2", svrob = "gk")
}
\arguments{
  \item{niter}{ number of iterations, usually 1 or 2 since iterations
    beyond the second do not lead to improvement. }
  \item{beta}{ coverage parameter for the final reweighted estimate }
  \item{mrob}{ function for computing the robust univariate location
    and dispersion - one could use the \code{tau scale} defined in
    Yohai and Zamar (1998), see \code{\link[robustbase]{scaleTau2}}.
    The C version of this function defined by \code{smrob} is the default.
  }
  \item{vrob}{ function for computing robust estimate
    of covariance between two random vectors - one could use the function
    proposed by Gnanadesikan and Kettenring (1972), see
    \code{\link[robustbase]{covOGK}()}. The C version of this function defined
    by \code{svrob} is the default.
  }
  \item{smrob}{ a string indicating the name of the function for computing
    the robust univariate location and dispersion - defaults to
    \code{scaleTau2}  - the scale tau function defined in Yohai and Zamar (1998)
  }
  \item{svrob}{ a string indicating the name of the function for computing
    robust estimate of covariance between two random vectors - defaults \code{gk},
    the one proposed by Gnanadesikan and Kettenring (1972)
  }
}
\details{
If the user does not specify a scale and covariance function to be used in
the computations or specifies one by using the arguments \code{smrob} and \code{svrob}
(i.e. the names of the functions as strings), a native code written in C will be called which
is by far faster than the R version.

If the arguments \code{mrob} and \code{vrob} are not NULL, the specified functions
will be used via the pure R implementation of the algorithm. This could be quite slow.
}
\value{
A \code{CovControlOgk} object
}
\references{
  Maronna, R.A. and Zamar, R.H. (2002)
  Robust estimates of location and dispersion of high-dimensional datasets;
  \emph{Technometrics} \bold{44}(4), 307--317.

  Yohai, R.A. and Zamar, R.H. (1998)
  High breakdown point estimates of regression by means of the minimization of efficient scale
  \emph{JASA} \bold{86}, 403--413.

  Gnanadesikan, R. and John R. Kettenring (1972)
  Robust estimates, residuals, and outlier detection with multiresponse data.
  \emph{Biometrics} \bold{28}, 81--124.

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
