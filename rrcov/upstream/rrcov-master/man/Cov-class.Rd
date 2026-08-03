\name{Cov-class}
\docType{class}
\alias{Cov-class}
\alias{getCenter-method}
\alias{getCenter,Cov-method}
\alias{getCov,Cov-method}
\alias{getCorr,Cov-method}
\alias{getData,Cov-method}
\alias{getDistance,Cov-method}
\alias{getEvals,Cov-method}
\alias{getDet,Cov-method}
\alias{getShape,Cov-method}
\alias{getFlag,Cov-method}
\alias{isClassic,method}
\alias{isClassic,Cov-method}
\alias{plot}
\alias{plot,Cov,missing-method}
\alias{show,Cov-method}
\alias{summary,Cov-method}
\alias{Uvector-class}
\alias{Ulist-class}
\alias{Utable-class}
\alias{Umatrix-class}
\alias{Ufunction-class}

\title{Class "Cov" -- a base class for estimates of multivariate
    location and scatter}
\description{The class \code{Cov} represents an estimate of the
    multivariate location and scatter of a data set. The objects of class \code{Cov}
    contain the classical estimates and serve as base for deriving other
    estimates, i.e. different types of robust estimates.
}
\section{Objects from the Class}{
    Objects can be created by calls of the form \code{new("Cov", ...)},
    but the usual way of creating \code{Cov} objects is a call to the function
    \code{Cov} which serves as a constructor.
}
\section{Slots}{
    \describe{
    \item{\code{call}:}{Object of class \code{"language"} }
    \item{\code{cov}:}{covariance matrix}
    \item{\code{center}:}{location}
    \item{\code{n.obs}:}{number of observations used for the computation of the estimates}
    \item{\code{mah}:}{mahalanobis distances}
    \item{\code{det}:}{determinant}
    \item{\code{flag}:}{flags (FALSE if suspected an outlier)}
    \item{\code{method}:}{a character string describing the method used to
        compute the estimate: "Classic"}
    \item{\code{singularity}:}{a list with singularity information for the
        covariance matrix (or \code{NULL} of not singular)}
    \item{\code{X}:}{data}
  }
}
\section{Methods}{
  \describe{
    \item{getCenter}{\code{signature(obj = "Cov")}: location vector }
    \item{getCov}{\code{signature(obj = "Cov")}: covariance matrix }
    \item{getCorr}{\code{signature(obj = "Cov")}: correlation matrix }
    \item{getData}{\code{signature(obj = "Cov")}: data frame }
    \item{getDistance}{\code{signature(obj = "Cov")}: distances }
    \item{getEvals}{\code{signature(obj = "Cov")}: Computes and returns
        the eigenvalues of the covariance matrix }
    \item{getDet}{\code{signature(obj = "Cov")}: Computes and returns 
    the determinant of the covariance matrix (or 0 if the covariance matrix is singular) }
    \item{getShape}{\code{signature(obj = "Cov")}: Computes and returns 
    the shape matrix corresponding to the covariance matrix (i.e. the covariance matrix scaled to have determinant =1) }
    \item{getFlag}{\code{signature(obj = "Cov")}: Flags observations as outliers if the corresponding mahalanobis distance is larger then \code{qchisq(prob, p)} where \code{prob} defaults to 0.975. }
    \item{isClassic}{\code{signature(obj = "Cov")}: returns TRUE by default. If necessary, the robust
            classes will override }
    \item{plot}{\code{signature(x = "Cov")}: plot the object }
    \item{show}{\code{signature(object = "Cov")}: display the object }
    \item{summary}{\code{signature(object = "Cov")}: calculate summary information }
     }
}
\references{
  Todorov V & Filzmoser P (2009), An Object Oriented Framework for Robust Multivariate Analysis. 
  \emph{Journal of Statistical Software}, \bold{32}(3), 1--47.
  \doi{10.18637/jss.v032.i03}.
}
\author{Valentin Todorov \email{valentin.todorov@chello.at}}
\examples{
showClass("Cov")
}
\keyword{classes}
\keyword{robust}
\keyword{multivariate}
