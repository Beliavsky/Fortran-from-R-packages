\name{PcaCov-class}
\Rdversion{1.1}
\docType{class}
\alias{PcaCov-class}
\alias{getQuan,PcaCov-method}

\title{Class "PcaCov" - Robust PCA based on a robust covariance matrix }
\description{
Robust PCA are obtained by replacing the classical covariance matrix
by a robust covariance estimator. This can be one of the available
in \code{rrcov} estimators, i.e. MCD, OGK, M, S or Stahel-Donoho estimator.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("PcaCov", ...)} but the
    usual way of creating \code{PcaCov} objects is a call to the function
    \code{PcaCov} which serves as a constructor.
}
\section{Slots}{
    \describe{
        \item{\code{quan}:}{Object of class \code{"numeric"}
            The quantile \code{h} used throughout the algorithm }
        \item{\code{call}, \code{center}, \code{rank}, \code{loadings},
            \code{eigenvalues}, \code{scores}, \code{k},
            \code{sd}, \code{od}, \code{cutoff.sd}, \code{cutoff.od},
            \code{flag}, \code{n.obs}, \code{eig0}, \code{totvar0}:}{
            from the \code{"\linkS4class{Pca}"} class.
        }
    }
}

\section{Extends}{
Class \code{"\linkS4class{PcaRobust}"}, directly.
Class \code{"\linkS4class{Pca}"}, by class "PcaRobust", distance 2.
}
\section{Methods}{
  \describe{
    \item{getQuan}{\code{signature(obj = "PcaCov")}: ... }
     }
}
\references{
  Todorov V & Filzmoser P (2009), An Object Oriented Framework for Robust Multivariate Analysis. 
  \emph{Journal of Statistical Software}, \bold{32}(3), 1--47.
  \doi{10.18637/jss.v032.i03}.
}
\author{ Valentin Todorov \email{valentin.todorov@chello.at} }
\seealso{
     \code{\link{PcaRobust-class}}, \code{\link{Pca-class}}, \code{\link{PcaClassic}}, \code{\link{PcaClassic-class}}
}

\examples{
showClass("PcaCov")
}
\keyword{robust}
\keyword{multivariate}
