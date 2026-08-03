\name{PcaGrid-class}
\Rdversion{1.1}
\docType{class}
\alias{PcaGrid-class}
\alias{getQuan,PcaGrid-method}

\title{Class "PcaGrid" - Robust PCA using PP - GRID search Algorithm  }
\description{
    Holds the results of an approximation of the PP-estimators for PCA using the grid search algorithm in the plane.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("PcaGrid", ...)} but the
    usual way of creating \code{PcaGrid} objects is a call to the function
    \code{PcaGrid()} which serves as a constructor.
}
\section{Slots}{
    \describe{
        \item{\code{call}, \code{center}, \code{scale}, \code{rank}, \code{loadings},
            \code{eigenvalues}, \code{scores}, \code{k},
            \code{sd}, \code{od}, \code{cutoff.sd}, \code{cutoff.od},
            \code{flag}, \code{n.obs}:}{
            from the \code{"\linkS4class{Pca}"} class.
        }
    }
}

\section{Extends}{
    Class \code{"\linkS4class{PcaRobust}"}, directly.
    Class \code{"\linkS4class{Pca}"}, by class \code{"\linkS4class{PcaRobust}"}, distance 2.
}
\section{Methods}{
  \describe{
    \item{getQuan}{\code{signature(obj = "PcaGrid")}: ... }
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
showClass("PcaGrid")
}
\keyword{robust}
\keyword{multivariate}
