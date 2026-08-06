\name{PcaProj-class}
\Rdversion{1.1}
\docType{class}
\alias{PcaProj-class}
\alias{getQuan,PcaProj-method}

\title{Class "PcaProj" - Robust PCA using PP - Croux and Ruiz-Gazen (2005) algorithm }
\description{
    Holds the results of an approximation of the PP-estimators for PCA by a fast and simple algorithm: Croux and Ruiz-Gazen (2005) algorithm.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("PcaProj", ...)} but the
    usual way of creating \code{PcaProj} objects is a call to the function
    \code{PcaProj()} which serves as a constructor.
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
    \item{getQuan}{\code{signature(obj = "PcaProj")}: ... }
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
showClass("PcaProj")
}
\keyword{robust}
\keyword{multivariate}
