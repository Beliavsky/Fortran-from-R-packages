\name{PcaLocantore-class}
\Rdversion{1.1}
\docType{class}
\alias{PcaLocantore-class}
\alias{getQuan,PcaLocantore-method}

\title{Class "PcaLocantore" Spherical Principal Components }
\description{
The Spherical Principal Components procedure was proposed by
Locantore et al., (1999) as a functional data analysis method.
The idea is to perform classical PCA on the the data, \
projected onto a unit sphere. The estimates of the eigenvectors are consistent
and the procedure is extremly fast. The simulations of Maronna (2005) show
that this method has very good performance.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("PcaLocantore", ...)} but the
    usual way of creating \code{PcaLocantore} objects is a call to the function
    \code{PcaLocantore} which serves as a constructor.
}
\section{Slots}{
    \describe{
        \item{\code{delta}:}{Accuracy parameter}
        \item{\code{quan}:}{Object of class \code{"numeric"} The quantile h used throughout the algorithm }
        \item{\code{call}, \code{center}, \code{scale}, \code{rank}, \code{loadings},
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
    \item{getQuan}{\code{signature(obj = "PcaLocantore")}: ... }
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
showClass("PcaLocantore")
}
\keyword{robust}
\keyword{multivariate}
