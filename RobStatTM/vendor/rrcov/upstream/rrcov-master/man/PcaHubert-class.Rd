\name{PcaHubert-class}
\Rdversion{1.1}
\docType{class}
\alias{PcaHubert-class}
\alias{getQuan,PcaHubert-method}
\title{Class "PcaHubert" - ROBust method for Principal Components Analysis }
\description{ The ROBPCA algorithm was proposed by Hubert et al (2005) and stays for
    'ROBust method for Principal Components Analysis'. It is resistant to
    outliers in the data. The robust loadings are computed using
    projection-pursuit techniques and the MCD method. Therefore ROBPCA
    can be applied to both low and high-dimensional data sets. In low
    dimensions, the MCD method is applied.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("PcaHubert", ...)} but the
    usual way of creating \code{PcaHubert} objects is a call to the function
    \code{PcaHubert} which serves as a constructor.
}
\section{Slots}{
    \describe{
        \item{\code{alpha}:}{Object of class \code{"numeric"} the fraction of outliers
            the algorithm should resist - this is the argument alpha}
        \item{\code{quan}:}{The quantile \code{h} used throughout the algorithm }
        \item{\code{skew}:}{Whether the adjusted outlyingness algorithm for skewed data was used}
        \item{\code{ao}:}{Object of class \code{"Uvector"} Adjusted outlyingness within the robust PCA subspace }
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
    \item{getQuan}{\code{signature(obj = "PcaHubert")}: Returns the quantile
    used throughout the algorithm }
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
showClass("PcaHubert")
}
\keyword{robust}
\keyword{multivariate}
