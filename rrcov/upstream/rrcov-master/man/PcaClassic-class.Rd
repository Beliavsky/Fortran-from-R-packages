\name{PcaClassic-class}
\docType{class}
\alias{PcaClassic-class}
\alias{getQuan,PcaClassic-method}

\title{Class "PcaClassic" - Principal Components Analysis }
\description{ Contains the results of a classical Principal Components Analysis }
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("PcaClassic", ...)} but the 
    usual way of creating \code{PcaClassic} objects is a call to the function
    \code{PcaClassic} which serves as a constructor.
}
\section{Slots}{
     \describe{
    \item{\code{call}:}{Object of class \code{"language"} }
    \item{\code{center}:}{Object of class \code{"vector"} the center of the data }
    \item{\code{scale}:}{Object of class \code{"vector"} the scaling applied to each variable }
    \item{\code{rank}:}{Object of class \code{"numeric"} the rank of the data matrix }
    \item{\code{loadings}:}{Object of class \code{"matrix"} the matrix 
    of variable loadings (i.e., a matrix whose columns contain the eigenvectors) }
    \item{\code{eigenvalues}:}{Object of class \code{"vector"} the eigenvalues }
    \item{\code{scores}:}{Object of class \code{"matrix"} the scores - the value 
    of the projected on the space of the principal components data (the centred 
    (and scaled if requested) data multiplied 
    by the \code{loadings} matrix) is returned.  Hence, \code{cov(scores)} 
    is the diagonal matrix \code{diag(eigenvalues)} }
    \item{\code{k}:}{Object of class \code{"numeric"} number of (choosen) principal components }
    \item{\code{sd}:}{Object of class \code{"Uvector"} Score distances within the robust PCA subspace }
    \item{\code{od}:}{Object of class \code{"Uvector"} Orthogonal distances to the robust PCA subspace }
    \item{\code{cutoff.sd}:}{Object of class \code{"numeric"} Cutoff value for the score distances}
    \item{\code{cutoff.od}:}{Object of class \code{"numeric"} Cutoff values for the orthogonal distances }
    \item{\code{flag}:}{Object of class \code{"Uvector"} The observations whose score distance is larger 
    than cutoff.sd or whose orthogonal distance is larger than cutoff.od can be considered 
    as outliers and receive a flag equal to zero.
    The regular observations receive a flag 1 }
    \item{\code{n.obs}:}{Object of class \code{"numeric"} the number of observations }
    \item{\code{eig0}:}{Object of class \code{"vector"} all eigenvalues }
    \item{\code{totvar0}:}{Object of class \code{"numeric"} the total variance explained (\code{=sum(eig0)}) }
  }
}
\section{Extends}{
Class \code{"\linkS4class{Pca}"}, directly.
}
\section{Methods}{
  \describe{
    \item{getQuan}{\code{signature(obj = "PcaClassic")}: returns the number of 
    observations used in the computation, i.e. n.obs }
     }
}
\references{
  Todorov V & Filzmoser P (2009), An Object Oriented Framework for Robust Multivariate Analysis. 
  \emph{Journal of Statistical Software}, \bold{32}(3), 1--47.
  \doi{10.18637/jss.v032.i03}.
}
\author{ Valentin Todorov \email{valentin.todorov@chello.at} }
\seealso{
     \code{\link{PcaRobust-class}}, \code{\link{Pca-class}}, \code{\link{PcaClassic}} 
}

\examples{
showClass("PcaClassic")
}
\keyword{robust}
\keyword{multivariate}
