\name{PcaRobust-class}
\docType{class}
\alias{PcaRobust-class}
\title{Class "PcaRobust" is a virtual base class for all robust PCA classes }
\description{ The class \code{PcaRobust} searves as a base class for deriving all other 
classes representing the results of the robust Principal Component Analisys methods}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Slots}{
     \describe{
    \item{\code{call}:}{Object of class \code{"language"} }
    \item{\code{center}:}{Object of class \code{"vector"} the center of the data }
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
  }
}
\section{Extends}{
Class \code{"\linkS4class{Pca}"}, directly.
}
\section{Methods}{
No methods defined with class "PcaRobust" in the signature.
}
\references{
  Todorov V & Filzmoser P (2009), An Object Oriented Framework for Robust Multivariate Analysis. 
  \emph{Journal of Statistical Software}, \bold{32}(3), 1--47.
  \doi{10.18637/jss.v032.i03}.
}
\author{ Valentin Todorov \email{valentin.todorov@chello.at} }
\seealso{
     \code{\link{Pca-class}}, \code{\link{PcaClassic-class}}, 
}

\examples{
showClass("PcaRobust")
}
\keyword{robust}
\keyword{multivariate}
