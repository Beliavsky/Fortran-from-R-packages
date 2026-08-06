\name{pca.distances}
\alias{pca.distances}
\title{
    Compute score and orthogonal distances for Principal Components (objects of class 'Pca')
}
\description{
  Compute score and orthogonal distances for an object (derived from)\code{\link{Pca-class}}.
}
\usage{
    pca.distances(obj, data, r, crit=0.975)
}

\arguments{
  \item{obj}{an object of class (derived from) \code{"Pca"}.}
  \item{data}{The data matrix for which the \code{"Pca"} object was computed.}
  \item{r}{rank of data}
  \item{crit}{Criterion to use for computing the cutoff values.}
}
\details{
This function calculates the score and orthogonal distances and the
appropriate cutoff values for identifying outlying observations. 
The computed values are used to 
create a vector a of flags, one for each observation, identifying the outliers.

}

\value{
  An S4 object of class derived from the virtual class \code{\link{Pca-class}} - 
  the same object passed to the function, but with the score and orthogonal 
  distances as well as their cutoff values and the corresponding flags appended to it. 
}

\references{
    M. Hubert, P. J. Rousseeuw, K. Vanden Branden (2005), ROBPCA: a new 
    approach to robust principal components analysis, \emph{Technometrics}, \bold{47}, 64--79.   
    
  Todorov V & Filzmoser P (2009), An Object Oriented Framework for Robust Multivariate Analysis. 
  \emph{Journal of Statistical Software}, \bold{32}(3), 1--47.
  \doi{10.18637/jss.v032.i03}.
}
\author{Valentin Todorov \email{valentin.todorov@chello.at}}
\examples{

## PCA of the Hawkins Bradu Kass's Artificial Data
##  using all 4 variables
data(hbk)
pca <- PcaHubert(hbk)
pca.distances(pca, hbk, rankMM(hbk))
}
\keyword{robust}
\keyword{multivariate}
