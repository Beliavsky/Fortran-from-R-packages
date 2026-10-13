\name{condinformation}
\alias{condinformation}

\title{conditional mutual information computation}
\usage{condinformation(X, Y, S=NULL, method="emp")}
\arguments{
  \item{X}{vector/factor denoting a random variable or a data.frame denoting a random vector where columns contain variables/features and rows contain outcomes/samples.}
  \item{Y}{another random variable or random vector (vector/factor or data.frame).}
  \item{S}{the conditioning random variable or random vector (vector/factor or data.frame).}
  \item{method}{The name of the entropy estimator. The package implements four estimators : 
			    "emp", "mm", "shrink", "sg" (default:"emp") - see details. 
               These estimators require discrete data values - see \code{\link{discretize}}.}
}
\value{ \code{condinformation} returns the conditional mutual information, I(X;Y|S), in nats.}
\description{
	  \code{condinformation} takes three random variables as input and computes the 
	  conditional mutual information in nats according to the entropy estimator \code{method}.
	  If S is not supplied the function returns the mutual information between X and Y - see \code{\link{mutinformation}}
}
\details{
\itemize{
      \item "emp" : This estimator computes the entropy of the empirical probability distribution.
      \item "mm" : This is the Miller-Madow asymptotic bias corrected empirical estimator.
      \item "shrink" : This is a shrinkage estimate of the entropy of a Dirichlet probability distribution.
	  \item "sg" : This is the Schurmann-Grassberger estimate of the entropy of a Dirichlet probability distribution.
}
}
\author{
  Patrick E. Meyer
}
\references{
Meyer,  P. E.  (2008). Information-Theoretic Variable Selection and Network Inference from Microarray Data. PhD thesis of the Universite Libre de Bruxelles.
 
Cover, T. M. and Thomas, J. A. (1990). Elements of Information Theory. John Wiley,
New York.
}
\seealso{\code{\link{mutinformation}}, \code{\link{multiinformation}}, \code{\link{interinformation}}, \code{\link{natstobits}}}
\examples{
  data(USArrests)
  dat<-discretize(USArrests)
  I <- condinformation(dat[,1],dat[,2],dat[,3],method="emp")
}
\keyword{misc}
