\name{mutinformation}
\alias{mutinformation}

\title{mutual information computation}
\usage{mutinformation(X, Y, method="emp")}
\arguments{
  \item{X}{vector/factor denoting a random variable or a data.frame denoting a random vector where columns contain variables/features and rows contain outcomes/samples.}
  \item{Y}{another random variable or random vector (vector/factor or data.frame).}
  \item{method}{The name of the entropy estimator. The package implements four estimators : 
			    "emp", "mm", "shrink", "sg" (default:"emp") - see details. 
               These estimators require discrete data values - see \code{\link{discretize}}.}
}
\value{ \code{mutinformation} returns the mutual information I(X;Y) in nats.}
\description{
	  \code{mutinformation} takes two random variables as input and computes the 
	   mutual information in nats according to the entropy estimator \code{method}.
	  If Y is not supplied and X is a matrix-like argument, the function returns a matrix of mutual 
	  information between all pairs of variables in the dataset X.
 
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
\seealso{\code{\link{condinformation}}, \code{\link{multiinformation}}, \code{\link{interinformation}}, \code{\link{natstobits}}}
\examples{
  data(USArrests)
  dat<-discretize(USArrests)
  #computes the MIM (mutual information matrix)
  I <- mutinformation(dat,method= "emp")
  I2<- mutinformation(dat[,1],dat[,2])
}
\keyword{misc}
