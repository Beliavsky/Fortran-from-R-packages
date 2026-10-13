\name{condentropy}
\alias{condentropy}

\title{conditional entropy computation}
\usage{condentropy(X, Y=NULL, method="emp")}
\arguments{
  \item{X}{data.frame denoting a random variable or random vector where columns contain variables/features and rows contain outcomes/samples.}
  \item{Y}{data.frame denoting a conditioning random variable or random vector where columns contain variables/features and rows contain outcomes/samples.}
  \item{method}{The name of the entropy estimator. The package implements four estimators : 
			    "emp", "mm", "shrink", "sg" (default:"emp") - see details. 
               These estimators require discrete data values - see \code{\link{discretize}}.}
}
\value{ \code{condentropy} returns the conditional entropy, H(X|Y), of X given Y in nats.}
\description{
	  \code{condentropy} takes two random vectors, X and Y, as input and returns the 
	  conditional entropy, H(X|Y), in nats (base e), according to the entropy estimator \code{method}. 
	  If Y is not supplied the function returns the entropy of X - see \code{\link{entropy}}.
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
\seealso{\code{\link{entropy}}, \code{\link{mutinformation}}, \code{\link{natstobits}}}
\examples{
  data(USArrests)
  dat<-discretize(USArrests)
  H <- condentropy(dat[,1], dat[,2], method = "mm")
}
\keyword{misc}
