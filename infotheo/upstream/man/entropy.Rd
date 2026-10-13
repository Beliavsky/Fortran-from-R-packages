\name{entropy}
\alias{entropy}

\title{entropy computation}
\usage{entropy(X, method="emp")}
\arguments{
  \item{X}{data.frame denoting a random vector  where columns contain variables/features and rows contain outcomes/samples.}
  \item{method}{The name of the entropy estimator. The package implements four estimators : 
			    "emp", "mm", "shrink", "sg" (default:"emp") - see details. 
               These estimators require discrete data values - see \code{\link{discretize}}.}
}
\value{ \code{entropy} returns the entropy of the data in nats.}
\description{
	  \code{entropy} takes the dataset as input and computes the 
	  entropy according 
	  to the entropy estimator \code{method}. 
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

 J. Beirlant, E. J. Dudewica, L. Gyofi, and E. van der Meulen (1997). Nonparametric 
 entropy estimation : An overview. Journal of Statistics.

 Hausser J. (2006). Improving entropy estimation and the inference of genetic regulatory networks.
 Master thesis of the National Institute of Applied Sciences of Lyon.
}
\seealso{\code{\link{condentropy}}, \code{\link{mutinformation}}, \code{\link{natstobits}}}
\examples{
  data(USArrests)
  H <- entropy(discretize(USArrests),method="shrink")
}
\keyword{misc}
