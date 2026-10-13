\name{multiinformation}
\alias{multiinformation}

\title{multiinformation computation}
\usage{multiinformation(X, method ="emp")}
\arguments{
  \item{X}{data.frame containing a set of random variables where columns contain variables/features and rows contain outcomes/samples.}
  \item{method}{The name of the entropy estimator. The package implements four estimators : 
			    "emp", "mm", "shrink", "sg" (default:"emp") - see details. 
               These estimators require discrete data values - see \code{\link{discretize}}.}
}
\value{ \code{multiinformation} returns the multiinformation (also called total correlation) among the variables in the dataset (in nats).}
\description{
	  \code{multiinformation} takes a dataset as input and computes the 
	  multiinformation (also called total correlation) among the random variables in the dataset.
	  The value is returned in nats using the entropy estimator \code{estimator}. 
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

Studeny, M. and Vejnarova, J. (1998). The multiinformation function as a tool for
measuring stochastic dependence. In Proceedings of the NATO Advanced Study Institute
on Learning in graphical models,
}
\seealso{\code{\link{condinformation}}, \code{\link{mutinformation}}, \code{\link{interinformation}}, \code{\link{natstobits}}}

\examples{
  data(USArrests)
  dat<-discretize(USArrests)
  M <- multiinformation(dat)
}
\keyword{misc}
