\name{interinformation}
\alias{interinformation}

\title{interaction information computation}
\usage{interinformation(X, method="emp")}
\arguments{
  \item{X}{data.frame denoting a random vector where columns contain variables/features and rows contain outcomes/samples.}
  \item{method}{The name of the entropy estimator. The package implements four estimators : 
			    "emp", "mm", "shrink", "sg" (default:"emp") - see details. 
               These estimators require discrete data values - see \code{\link{discretize}}.}
}
\value{ \code{interinformation} returns the interaction information (also called synergy or complementarity), in nats, among the random variables (columns of the data.frame).}
\description{
	  \code{interinformation} takes a dataset as input and computes the 
	  the interaction information among the random variables in the dataset using   
	  the entropy estimator \code{method}. This measure is also called synergy or complementarity.
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

Jakulin, A. and Bratko, I. (2004). Testing the significance of attribute interactions.
In Proc. of 21st International Conference on Machine Learning (ICML).

McGill, W. J. (1954). Multivariate information transmission. Psychometrika, 19.
}
\seealso{\code{\link{condinformation}}, \code{\link{multiinformation}}, \code{\link{mutinformation}}, \code{\link{natstobits}}}

\examples{
  data(USArrests)
  dat<-discretize(USArrests)
  ii <- interinformation(dat, method = "sg")
}
\keyword{misc}
