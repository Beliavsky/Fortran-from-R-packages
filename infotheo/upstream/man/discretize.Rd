\name{discretize}
\alias{discretize}

\title{Unsupervized Data Discretization}
\usage{discretize( X, disc="equalfreq", nbins=NROW(X)^(1/3) )}
\arguments{
  \item{X}{ A data.frame containing data to be discretized. The columns contains variables and the rows samples.}
  \item{disc}{ The name of the discretization method to be used :"equalfreq", "equalwidth" or "globalequalwidth" (default : "equalfreq") - see references.}
  \item{nbins}{ Integer specifying the number of bins to be used for the discretization. By default the number of bins is set to 
			  \eqn{{N}^(1/3)}{(N)^(1/3)} where N is the number of samples.}
}
\value{\code{discretize} returns the discretized dataset.}
\description{
      \code{discretize} discretizes \code{data} using the equal frequencies or equal width binning algorithm. 
      "equalwidth" and "equalfreq" discretizes each random variable (each column) of the data into \code{nbins}. 
      "globalequalwidth" discretizes the range of the random vector \code{data} into \code{nbins}. 
}
\references{
Meyer,  P. E.  (2008). Information-Theoretic Variable Selection and Network Inference from Microarray Data. PhD thesis of the Universite Libre de Bruxelles.

Dougherty, J., Kohavi, R., and Sahami, M. (1995). Supervised and unsupervised
discretization of continuous features. In International Conference on Machine Learning.

Yang, Y. and Webb, G. I. (2003). Discretization for naive-bayes learning: managing
discretization bias and variance. Technical Report 2003/131 School of Computer Science
and Software Engineering, Monash University.
}
\author{
  Patrick E. Meyer, Frederic Lafitte, Gianluca Bontempi, Korbinian Strimmer
}
\examples{
data(USArrests)
nbins<- sqrt(NROW(USArrests))
ew.data <- discretize(USArrests,"equalwidth", nbins)
ef.data <- discretize(USArrests,"equalfreq", nbins)
gew.data <- discretize(USArrests,"globalequalwidth", nbins)
}
\keyword{misc}
