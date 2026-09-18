\name{generator.otrimle}
\alias{generator.otrimle}


\title{Generates random data from OTRIMLE output model}

\description{
  This uses data and the output of \code{\link{otrimle}} or
  \code{\link{rimle}} to generate a new artificial dataset of the size
  of the original data using noise and cluster proportions from the
  clustering output. The clusters are then generated from multivariate
  normal distributions with the parameters estimated by
  \code{\link{otrimle}}, the noise is generated resampling from what is
  estimated as moise component with weights given by posterior
  probabilities of all observations to be noise. See Hennig and Coretto
  (2021).   
} 

\usage{
  generator.otrimle(data, fit)
}


\arguments{
  \item{data}{something that can be coerced into a matrix. Dataset.}
  \item{fit}{output object of \code{\link{otrimle}} or \code{\link{rimle}}.}
}
  
\value{
  A list with components \code{data, clustering}.

  \item{data}{matrix of generated data.}
  \item{cs}{vector of integers. Clustering indicator.}
}

\references{
  Hennig, C. and P.Coretto (2021). An adequacy approach for deciding the
  number of clusters for OTRIMLE robust Gaussian mixture based
  clustering. To appear in \emph{Australian and New Zealand Journal of
  Statistics}, \url{https://arxiv.org/abs/2009.00921}.
}

\author{Christian Hennig
  \email{christian.hennig@unibo.it}
  \url{https://www.unibo.it/sitoweb/christian.hennig/en/}
}

\seealso{
  \code{\link{kerndensp}}, \code{\link{kerndensmeasure}},
  \code{\link{otrimle}},  \code{\link{rimle}}
}

\examples{
   data(banknote)
   selectdata <- c(1:30,101:110,117:136,160:161)
   set.seed(555566)
   x <- banknote[selectdata,5:7]
   ox <- otrimle(x, G=2 , ncores = 1)
   str(generator.otrimle(x, ox))
}

\keyword{cluster}
\keyword{datagen}


