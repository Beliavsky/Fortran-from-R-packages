\name{fitSameMeanDifferentVarianceModel}
\alias{fitSameMeanDifferentVarianceModel}
\title{Fit Same Mean Different Variance Model}
\description{
  Fits a model of the same mean, but different variances model to a set
  of several multivariate normal groups by maximum likelihood.
}
\usage{
fitSameMeanDifferentVarianceModel(x)
          }
\arguments{
  \item{x}{list of rmult type datasets}
}
\value{
  \item{mean}{the estimated mean}
  \item{vars}{a list of estimated variance-covariance matrices}
  \item{N}{a vector containing the sizes of the groups}
}
\details{
  The function tries to fit a normal model with different variances but the
  same mean between different groups.
}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{acompNormalLocation.test}} 
}
\references{
Aitchison, J. (1986) \emph{The Statistical Analysis of Compositional
    Data} Monographs on Statistics and Applied Probability. Chapman &
  Hall Ltd., London (UK). 416p.\cr
}
\examples{
fitSameMeanDifferentVarianceModel
}
\keyword{htest}
