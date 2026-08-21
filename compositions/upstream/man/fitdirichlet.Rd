\name{fitdirichlet}
\alias{fitDirichlet}
\title{Fitting a Dirichlet distribution}
\description{ Fits a Dirichtlet Distribution to a dataset
  by maximum likelihood.  }
\usage{
fitDirichlet(x,elog=mean(ult(x)),alpha0=rep(1,length(elog)),maxIter=20,n=nrow(x))
}
\arguments{
  \item{x}{a dataset of compositions (acomp)}
  \item{elog}{the expected log can provided instead of the dataset
    itself.}
  \item{alpha0}{the start value for alpha parameter in the
    iteration}
  \item{maxIter}{The maximum number of iterations in the
    Fischer scoring method. }
  \item{n}{the number of datapoints used to estimate elog}
}
\value{
  \item{alpha}{the estimated parameter}
  \item{loglikelihood}{the likelihood}
  \item{df}{The dimension of the dataset minus the dimension of the
    parameter}
}
\details{

  The fitting is done using a modified version of the Fisher-Scoring
method using analytiscal expressions for log mean and log variance. The
modification is introducted to prevent the algorithm from leaving the
admissible parameter set. It reduced the stepsize to at most have of
distance to the limit of the admissible parameter set.
  
}
\section{Missing Policy}{
  Up to now the fitting can not handle missings. 
  }

\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}

\seealso{ \code{\link{rDirichlet}}, \code{acompDirichletGOF.test},
  \code{\link{runif.acomp}}, \code{\link{rnorm.acomp}}, }

\references{
Aitchison, J. (1986) \emph{The Statistical Analysis of Compositional
    Data} Monographs on Statistics and Applied Probability. Chapman &
  Hall Ltd., London (UK). 416p.\cr
}
\examples{
x <- rDirichlet.acomp(100,c(1,2,3,4))
fitDirichlet(x)
}
\keyword{htest}
