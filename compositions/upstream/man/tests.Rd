\name{NormalTests}
\alias{acompNormalLocation.test}
\title{Compositional Goodness of fit test}
\description{
  Tests for several groups of additive lognormally distributed compositions.
}
\usage{
  acompNormalLocation.test(x, g=NULL, var.equal=FALSE, paired=FALSE, 
                                R=ifelse(var.equal,999,0))
          }
\arguments{
  \item{x}{a dataset of compositions (acomp) or a list of such}
  \item{g}{a factor grouping the data, not used if x is a list already. 
  Alternatively, \code{g} can be a second compositional data set. }
  \item{var.equal}{a boolean telling wether the variance of the groups
    should be considered equal}
  \item{paired}{true if a paired test should be performed}
  \item{R}{number of replicates that should be used to compute
    p-values. 0 means comparing the likelihood statistic with the
    correponding asymptotic chisq-distribution. }
}
\value{
  A classical \code{"htest"} object
  \item{data.name}{The name of the dataset as specified}
  \item{method}{a name for the test used}
  \item{alternative}{an empty string}
  \item{replicates}{a dataset of p-value distributions under the
    Null-Hypothesis got from
    nonparametric bootstrap}
  \item{p.value}{The p.value computed for this test}
}
\details{
  The tests are based on likelihood ratio statistics. 
}
\note{
  Do not trust the p-values obtained forcing \code{var.equal=TRUE} and \code{R=0}.
  This will include soon equivalent spread tests.
}
\section{Missing Policy}{
  Up to now the tests cannot handle missings. 
  }

\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}

\seealso{
  \code{\link{fitDirichlet}},\code{\link{rDirichlet}}, \code{\link{runif.acomp}},
  \code{\link{rnorm.acomp}}, 
}

\references{
Aitchison, J. (1986) \emph{The Statistical Analysis of Compositional
    Data} Monographs on Statistics and Applied Probability. Chapman &
  Hall Ltd., London (UK). 416p.\cr
}
\examples{
x <- runif.acomp(100,4)
y <- runif.acomp(100,4)
acompNormalLocation.test(list(x,y))
}
\keyword{htest}
