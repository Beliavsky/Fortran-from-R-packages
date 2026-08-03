\name{HP.loc.test}
\Rdversion{1.1}
\alias{HP.loc.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Hallin and Paindaveine Signed-Rank Tests
}
\description{
This function implements the signed-rank location tests as suggested by Hallin and Paindaveine (2002a, 2002b).
}
\usage{
HP.loc.test(X, mu = NULL, score = "rank", angles = "tyler", 
            method = "approximation", n.perm = 1000, 
            na.action = na.fail)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
a numeric data frame or matrix.
}
  \item{mu}{
a vector indicating the hypothesized value of the location. NULL represents the origin.
}
  \item{score}{
score for the pseudo mahalanobis distance. Options are 'rank', 'sign' and 'normal' scores.
}
  \item{angles}{
which angle to use. Possible are 'tyler' for spatial sign type anlges or 'interdirections'. Note however that currently only 'tyler' is implemented.
}
  \item{method}{
defines the method used for the computation of the p-value. The possibilites are
        'approximation' or 'permutation'.
}
  \item{n.perm}{
if \code{method="permutation"} specifies this the number of replications used in the permutation procedure.
}
  \item{na.action}{
a function which indicates what should happen when the data
          contain 'NA's.  Default is to fail.
}
}
\details{
The test based on interdirections is described in Hallin and Paindaveine (2002a) and the test based on Tyler's angles is described in Hallin and Paindaveine (2002b).
The two different tests are asymptotically equivalent and
in both cases is assumed that the data comes from an elliptic distribution.
}
\value{
A list with class 'htest' containing the following components:
  \item{statistic }{the value of the Q-statistic.}
  \item{parameter}{the degrees of freedom for the Q-statistic.}
  \item{p.value}{the p-value for the test.}
  \item{null.value}{the specified hypothesized value of the location.}
  \item{alternative}{a character string with the value 'two.sided'.}
  \item{method}{a character string indicating what type of test was performed.}
  \item{data.name}{a character string giving the name of the data.}
}
\references{
\cite{Hallin, M. and Paindaveine, D. (2002a), Optimal tests for multivariate location based on interdirections and pseudo-Mahalanobis ranks, \emph{Annals of Statistics}, \bold{30},  1103--1133. } 

\cite{Hallin, M. and Paindaveine, D. (2002b), Randles' interdirections or Tyler's angles?, In Y. Dodge, Ed. \emph{Statistical data analysis based on the L1-norm and related methods},  271--282. } 
}

\author{
Klaus Nordhausen
}

\seealso{
 \code{\link{tyler.shape}}, \code{\link{spatial.sign}}
}
\examples{
X <- rmvnorm(100, c(0,0,0.1)) 
HP.loc.test(X)
HP.loc.test(X, score="s")
HP.loc.test(X, score="n")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}
\keyword{ multivariate}
\keyword{nonparametric}

