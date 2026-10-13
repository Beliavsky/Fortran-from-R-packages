\name{BartelsRank}
\alias{dbartelsrank}
\alias{pbartelsrank}
%\alias{qbartelsrank}
%\alias{rbartelsrank}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Distribution of the Bartels Rank Test Statistic NM
}
\description{
Probability function, distribution function%, quantile function and random generation 
for the distribution of the Bartels Rank statistic NM, for a sample of size \eqn{n}{n}.
}
\usage{
dbartelsrank(x, n, log = FALSE)
pbartelsrank(q, n, lower.tail = TRUE, log.p = FALSE)
%qbartelsrank(p, n1, n2, lower.tail = TRUE, log.p = FALSE)
%rbartelsrank(n, n1, n2)
}
\arguments{
  \item{x, q}{a numeric vector of quantiles.}
  %q, \item{p}{a numeric vector of probabilities.}
  \item{n}{number of observations to return.}
  \item{log, log.p}{logical; if TRUE, probabilities p are given as log(p).}
  \item{lower.tail}{logical; if TRUE (default), probabilities are P[X \eqn{\le} x] otherwise, P[X > x].}  
}
%\details{
%}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%
%The NM statistic is 
%\deqn{NN=\sum_{i=1}^{n-1}(R_i-R_{i+1})^2}{NM=\sum(R_i-R_{i+1})^2}
%where \eqn{R_i=rank(X_i)}, \eqn{i=1,\dots, n}{i=1,...,n}. 
%}
\value{
\code{dbartelsrank} gives the probability function and \code{pbartelsrank} gives the distribution function.
}
\section{Warning}{
This function can use large amounts of memory and stack (and even crash R if the stack limit is exceeded) if the sample size \eqn{n} is large.
}
\seealso{
\code{\link{bartels.rank.test}} to calculate the value of the statistic NM from data.
}
\author{Frederico Caeiro}
\references{
Bartels, R. (1982). The Rank Version of von Neumann's Ratio Test for Randomness, \emph{Journal of the American Statistical Association}, \bold{77}(377), 40--46.

Gibbons, J.D. and Chakraborti, S. (2003). \emph{Nonparametric Statistical Inference}, 4th ed. (pp. 97--98).\cr 
URL: \url{http://books.google.pt/books?id=dPhtioXwI9cC&lpg=PA97&ots=ZGaQCmuEUq}
}
%\examples{
%}
\keyword{ distribution }
%\keyword{ randomness test }% __ONLY ONE__ keyword per line
