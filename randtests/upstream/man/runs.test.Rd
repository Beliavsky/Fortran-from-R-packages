\name{runs.test}
\alias{runs.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Wald-Wolfowitz Runs Test
}
\description{
Performs the Wald-Wolfowitz runs test of randomness for continuous data.
}
\usage{
runs.test(x, alternative, threshold, pvalue, plot)
}
\arguments{
  \item{x}{a numeric vector containing the observations}
  \item{alternative}{a character string with the alternative hypothesis. Must be one of "\code{two.sided}" (default), "\code{left.sided}" or "\code{right.sided}". You can specify just the initial letter.}
  \item{threshold}{the cut-point to transform the data into a dichotomous vector}
  \item{pvalue}{a character string specifying the method used to compute the p-value. Must be one of normal (default), or exact.} 
  \item{plot}{a logic value to select whether a plot should be created. If 'TRUE', then the graph will be plotted.}
}
\details{
Data is transformed into a dichotomous vector according as each values is above or below a given \code{threshold}. Values equal to the level are removed from the sample. 

The default \code{threshold} value used in applications is the sample median which give us the special case of this test with \eqn{n_1=n_2}{n1 = n2}, the runs test above and below the median.

The possible \code{alternative} values are "\code{two.sided}", "\code{left.sided}" and "\code{right.sided}" define the alternative hypothesis. By using the alternative "\code{left.sided}" the null of randomness is tested against a trend. By using the alternative "\code{right.sided}" the null hypothesis of randomness is tested against a first order negative serial correlation.
}
\value{
A list with class "htest" containing the components:
\item{statistic}{the value of the normalized statistic test.}
\item{parameter}{a vector with the sample size, and the values of \eqn{n_1}{n1} and \eqn{n_2}{n2}.}
\item{p.value}{the p-value of the test.}
\item{alternative}{a character string describing the alternative hypothesis.}
\item{method}{a character string indicating the test performed.}
\item{data.name}{a character string giving the name of the data.}
\item{runs}{the total number of runs (not shown on screen).}
\item{mu}{the mean value of the statistic test (not shown on screen).}
\item{var}{the variance of the statistic test (not shown on screen).}
}
\references{
Brownlee, K. A. (1965). \emph{Statistical Theory and Methodology in Science and Engineering}, 2nd ed. New York: Wiley.

Gibbons, J.D. and Chakraborti, S. (2003). \emph{Nonparametric Statistical Inference}, 4th ed. (pp. 78--86).
URL: \url{http://books.google.pt/books?id=dPhtioXwI9cC&lpg=PA97&ots=ZGaQCmuEUq}

Wald, A. and Wolfowitz, J. (1940). On a test whether two samples are from the same population, \emph{The Annals of Mathematical Statistics} \bold{11}, 147--162. doi:10.1214/aoms/1177731909. \url{https://projecteuclid.org/journals/annals-of-mathematical-statistics/volume-11/issue-2/On-a-Test-Whether-Two-Samples-are-from-the-Same/10.1214/aoms/1177731909.full}
}
\author{Frederico Caeiro}
\examples{
##
## Example 1
## Data from example in Brownlee (1965), p. 223.
## Results of 23 determinations, ordered in time, of the density of the earth.
##
earthden <- c(5.36, 5.29, 5.58, 5.65, 5.57, 5.53, 5.62, 5.29, 5.44, 5.34, 5.79, 
5.10, 5.27, 5.39, 5.42, 5.47, 5.63, 5.34, 5.46, 5.30, 5.75, 5.68, 5.85)
runs.test(earthden)


##
## Example 2
## Sweet potato yield per acre, harvested in the United States, between 1868 and 1937.
## Data available in this package.
##
data(sweetpotato)
runs.test(sweetpotato$yield)
}
%\keyword{ runs test }
\keyword{ randomness test }% __ONLY ONE__ keyword per line
