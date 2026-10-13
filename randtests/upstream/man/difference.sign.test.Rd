\name{difference.sign.test}
\alias{difference.sign.test}
%
\title{Difference Sign Test}
%%  ~~function to do ... ~~

\description{
Performs the nonparametric Difference-sign test of randomness.
}
\usage{
difference.sign.test(x, alternative)
}

\arguments{
  \item{x}{a numeric vector containing the data}
  \item{alternative}{a character string specifying the alternative hypothesis. Must be one of "\code{two.sided}" (default), "\code{left.sided}" or "\code{right.sided}".}
}
\details{
Consecutive equal values are eliminated.

The possible values "\code{two.sided}", "\code{left.sided}" and "\code{right.sided}" define the alternative hypothesis.
By using the alternative "\code{two.sided}" the null hypothesis of randomness is tested against either an increasing or decreasing trend. By using the alternative "\code{left.sided}" the null hypothesis of randomness is tested against an decreasing trend. By using the alternative "\code{right.sided}" the null hypothesis of randomness is tested against an increasing trend 
}
\value{
A list with class "htest" containing the components:
\item{statistic}{the (normalized) value of the statistic test.}
\item{parameter}{the size of the data, after the remotion of consecutive duplicate values.}
\item{p.value}{the p-value of the test.}
\item{alternative}{a character string describing the alternative hypothesis.}
\item{method}{a character string indicating the test performed.}
\item{data.name}{a character string giving the name of the data.}
\item{ds}{the total number of positive diferences (not shown on screen).}
\item{mu}{the mean value of the statistic DS (not shown on screen).}
\item{var}{the variance of the statistic DS (not shown on screen).}
}
\references{
Brockwell, P.J. and Davis, R.A. (2002). \emph{Introduction to Time Series and Forecasting}, 2nd edition, Springer (p. 37).

Mateus, A. and Caeiro, F. (2013). Comparing several tests of randomness based on the difference of observations. In T. Simos, G. Psihoyios and Ch. Tsitouras (eds.), \emph{AIP Conf. Proc.} \bold{1558}, 809--812.

Moore, G. H. and Wallis, W. A. (1943). Time Series Significance
Tests Based on Signs of Differences, \emph{Journal of the
American Statistical Association}, \bold{38}, 153--154.
}
\author{
Ayana Mateus and Frederico Caeiro
}
\examples{
##
## Example 1
## Annual Canadian Lynx trappings 1821-1934 in Canada.
## Available in datasets package
##
\dontrun{plot(lynx)}
difference.sign.test(lynx)

##
## Example 2
## Sweet potato production, harvested in the United States, between 1868 and 1937.
## Available in this package.
##
data(sweetpotato)
difference.sign.test(sweetpotato$production)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ difference sign test }
\keyword{ randomness test }% __ONLY ONE__ keyword per line
