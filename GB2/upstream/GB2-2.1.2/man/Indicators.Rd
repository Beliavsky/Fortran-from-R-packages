\name{Indicators}
\Rdversion{2.1}
\alias{Indicators}
\alias{arpt.gb2}
\alias{arpr.gb2}
\alias{rmpg.gb2}
\alias{qsr.gb2}
\alias{main.gb2}
\alias{main2.gb2}
\title{Monetary Laeken Indicators under the GB2
}
\description{
Functions to calculate four primary social welfare indicators under the GB2, i.e. the at-risk-of-poverty threshold,
 the at-risk-of-poverty rate, the relative median at-risk-of-poverty gap, and the income quintile share ratio.
}
\usage{
arpt.gb2(prop, shape1, scale, shape2, shape3)
arpr.gb2(prop, shape1, shape2, shape3)
rmpg.gb2(arpr, shape1, shape2, shape3)
qsr.gb2(shape1, shape2, shape3)
main.gb2(prop, shape1, scale, shape2, shape3)
main2.gb2(prop, shape1, scale, shape12, shape13)
}
\arguments{
  \item{prop}{numeric; proportion (in general is set to 0.6).}
  \item{arpr}{numeric; the value of the at-risk-of-poverty rate.}
  \item{shape1}{numeric; positive parameter.}
  \item{scale}{numeric; positive parameter.}
  \item{shape2, shape3}{numeric; positive parameters of the Beta distribution.}
  \item{shape12}{numeric; the product of the two parameters \code{shape1} and \code{shape2}.}
  \item{shape13}{numeric; the product of the two parameters \code{shape1} and \code{shape3}.}
}
\details{
In June 2006, the  Social Protection Committee, which is a group of officials of the European Commisiion, adopts a set of common indicators for the social protection and social inclusion process.
It consists of a portfolio of 14 overarching indicators (+11 context indicators) meant to reflect the overarching objectives (a) "social cohesion" and (b) "interaction with the Lisbon strategy for growth and jobs (launched in 2000) objectives";
and of three strand portfolios for social inclusion, pensions, and health and long-term care.

The at-risk-of-poverty threshold (or ARPT) is defined as 60\% of the median national equivalized income. 

The at-risk-of-poverty rate (or ARPR) is defined as the share of persons with an equivalised disposable income below the ARPT.

The relative median at-risk-of-poverty gap (or RMPG) is defined as the difference between the median equivalised income of persons
below the ARPT and the ARPT itself, expressed as a percentage of the ARPT.

The income quintile share ratio (or QSR) is defined as the ratio of total income received by the 20\% of the country's
population with the highest income (top quintile) to that received by the 20\% of the country's population with the lowest
income (lowest quintile).

Let \eqn{x_{0.5}}{x_0.5} be the median of a GB2 with parameters \code{shape1} \eqn{= a}, \code{scale} \eqn{= b},  \code{shape2} \eqn{= p} and \code{shape3} \eqn{= q}. Then, 
\deqn{ARPT(a,b,p,q)=0.6 x_{0.5}}{ARPT(a,b,p,q)=0.6 x_0.5.}

The ARPR being scale-free, \eqn{b} can be chosen arbitrarily and can be fixed to 1.

The QSR is calculated with the help of the incomplete moments of order 1.

\code{main.gb2} and \code{main2.gb2} return a vector containing the following set of GB2 indicators: the median, the mean, the ARPR, the RMPG, the QSR and the Gini coefficient.  
The only difference is in the input parameters.
}
\value{
\code{arpt.gb2} gives the ARPT, \code{arpr.gb2} the ARPR, \code{rmpg.gb2} the RMPG, and \code{qsr.gb2}
calculates the QSR. \code{main.gb2} returns a vector containing the median of the distribution, the mean of the distribution, the ARPR, the RMPG, the QSR and the Gini coefficient. 
\code{main2.gb2} produces the same output as \code{main.gb2}.
}
\references{
\url{https://ec.europa.eu/social/main.jsp?langId=en&catId=750}
}
\author{
Monique Graf
}
\seealso{
\code{\link{qgb2}}, \code{\link{incompl.gb2}}
}
\examples{
a <- 3.9
b <- 18873
p <- 0.97
q <- 1.03
ap <- a*p
aq <- a*q

arpt <- arpt.gb2(0.6, a, b, p, q)
arpr <- arpr.gb2(0.6, a, p, q)
rmpg <- rmpg.gb2(arpr, a, p, q)
qsr <- qsr.gb2(a, p, q)

ind1 <- main.gb2(0.6, a, b, p, q)
ind2 <- main2.gb2(0.6, a, b, ap, aq)
}
\keyword{distribution}

