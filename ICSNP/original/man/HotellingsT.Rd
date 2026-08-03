\name{HotellingsT2}
\alias{HotellingsT2}
\alias{HotellingsT2.default}
\alias{HotellingsT2.formula}

\title{Hotelling's T2 Test}
\description{
  Hotelling's T2 test for the one and two sample case.
}

\usage{ HotellingsT2(X, \dots)

\method{HotellingsT2}{default}(X, Y = NULL, mu = NULL, test = "f",
             na.action = na.fail, \dots)

\method{HotellingsT2}{formula}(formula, na.action = na.fail, \dots) }

\arguments{
  \item{X}{a numeric data frame or matrix. }
  \item{Y}{an optional numeric data frame or matrix for the two sample test. If NULL a one sample test is performed.}
  \item{mu}{a vector indicating the hypothesized value of the mean (or difference
          in means if a two sample test is performed). NULL represents origin or no difference between the groups. }
  \item{test}{if 'f', the decision is based on the F-distribution, if 'chi' a chi-squared approximation is used. }
  \item{formula}{a formula of the form \code{X ~ g} where \code{X}
    is a numeric matrix giving the data values and \code{g} a factor
    with two levels giving the corresponding groups.}
  \item{na.action}{a function which indicates what should happen when the data
          contain 'NA's. Default is to fail. }
  \item{\dots}{further arguments to be passed to or from methods.}

}

\details{

  The classical test for testing the location of a multivariate population or for testing the mean
  difference for two multivariate populations. When \code{test = "f"} the F-distribution is used for
  the test statistic and it is assumed that the data are normally distributed. If the chisquare
  approximation is used, the normal assumption can be relaxed to existence of second moments.
  In the two sample case both populations are assumed to have the same covariance matrix.

  The formula interface is only applicable for the 2-sample tests.

}

\value{
  A list with class 'htest' containing the following components:
  \item{statistic }{the value of the T2-statistic. (That is the scaled value of the statistic that has an 
                    F distribution or a chisquare distribution depending on the value of \code{test}).}
  \item{parameter}{the degrees of freedom for the T2-statistic.}
  \item{p.value}{the p-value for the test.}
  \item{null.value}{the specified hypothesized value of the mean or mean difference
   depending on whether it was a one-sample test or a two-sample test.}
  \item{alternative}{a character string with the value 'two.sided'.}
  \item{method}{a character string indicating what type of test was performed.}
  \item{data.name}{a character string giving the name of the data (and grouping vector).}

}

\references{ \cite{Anderson, T.W. (2003), \emph{An introduction to
multivariate analysis}, New Jersey: Wiley. }  }

\author{Klaus Nordhausen}

\examples{ 
# one sample test:

data(pulmonary)

HotellingsT2(pulmonary) 
HotellingsT2(pulmonary, mu = c(0,0,2), test = "chi")

# two sample test:

set.seed(123456)
X <- rmvnorm(20, c(0, 0, 0, 0), diag(1:4))
Y <- rmvnorm(30, c(0.5, 0.5, 0.5, 0.5), diag(1:4))
Z <- rbind(X, Y)
g <- factor(rep(c(1,2),c(20,30)))

HotellingsT2(X, Y)
HotellingsT2(Z ~ g, mu = rep(-0.5,4))

rm(.Random.seed) }

\keyword{htest } 
\keyword{ multivariate }
