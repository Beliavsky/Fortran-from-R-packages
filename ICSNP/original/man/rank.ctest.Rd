\name{rank.ctest}
\alias{rank.ctest}
\alias{rank.ctest.default}
\alias{rank.ctest.formula}
\alias{rank.ctest.ics}

\title{ One, Two and C Sample Rank Tests for Location based on Marginal Ranks }
\description{
  Performs the one, two or c sample location test based on marginal ranks. Three different score functions are available.
}

\usage{rank.ctest(X, \dots)

\method{rank.ctest}{default}(X, Y = NULL, mu = NULL, scores = "rank", 
           na.action = na.fail, \dots)

\method{rank.ctest}{formula}(formula, na.action = na.fail, \dots)

\method{rank.ctest}{ics}(X, g = NULL, index = NULL, na.action = na.fail, \dots)
}

\arguments{
  \item{X}{a numeric data frame or matrix or an ics object. }
  \item{Y}{an optional numeric data frame or matrix for the two sample test. If NULL a one sample test is performed.}
  \item{mu}{a vector indicating the hypothesized value of the mean (or difference
          in means if you are performing a two sample test). NULL represents origin or no difference between the groups.
          For more than two groups \code{mu} should be 0 or not be specified at all.}
  \item{scores}{ if 'sign', a sign test is performed, if 'rank' a signed rank test is performed or if 'normal'
        a normal score test is performed. }
  \item{formula}{a formula of the form \code{X ~ g} where \code{X}
    is a numeric matrix giving the data values and \code{g} a factor
    with at least two levels giving the corresponding groups.}
  \item{g}{a grouping factor with at least two levels.}
  \item{index}{an integer vector that gives the columns to choose the invariant coordinates form the 'ics' object. The default uses all columns.}
  \item{na.action}{  a function which indicates what should happen when the data
          contain 'NA's.  Default is to fail. }
  \item{\dots}{further arguments to be passed to or from methods.}
}
\details{
  These tests are well described in Puri and Sen (1971). The tests are based on the marginal ranks for which three score functions are available. The scores are
  also used to estimate the covariance matrices. In the multisample case it is assumed that the distribution of the different populations differs only in their location.

  The ics interface provides an invariant test based on the invariant coordinate selection. The assymptotic distribution is however still an open question
  when more than one component is used, though the chi-square approximation works well also for several components as shown in Nordhausen, Oja and Tyler (2006).}

\value{
  A list with class 'htest' containing the following components:
  \item{statistic }{the value of the T-statistic.}
  \item{parameter}{the degrees of freedom for the T-statistic.}
  \item{p.value}{the p-value for the test.}
  \item{null.value}{the specified hypothesized value of the mean or mean difference depending on whether it was a one-sample test or a two-sample test.}
  \item{alternative}{a character string with the value 'two.sided'.}
  \item{method}{a character string indicating what type of test was performed.}
  \item{data.name}{a character string giving the name of the data (and grouping vector).}
}

\references{ \cite{Puri , M.L. and Sen, P.K. (1971), \emph{Nonparametric Methods in Multivariate Analysis}, New York: Wiley.  }

            \cite{Nordhausen, K., Oja, H. and Tyler, D.E. (2006), On the Efficiency of Invariant Multivariate Sign and Rank Tests, in 
            \emph{Festschrift of Tarmo Pukkila on his 60th Birthday}, 217--231.}}
            
\author{Klaus Nordhausen}

\examples{
# one sample tests:

data(pulmonary)

rank.ctest(pulmonary, scores = "sign")
rank.ctest(pulmonary, mu = c(0,0,2))

# two sample tests:

set.seed(123456)
X <- rmvnorm(20, c(0,0,0,0), diag(1:4))
Y <- rmvnorm(30, c(0.5,0.5,0.5,0.5), diag(1:4))
Z <- rbind(X,Y)
g <- factor(rep(c(1,2), c(20,30)))

rank.ctest(X, Y, scores = "normal")
rank.ctest(Z~g, scores = "sign", mu = rep(-0.5,4))

# c sample test:

W <- rmvnorm(30, c(0,0,0,0), diag(1:4))
Z2 <- rbind(X,Y,W)
g2 <- factor(rep(1:3, c(20,30,30)))

rank.ctest(Z2~g2, scores = "normal")

# in an invariant coordinate system

rank.ctest(ics(Z2,covOrigin, cov4, S2args=list(location =
           "Origin")), index = c(1,4), scores = "sign")

rank.ctest(ics(Z), g, index = 4)

rank.ctest(ics(Z2), g2, scores = "normal",index = 4)

rm(.Random.seed)
         }


\keyword{htest }
\keyword{ multivariate }
\keyword{nonparametric}
