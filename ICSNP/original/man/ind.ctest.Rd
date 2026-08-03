\name{ind.ctest}
\alias{ind.ctest}

\title{Test of Independece based on Marginal Ranks}
\description{
  Performs the test that a group of variables is independent of an other based on marginal ranks. Three different score functions are available.
}
\usage{
ind.ctest(X, index1, index2 = NULL, scores = "rank", 
          na.action = na.fail)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a data frame or matrix. }
  \item{index1}{integer vector that selects the columns of \code{X} that form group one. Only numeric columns can be selected. }
  \item{index2}{integer vector that selects the columns of \code{X} that form group two. Only numeric columns can be selected.
                If \code{NULL}, all remaining columns of X will be selected.}
  \item{scores}{if 'sign', a sign test is performed, if 'rank' a rank test is performed or if 'normal'
        a normal score test is performed.}
  \item{na.action}{a function which indicates what should happen when the data
          contain 'NA's.  Default is to fail.}
}
\details{
  The test tests if X[ , index1] is independent of X[ , index2] and is described in great detail in Puri and Sen (1971).
}
\value{
  A list with class 'htest' containing the following components:
  \item{statistic }{the value of the W-statistic.}
  \item{parameter}{the degrees of freedom for the W-statistic.}
  \item{p.value}{the p-value for the test.}
  \item{method}{a character string indicating what type of test was performed.}
  \item{data.name}{a character string giving the name of the data.}
}
\references{\cite{Puri , M.L. and Sen, P.K. (1971), \emph{Nonparametric Methods in Multivariate Analysis}, New York: Wiley.   }}

\author{Klaus Nordhausen}

\examples{
A1 <- matrix(c(4, 4, 5, 4, 6, 6, 5, 6, 7), ncol = 3)
A2 <- matrix(c(0.5, -0.3, -0.3, 0.7), ncol = 2)
X <- cbind(rmvnorm(100, c(-1, 0, 1), A1), rmvnorm(100, c(0, 0), A2))
ind.ctest(X,1:3)
ind.ctest(X, c(1, 5), c(2, 3), scores = "normal")
}
\keyword{htest }
\keyword{ multivariate }
\keyword{nonparametric}
