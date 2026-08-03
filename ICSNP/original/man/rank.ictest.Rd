\name{rank.ictest}
\alias{rank.ictest}
\alias{rank.ictest.default}
\alias{rank.ictest.ics}

\title{ One Sample Location Test based on Marginal Ranks in the Independent Component Model }
\description{
  marginal rank test for the location problem in the one sample case when the margins are assumed independent.
}
\usage{rank.ictest(X, \dots)

\method{rank.ictest}{default}(X, mu = NULL, scores = "rank", method = "approximation",
            n.simu = 1000, na.action = na.fail, \dots)

\method{rank.ictest}{ics}(X, index = NULL, na.action = na.fail, \dots)}

\arguments{
  \item{X}{ a numeric data frame or matrix or an ics object.}
  \item{mu}{ a vector indicating the hypothesized value of the location. NULL represents the origin.}
  \item{scores}{ options are 'rank' for the signed rank test, 'sign' for the sign test and 'normal' for the normal score test. }
  \item{method}{defines the method used for the computation of the p-value. The possibilites are
        "approximation" (default), "simulation" or "permutation". Details below.}
  \item{n.simu}{if '\code{method}=simulation' or '\code{method}=permutation' this specifies the number of replications used in the
                simulation or permutation procedure. }
  \item{index}{an integer vector that gives the columns to choose from invariant coordinates form the 'ics' object. The default uses all columns.}
  \item{na.action}{a function which indicates what should happen when the data
          contain 'NA's.  Default is to fail. }
   \item{\dots}{further arguments to be passed to or from methods.}
}
\details{
  The test is normally used to test for location in the symmetric independent component model.
  
  By default the limiting distribution is used to compute the p-values. However for moderate sample sizes (N=50) was observed in 
  Nordhausen et al. (2009) that the normal score test can be sometimes slightly biased. Therefore the argument \code{method}
  can be used to get p-values based on simulations from a multivariate normal under the null or by permuting the signs of the centered
  observations.  
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
\references{ \cite{Nordhausen, K., Oja, H. and Paindaveine, D. (2009), Signed-rank tests for location in the symmetric independent component model, \emph{Journal of Multivariate Analysis}, \bold{100},  821--834. }  }
\author{Klaus Nordhausen}

\examples{
set.seed(555)
X <- cbind(rt(30,8), rnorm(30,0.5), runif(30,-3,3))
mix.matrix <- matrix(c(3,2,1,2,4,-0.5,1,-0.5,2), ncol=3)
X.mixed <- X \%*\% t(mix.matrix)
ica.X <- ics(X, covOrigin, cov4, S2args = list(location = "Origin"))
rank.ictest(ica.X)
rank.ictest(ica.X, scores = "normal", method = "simu")
rank.ictest(ics.components(ica.X), scores = "normal", method = "perm")
rm(.Random.seed)
         }


\keyword{htest}
\keyword{ multivariate}
\keyword{nonparametric}
