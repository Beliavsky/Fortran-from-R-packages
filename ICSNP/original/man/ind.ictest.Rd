\name{ind.ictest}
\alias{ind.ictest}
\title{Test of Independence based on Marginal Ranks in a Symmetric IC Model}
\description{
  Performs the test that a group of variables is independent of an other based on marginal ranks. It is assumed that the 
  data follows a symmetric IC model. Three different score functions are available.
}

\usage{
ind.ictest(X, index1, index2 = NULL, scores = "rank", 
           method = "approximation", n.simu = 1000, 
           ..., na.action = na.fail)
}

\arguments{
  \item{X}{a data frame or matrix.}
  \item{index1}{integer vector that selects the columns of \code{X} that form group one. Only numeric columns can be selected. }
  \item{index2}{integer vector that selects the columns of \code{X} that form group two. Only numeric columns can be selected.
                If \code{NULL}, all remaining columns of X will be selected.}
  \item{scores}{if 'sign', a sign test is performed, if 'rank' a signed rank test is performed or if 'normal'
        a normal score test is performed.}
  \item{method}{ defines the method used for the computation of the p-value. The possobilites are
        "approximation" (default), "simulation" or "permutation". Details below.}
  \item{n.simu}{if '\code{method} = "simulation"' or '\code{method} = "permutation"' this specifies the number of replications used in the
                simulation or permutation procedure. }
  \item{\dots}{further arguments to be passed to the function \code{\link[ICS]{ics}}}
  \item{na.action}{a function which indicates what should happen when the data
          contain 'NA's.  Default is to fail.}
}
\details{
  Assumed is here that \code{X[ , index1]} comes from a symmetric independent component model which in turn is independent from \code{X[ , index2]} which has also
  an underlying symmetric independent component model. This function recovers the independent components using the function \code{\link[ICS]{ics}}, centers them by a marginal
  loaction estimate based on the same scores that will be used in the actual test. The test is described in Oja, Paindaveine and Taskinen (2009). 
  The asymptotic chi-square distibution is however even for large sample sizes inadequat and therefore p-values can be simulated by resampling the test statistic under the null
  hypothesis or by permuting the rows of the independent components of \code{X[ , index2]}. Both alternatives are also described in Oja, Paindaveine and Taskinen (2009).
}
\value{
  A list with class 'htest' containing the following components:
  \item{statistic }{the value of the Q-statistic.}
  \item{parameter}{the degrees of freedom for the Q-statistic or the number of replications depending on the chosen method.}
  \item{p.value}{the p-value for the test.}
  \item{method}{a character string indicating what type of test was performed.}
  \item{data.name}{a character string giving the name of the data.}
}
\references{ \cite{Oja, H. and Paindaveine, D. and Taskinen, S. (2016), Affine-invariant rank tests for multivariate independence in independent component models, \emph{Electronic Journal of Statistics}, \bold{10},  2372--2419. }  }

\author{Klaus Nordhausen}
\examples{
Z1<-cbind(rt(500,5),rnorm(500),runif(500))
Z2<-cbind(rt(500,8),rbeta(500,2,2))
A1 <- matrix(c(4, 4, 5, 4, 6, 6, 5, 6, 7), ncol = 3)
A2 <- matrix(c(0.5, -0.3, -0.3, 0.7), ncol = 2)

X <- cbind(Z1 \%*\% t(A1), Z2 \%*\% t(A2))

ind.ictest(X,1:3)
ind.ictest(X,1:3,method="simu")

ind.ictest(X,1:2,3:5,method="perm", S1=tyler.shape,S2=cov)

}

\keyword{htest}
\keyword{ multivariate}
\keyword{nonparametric}
