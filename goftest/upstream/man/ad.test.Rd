\name{ad.test}
\alias{ad.test}
\title{
  Anderson-Darling Test of Goodness-of-Fit
}
\description{
  Performs the Anderson-Darling test
  of goodness-of-fit to a specified continuous univariate
  probability distribution.
}
\usage{
ad.test(x, null = "punif", ..., estimated=FALSE, nullname)
}
\arguments{
  \item{x}{
    Numeric vector of data values.
  }
  \item{null}{
    A function, or a character string giving the name of a function,
    to compute the cumulative distribution function for the
    null distribution.
  }
  \item{\dots}{
    Additional arguments for the cumulative distribution function.
  }
  \item{estimated}{
    Logical value indicating whether the parameters of the distribution
    were estimated using the data \code{x} (composite null hypothesis),
    or were fixed in advance (simple null hypothesis, the default).
  }
  \item{nullname}{
    Optional character string describing the null distribution.
    The default is \code{"uniform distribution"}.
  }
}
\details{
  This command performs the Anderson-Darling test
  of goodness-of-fit to the distribution specified by the argument
  \code{null}. It is assumed that the values in \code{x} are
  independent and identically distributed random values, with some
  cumulative distribution function \eqn{F}.
  The null hypothesis is that \eqn{F} is the function
  specified by the argument \code{null}, while the alternative
  hypothesis is that \eqn{F} is some other function.

  By default, the test assumes that all the parameters of the null
  distribution are known in advance (a \emph{simple} null hypothesis).
  This test does not account for the effect of estimating the parameters.

  If the parameters of the distribution were estimated (that is,
  if they were calculated from the same data \code{x}),
  then this should be indicated by setting the argument \code{estimated=TRUE}.
  The test will then use the method of Braun (1980)
  to adjust for the effect of parameter estimation.

  Note that Braun's method involves randomly
  dividing the data into two equally-sized subsets, so the \eqn{p}-value
  is not exactly the same if the test is repeated.
  This technique is expected to work well when the number of
  observations in \code{x} is large.
}
\value{
  An object of class \code{"htest"} representing the result of
  the hypothesis test.
}
\references{
  Anderson, T.W. and Darling, D.A. (1952)
  Asymptotic theory of certain 'goodness-of-fit' criteria based
  on stochastic processes.
  \emph{Annals of Mathematical Statistics} \bold{23}, 193--212.

  Anderson, T.W. and Darling, D.A. (1954)
  A test of goodness of fit.
  \emph{Journal of the American Statistical Association} \bold{49}, 765--769.

  Braun, H. (1980)
  A simple method for testing goodness-of-fit in the presence of
  nuisance parameters.
  \emph{Journal of the Royal Statistical Society} \bold{42}, 53--63.
  
  Marsaglia, G. and Marsaglia, J. (2004)
  Evaluating the Anderson-Darling Distribution.
  \emph{Journal of Statistical Software} \bold{9} (2), 1--5.
  February 2004. 
  \doi{10.18637/jss.v009.i02}
}
\author{
  Original C code by George Marsaglia and John Marsaglia. 
  \R interface by Adrian Baddeley.
}
\seealso{
  \code{\link{pAD}} for the null distribution of the test statistic.
}
\examples{
x <- rnorm(10, mean=2, sd=1)
ad.test(x, "pnorm", mean=2, sd=1)
ad.test(x, "pnorm", mean=mean(x), sd=sd(x), estimated=TRUE)
}
\keyword{htest}
