\name{ad.test}
\alias{ad.test}
\title{ Anderson-Darling GoF test }
\description{
    Implementation of the Anderson-Darling goodness of fit test. 
}
\usage{
ad.test(x, distr.fun, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a random sample from a possibly unknown continuous distribution }
  \item{distr.fun}{ a named CDF, such as \code{pnorm}, \code{punif}, etc. }
  \item{\dots}{ extra parameters for the distribution function above, such as location and scale parameters, etc. }
}
\details{
  If the \code{distr.fun} is provided, the function checks whether \code{x} is a iid sample from the distribution described by such CDF.
  Otherwise, whether they follow a uniform law.
}
\value{
    The output is an object of the class \code{htest} exactly like for the Kolmogorov-Smirnov test, \code{\link{ks.test}}. 
    The \code{statistic} and \code{p.value} fields are the most relevant ones.
}
\references{ G. and J. Marsaglia, "Evaluating the Anderson-Darling Distribution", Journal of Statistical Software, 2004 }
\author{ Carlos J. Gil Bellosta }
%\note{ }
%\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
    set.seed( 123 )
    x <- runif( 100 )

    ad.test( x )$p.value

    ad.test( x, pnorm, 0, 1 )$p.value

    replicate( ad.test( rnorm( 100 ), pnorm )$p.value, 100 )
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
\keyword{ htest }
