\name{Random}
\alias{urand.bigz}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Generate a random number}
\description{
  Generate a uniformly distributed random number in the range 0 to
  \eqn{2^{size} -1}{2^size - 1}, inclusive.

%  If seed is precise, a new seed is initialised.
}
\usage{
urand.bigz(nb=1,size=200, seed = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nb}{Integer: number of random numbers to be generated (size of
    vector returned)}
  \item{size}{Integer: number will be generated in the range 0 to
  \eqn{2^{size} -1}{2^size - 1}}
  \item{seed}{Bigz: random seed initialisation}
}

\value{
  A biginteger of class bigz.
}
\references{
  \samp{mpz\_urandomb} from the GMP Library, see \url{https://gmplib.org}
}
\author{Antoine Lucas}

\examples{
# Integers are differents
urand.bigz()
urand.bigz()
urand.bigz()

# Integers are the same
urand.bigz(seed="234234234324323")
urand.bigz(seed="234234234324323")

# Vector
urand.bigz(nb=50,size=30)
}
\keyword{arith }% at least one, from doc/KEYWORDS

