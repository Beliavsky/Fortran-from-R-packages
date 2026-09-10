\name{Oakley}
\alias{Oakley}
\alias{Oakley1}
\alias{Oakley2}
\title{RFC 2409 Oakley Groups - Parameters for Diffie-Hellman Key Exchange}
\description{
  RFC 2409 standardizes global unique prime numbers and generators for
  the purpose of secure asymmetric key exchange on the Internet.
}
\usage{
data(Oakley1)
data(Oakley2)
}
\value{
  Oakley1 returns an object of class \code{\link{bigz}} for a 768 bit Diffie-Hellman group.  The generator is stored as value with the respective prime number as modulus attribute.

  Oakley2 returns an object of class \code{\link{bigz}} for a 1024 bit Diffie-Hellman group.  The generator is stored as value with the respective prime number as modulus attribute.
}
\references{The Internet Key Exchange (RFC 2409), Nov. 1998}

\examples{
packageDescription("gmp") # {possibly useful for debugging}

data(Oakley1)
(M1 <- modulus(Oakley1))
isprime(M1)# '1' : "probably prime"
sizeinbase(M1)#  232 digits (was 309 in older version)
}
\keyword{data}

