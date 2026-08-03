\name{idm}
\alias{idm}
\title{
Net Present Value (NPV) and Internal Rate of Return (IRR) of an investment
}

\description{
Basic investment decision methods: (1) Net Present Value (NPV), (2) Internal Rate of Return (IRR).
}

\usage{
idm(CFS,COST)
}

\arguments{
\item{CFS}{
A numerical vector with the investment's expected cash flows.
}
\item{COST}{
A numerical vector with the investment's expected cost of capital (interest rates).
}
}

\details{
The first cash flow must be the initial cost of investment, thus a negative value. Cash flows and interest rates must have the same length.
}

\value{
A data.frame with the investment decision methods and their respective values.
}

\references{
Berk, J. B. and DeMarzo, P. M. (2017). Corporate finance. Pearson Education. p.100-103 & p.248-251 - ISBN: 1292160160

Ehrhardt, M. C. and Brigham, E. F. (2010). Corporate finance: A focused approach. South-Western Cengage Learning. p. 383-389 - ISBN: 1439078084

}

\author{
Pavlos Pantatosakis.

R implementation and documentation: Pavlos Pantatosakis \email{pantatosakisp@yahoo.com}.
}

\seealso{
\code{\link{uniroot}
}
}

\examples{
##
# Example usage
# Initial Investment = 100 currency units
# Expected to bring 120 currency units the next period
# The cost of capital is 10%

example <- idm(
  CFS = c(-100,120),
  COST = c(0,0.1)
)

print(example)
}

