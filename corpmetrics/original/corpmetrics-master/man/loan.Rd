\name{loan}
\alias{loan}
\title{
Loan Payments & Debt Management
}

\description{
Calculate the installments, interest, principal and debt balance for each period of a loan.
}

\usage{
loan(AMOUNT,RATE,PER)
}

\arguments{
\item{AMOUNT}{
Amount of the loan (Numeric Variable).
}
\item{RATE}{
Loan's periodic interest rate (Numeric Variable).
}
\item{PER}{
Periods to maturity (Numeric Variable).
}
}

\value{
A list with 2 data.frames: (1) Summary, presenting the installment and the repayment amount. (2) Amortization Table for all periods of the loan until maturity.
}

\references{
Ehrhardt, M. C. and Brigham, E. F. (2011). Corporate finance: A focused approach. South-Western Cengage Learning. p. 156-160 - ISBN: 1439078084
}

\author{
Pavlos Pantatosakis.

R implementation and documentation: Pavlos Pantatosakis \email{pantatosakisp@yahoo.com}.
}


\seealso{
\code{\link{idm}
}
}

\examples{
##
# Example usage
example <- loan(
  AMOUNT = 100000, # 100.000 currency units loan
  RATE = 0.05, # 5 % periodic interest
  PER = 4 # 4 periods to maturity
)

print(example)

# Another example

example2 <- loan(
  AMOUNT = 1000, # 1.000 currency units loan
  RATE = 0.20, # 20 % periodic interest
  PER = 3 # 3 periods to maturity
)

print(example2$AmortizationTable) # For amortization table
print(example2$AmortizationTable[5]) # Balance for each period
}

