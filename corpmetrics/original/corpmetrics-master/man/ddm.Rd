\name{ddm}
\alias{ddm}
\title{
Dividend Discount Models (DDM)
}

\description{
Calculate the value of a common stock from discounted dividends, by employing (1) Zero Growth Model, (2) Gordon's Model, (3) Differential Growth Model.
}

\usage{
ddm(DIV,RETURN,G1,G2,PER)
}

\arguments{
\item{DIV}{
Dividend at period 0 (Numeric variable).
}
\item{RETURN}{
Required return of the investor (Numeric variable).
}
\item{G1}{
Expected growth rate (Numeric variable) - Optional (Essential for Gordon's model & the differential growth model).
}
\item{G2}{
Expected growth rate after the period of change (Numeric variable) - Optional (Essential for the differential growth model).
}
\item{PER}{
Period at which the growth rate changes (Numeric variable) - Optional (Essential for the differential growth model).
}
}

\details{
For the Zero Growth Model, fill in DIV and RETURN; for Gordon's Model, include DIV, RETURN, and G1; and for the Differential Growth Model, provide DIV, RETURN, G1, G2, and PER.
}

\value{
A data.frame presenting the model employed and the stock's value based on discounted dividents.
}

\references{
Jordan, B. D., Ross, S. A., and Westerfield, R. W. (2010). Fundamentals of corporate finance. McGraw Hill. p. 234-240 - ISBN: 9780073382395
}

\author{
Pavlos Pantatosakis.

R implementation and documentation: Pavlos Pantatosakis \email{pantatosakisp@yahoo.com}.
}


\examples{
##
# Example usage

#Company pays a dividend of 3 currency units per share
#Investors require a return of 8%

example <- ddm(
  DIV = 3, # Dividend Amount in currency units
  RETURN = 0.08 # Required Return of the investor
)

print(example)

#Company pays a dividend of 0.8 currency units per share
#Investors require a return of 10%
#The dividend is expected to grow at a constant rate of 4%

example2 <- ddm(
  DIV = 0.8, # Dividend Amount in currency units
  RETURN = 0.10, # Required Return of the investor
  G1 = 0.04 # Growth rate
)

print(example2)

#Company pays a dividend of 2 currency units per share
#Investors expect a return of 12%
#The dividend is projected to grow at 8% for the first 3 years
#Then at 4%

example3 <- ddm(
  DIV = 2, # Dividend Amount in currency units
  RETURN = 0.12, # Required Return of the investor
  G1 = 0.08, # Growth rate
  G2 = 0.04, # Growth rate after PER
  PER = 3 # Growth rate change happens in Period 3
)

print(example3)
}

