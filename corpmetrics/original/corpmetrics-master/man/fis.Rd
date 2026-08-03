\name{fis}
\alias{fis}
\title{
Price & Macaulay Duration of a Fixed Income Security
}

\description{
Calculate discounted cash flows of a bond's coupon payments and its sensitivity to interest rate change.
}

\usage{
fis(FV,CR,YTM,MAT,SEMI)
}

\arguments{
\item{FV}{
Face value of the bond (Numeric Variable).
}
\item{CR}{
Coupon rate of the bond (Numeric Variable).
}
\item{YTM}{
Yield to maturity (Numeric Variable).
}
\item{MAT}{
Maturity in years (Numeric Variable).
}
\item{SEMI}{
Select between annual or semi-annual coupon payments (default is annual).
}
}

\details{
The default option is annual coupon payments. To select semi-annual payments, set SEMI = TRUE.
}

\value{
A data.frame with the results of the financial instrument's price in currency units, its duration and modified duration in years.
}

\references{
Jordan, B. D., Ross, S. A., and Westerfield, R. W. (2010). Fundamentals of corporate finance. McGraw Hill. p. 193-202. - ISBN: 9780073382395

Berk, J. B. and DeMarzo, P. M. (2017). Corporate finance. Pearson Education. p.205-220 & p.1073-1074 -  ISBN: 1292160160
}

\author{
Pavlos Pantatosakis.

R implementation and documentation: Pavlos Pantatosakis \email{pantatosakisp@yahoo.com}.
}

\examples{
##
# Example usage

# Face value = 1,000 currency units
# Coupon rate = 8%
# Yield to maturity = 8%
# Maturity = 6 years

example <- fis(
  FV = 1000, # Bond with face value of 1.000 currency units
  CR = 0.08, # 8% Coupon rate
  YTM = 0.08, # 8% Yield to maturity
  MAT = 6 # 6 periods to maturity
)

print(example)

# Face value = 1,000 currency units
# Coupon rate = 8%
# Yield to maturity = 12%
# Maturity = 2 years
# Coupons pay semi-annually

example2 <- fis(
  FV = 1000, # Bond with face value of 1.000 currency units
  CR = 0.08, # 8% Coupon rate
  YTM = 0.12, # 8% Yield to maturity
  MAT = 2, # 6 periods to maturity
  SEMI = TRUE
)

print(example2)
}

