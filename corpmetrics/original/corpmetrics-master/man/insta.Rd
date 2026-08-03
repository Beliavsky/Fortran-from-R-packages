\name{insta}
\alias{insta}
\title{
P & L Key Performance Indicators (KPIs)
}

\description{
Calculate income statement ratios: (1) Gross Profit Margin, (2) Net Profit Margin, (3) Earnings Per Share (EPS), (4) Price to Earnings (P/E) Ratio.
}

\usage{
insta(REV,COS,NET,PREF,SHARES,PPS)
}

\arguments{
\item{REV}{
Revenue (Numeric Variable).
}
\item{COS}{
Cost of Sales (Numeric Variable).
}
\item{NET}{
Net Income (Numeric Variable).
}
\item{PREF}{
Amount of Preferred Stock Dividend (Numeric Variable) - Optional (Essential for the calculation of EPS & P/E Ratio).
}
\item{SHARES}{
Number of Shares (Numeric Variable) - Optional (Essential for the calculation of EPS & P/E Ratio).
}
\item{PPS}{
Price per Share (Numeric Variable) - Optional (Essential for the calculation of P/E Ratio).
}
}

\value{
A data.frame with the 4 ratios and their respective values.
}

\references{
Hayes A. (2024). Profitability Ratios: What They Are, Common Types, and How Businesses Use Them. https://www.investopedia.com/terms/p/profitabilityratios.asp

Picardo E. (2024). P/E Ratio vs. EPS vs. Earnings Yield: What's the Difference? https://www.investopedia.com/articles/investing/120513/comparing-pe-eps-and-earnings-yield.asp
}

\author{
Pavlos Pantatosakis.

R implementation and documentation: Pavlos Pantatosakis \email{pantatosakisp@yahoo.com}.
}


\examples{
##
# Example usage
example <- insta(
  REV = 25000000, # Revenue
  COS = 19850000, # Cost of Sales
  NET = 1000000, # Net Income
  PREF = 100000, # Preferred Stock Dividend
  SHARES = 100000, # Number of Shares
  PPS = 120 # Price per Share
)

print(example)
}



