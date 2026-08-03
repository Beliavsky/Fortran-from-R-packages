\name{balsh}
\alias{balsh}
\title{
Balance Sheet Key Performance Indicators (KPIs)
}

\description{
Calculate balanace sheet ratios: (1) Working Capital, (2) Current Ratio, (3) Acid Test Ratio, (4) Leverage (L/A), (5) Debt-to-Equity.
}

\usage{
balsh(FA,CA,INV,FL,CL)
}

\arguments{
\item{FA}{
Fixed Assets (Numeric Variable).
}
\item{CA}{
Current Assets (Numeric Variable).
}
\item{INV}{
Inventory (Numeric Variable).
}
\item{FL}{
Fixed (Long-term) Liabilities (Numeric Variable).
}
\item{CL}{
Current (Short-term) Liabilities (Numeric Variable).
}
}

\details{
Total Assets, Total Liabilities and Total Equity are computed in function.
}

\value{
A data.frame with the 5 metrics and their respective values.
}

\references{
The Investopedia Team (2024. Useful Balance Sheet Metrics. https://www.investopedia.com/financial-edge/1012/useful-balance-sheet-metrics.aspx)
}

\author{
Pavlos Pantatosakis.

R implementation and documentation: Pavlos Pantatosakis \email{pantatosakisp@yahoo.com}.
}

\examples{
##
# Example usage
example <- balsh(
  FA = 2450000, # Fixed Assets
  CA = 770000,  # Current Assets
  INV = 450000, # Inventory
  FL = 1180000, # Fixed Liabilities
  CL = 490000   # Current Liabilities
)

print(example)
}

