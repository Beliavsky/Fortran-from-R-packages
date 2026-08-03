\name{capm}
\alias{capm}
\title{
Capital Asset Pricing Model (CAPM)
}

\description{
Calculate if a stock is fairly valued based on Capital Asset Pricing Model (CAPM) returns.
}

\usage{
capm(Rf,Ri,Rm)
}

\arguments{
\item{Rf}{
A numerical vector with the risk-free rates.
}
\item{Ri}{
A numerical vector with the yields of the asset.
}
\item{Rm}{
A numerical vector with the yields of the market.
}
}

\details{
Use the same length for the three vectors.
}

\value{
A data.frame including the required return based on Capital Asset Pricing Model (CAPM), the expected (mean) return of the asset and the asset's beta. Lastly, the stock's valuation is assessed: the stock is considered overvalued if the required return exceeds the expected return, and undervalued if the expected return is higher than the required return.
}

\references{
Will Kenton (2024). Capital Asset Pricing Model (CAPM): Definition, Formula, and Assumptions. https://www.investopedia.com/terms/c/capm.asp
}

\author{
Pavlos Pantatosakis.

R implementation and documentation: Pavlos Pantatosakis \email{pantatosakisp@yahoo.com}.
}

\examples{
##
# Example usage

Rf <- rnorm(250, 0.03, 0) # Constant 3% risk free rate
Ri <- rnorm(250, 0.13, 0.10) # Asset under study
Rm <- rnorm(250, 0.09, 0.04) # Market's (can be an index) returns

example <- capm(Rf,Ri,Rm)
print(example)

#You can add a data.frame with real data
#Choose a vector with the risk free rate (Rf)
#Choose a vector with the stock returns (Ri)
#Choose a vector with the market returns (Rm)
}

