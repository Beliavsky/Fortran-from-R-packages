\name{Black_Scholes}
\alias{Black_Scholes}
\alias{Black_Scholes_Greeks}
\title{Black--Scholes formula and the Greeks}
\description{
  Compute the Black--Scholes formula and the Greeks.
}
\usage{
Black_Scholes(t, S, r, sigma, K, T, type = c("call", "put"))
Black_Scholes_Greeks(t, S, r, sigma, K, T, type = c("call", "put"))
}
\arguments{
  \item{t}{initial or current time \eqn{t} (in years).}
  \item{S}{stock price at time \eqn{t}.}
  \item{r}{risk-free annual interest rate.}
  \item{sigma}{annual volatility (standard deviation).}
  \item{K}{strike.}
  \item{T}{maturity (in years).}
  \item{type}{\code{\link{character}} string indicating whether
    a call (the default) or a put option is considered.}
}
\value{
  \code{Black_Scholes()} returns the value of a European-style call or put
  option (depending on the chosen \code{type}) on a non-dividend paying stock.

  \code{Black_Scholes_Greeks()} returns the first-order derivatives
  delta, theta, rho, vega and the second-order derivatives gamma, vanna
  and vomma (depending on the chosen \code{type}) in this order.
}
\details{
Note again that \code{t} is time in years. In the context of McNeil et
  al. (2015, Chapter 9), this is \eqn{\tau_t = t/250}{tau_t = t/250}.
}
\author{Marius Hofert}
\references{
   McNeil, A. J., Frey, R., and Embrechts, P. (2015).
   \emph{Quantitative Risk Management: Concepts, Techniques, Tools}.
   Princeton University Press.
}
%% \examples{
%% }
\keyword{models}