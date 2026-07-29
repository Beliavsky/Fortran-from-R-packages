\name{BKg3}
\alias{BKg3}
\title{Bronshtein And Kurelenkova (2009)'s Third Risk Measure}
\description{Computes the third risk measure due to Bronshtein and Kurelenkova (2009)}
\usage{
BKg3(spec, alpha, a, b, beta, ...)
}
\arguments{
  \item{spec}{a character string specifying the distribution (for example, "norm" corresponds to the standard normal)}
  \item{alpha}{a real valued parameter taking values in (0, 1), see Chan and Nadarajah for details}
  \item{a}{the lower end point of the distribution specified by \code{spec}}
  \item{b}{the upper end point of the distribution specified by \code{spec}}
  \item{beta}{a non-negative real valued parameter, see Chan and Nadarajah for details}
  \item{...}{other parameters}}
\value{An object of the same length as \code{alpha}, giving Bronshtein and Kurelenkova (2009)'s third risk measure of the distribution specified by \code{spec}}
\references{S. Chan and S. Nadarajah, Risk: An R package for risk measures, submitted

E. Bronshtein and J. Kurelenkova, Complex risk measures in portfolio optimization,
Ufa State Aviation Technical University, Russia, 2009}
\author{Stephen Chan, Saralees Nadarajah}
\examples{BKg3("norm", 0.9, -Inf, Inf, 1)}
\keyword{Bronshtein and Kurelenkova (2009)'s third risk measure}
