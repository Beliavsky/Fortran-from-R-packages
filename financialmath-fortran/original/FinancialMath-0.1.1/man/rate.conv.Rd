\arguments{
  \item{rate}{current rate}
  \item{conv}{how many times per year the current rate is convertible}
  \item{type}{current rate as one of "interest", "discount" or "force"}
  \item{nom}{desired number of times the calculated rates will be convertible}
}
\name{rate.conv}
\alias{rate.conv}
\title{Interest, Discount, and Force of Interest Converter}
\description{Converts given rate to desired nominal interest, discount, and force of interest rates.}
\usage{rate.conv(rate, conv=1, type="interest", nom=1)}
\author{Kameron Penn and Jack Schmidt}
\keyword{interest}
\value{A matrix of the interest, discount, and force of interest conversions for effective, given and desired conversion rates.

The row named 'eff' is used for the effective rates, and the nominal rates are in a row named 'nom(\emph{x})' where the rate is convertible \emph{x}  times per year.}
\examples{rate.conv(rate=.05,conv=2,nom=1)

rate.conv(rate=.05,conv=2,nom=4,type="discount")

rate.conv(rate=.05,conv=2,nom=4,type="force")}
\details{
\eqn{1+i=(1+\frac{i^{(n)}}{n})^n=(1-d)^{-1}=(1-\frac{d^{(m)}}{m})^{-m}=e^\delta}
}
