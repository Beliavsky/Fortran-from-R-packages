\name{star}
\alias{star}
\docType{data}
\title{Star Brightness}
\description{Brightness (magnitude)
of a particular star at midnight on 600 consecutive nights.
}
\usage{data(star)}
\source{Whittaker, E. T. and Robinson, G., (1924). The Calculus of Observations. London:
Blackie and Son.}
\examples{
data(star)
## maybe str(star) ; plot(star) ...
data(star)
plot(star,xlab='Day',ylab='Brightness')
}
\keyword{datasets}
