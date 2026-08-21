\name{rdata}
\alias{rdata}
\docType{data}
\title{Survival Data}
\description{
  Survival data.
}
\usage{data(rdata)}
\format{
  A data frame with 1040 observations on the following 6 variables:
 \describe{
    \item{time}{survival time (in days).}
    \item{cens}{censoring indicator (0=censoring, 1=death).}
    \item{age}{age (in years).}
    \item{sex}{sex (1=male, 2=female).}
    \item{year}{date of diagnosis (in date format).}
    \item{agegr}{age group.}
 }
}

\references{
Pohar M., Stare J. (2006) "Relative survival analysis in R." Computer Methods and Programs in Biomedicine, \bold{81}: 272-278.
}


\keyword{datasets}
