\name{cqc}
\alias{cqc}
\alias{cqc}
\docType{package}
\title{Composite Quantile Correlation
}
\description{
cqc is used to compute  the composite quantile correlation.
}
\usage{
cqc(x, y)
}

\arguments{
\item{x}{The covariate variable.}
\item{y}{The response variable.}
}

\value{
\item{cqc}{The value of composite quantile correlation.}
}

\author{
\packageAuthor{QCSIS}
}
\references{
 Xuejun Ma et al.. Robust feature screening via composite quantile correlation learning.   In submission. \cr
}

\examples{
x <- rnorm(100)
y <- rnorm(100)
cqc(x = x, y = y)

}
