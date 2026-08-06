\name{CQCSIS}
\alias{CQCSIS}
\alias{CQCSIS}
\docType{package}
\title{ Compsote Quantile Correlation-Sure Independence Screening (CQC-SIS)
}
\description{
The function implemrnts the  composite quantile correlation-sure independence screening (CQC-SIS).
}
\usage{
CQCSIS(x, y, d)
}

\arguments{
\item{x}{The design matrix, of dimensions n * p, without an intercept.}
\item{y}{The response vector of dimension n * 1.}
\item{d}{The tuning parameter used to covarites had significant effect on the response variable, such as [n/log(n)], or n-1.}
}

\value{
\item{w}{The estimate of w.}
\item{M}{The subscript of x recuited by CQC-SIS.}
}

\author{
\packageAuthor{QCSIS}
}

\references{
 Xuejun Ma et al.. Robust feature screening via composite quantile correlation learning.   In submission. \cr
}


\examples{
n <- 20
p <- 200
x <- matrix(rnorm(n * p), n, p)
e <-  rnorm(n, 0, 1)
beta1 <- 3 - runif(1)
beta2 <- 3 - runif(1)
beta3 <- 3 - runif(1)
y <- beta1 * x[, 1] + beta2 * x[, 2] + beta3 * x[, 3] + e
d <- 19
fit.CQCSIS <- CQCSIS(x = x, y = y, d = d)
fit.CQCSIS$M


}
