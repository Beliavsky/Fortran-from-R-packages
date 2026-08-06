\name{QCSIS}
\alias{QCSIS}
\alias{QCSIS}
\docType{package}
\title{ Quantile Correlation-Sure Independence Screening (QC-SIS)
}
\description{
The function implemrnts the  quantile correlation-sure independence  screening (QC-SIS).
}
\usage{
QCSIS(x, y, tau, d)
}

\arguments{
\item{x}{The design matrix, of dimensions n * p, without an intercept.}
\item{y}{The response vector of dimension n * 1.}
\item{tau}{The quantile(s) to be estimated. By default, tau=1:(n-1)/n. }
\item{d}{The tuning parameter used to covarites had significant effect on the response variable, such as [n/log(n)],or n-1}
}

\value{
\item{w}{The estimate of w.}
\item{M}{The subscript of x recuited by QC-SIS.}
}

\author{
\packageAuthor{QCSIS}
}

\references{
Xuejun Ma  and Jingxiao Zhang. Robust model-free feature screening via quantile correlation. Journal of Multivariate Analysis. Online, 2015.\cr
}


\examples{
n <- 20
p <- 200
r <- 0.05
x <- matrix(rnorm(n * p), n, p)
e <- rnorm(n, 0, 1)
inde <- sample(n, r * n)
x[inde, 1] <- 2 * sqrt(rchisq(r * n, df = p))
y <- 5 * x[, 1] + 5 * x[, 2] + 5 * x[, 3] + e
d <- 19
fit.QCSIS <- QCSIS(x = x, y = y, d = d)
fit.QCSIS$M
}


