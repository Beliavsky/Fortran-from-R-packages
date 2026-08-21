\name{Goodness of fit test for the bivariate Poisson distribution}
\alias{bp.gof}
\alias{bp.gof2}
\title{
Goodness of fit test for the bivariate Poisson distribution
}

\description{
Goodness of fit test for the bivariate Poisson distribution.
}

\usage{
bp.gof(x1, x2 = NULL, R = 999)
bp.gof2(x1, x2 = NULL, R = 999)
}

\arguments{
\item{x1}{
Either a numerical vector with the values of the first variable or a matrix with 2 columns containing both variables. In the latter case, x2 must be NULL.
}
\item{x2}{
A numerical vector with the values of the second. If x1 is a matrix with 2 columns containing both variables, x2 must be NULL.}
\item{R}{
The number of Monte Carlo replicates to use.
}
}

\details{
Kocherlakota and Kocherlakota (1992) mention the following a goodness of fit test for the
bivariate Poisson distribution, the index of dispersion test. They mention that Loukas and Kemp
(1986) developed this test as an extension of the univariate dispersion test. They test for
departures from the bivariate Poisson againsta alternatives which involve an increase in the
generalised variance, the determinant of the covariance matrix of the two variables.

Rayner, Thas and Best (2009) mentions a revised version of this test whose test statistic is now given by
\deqn{
I_{B^*}=\frac{n}{1-r^2}\left(\frac{S_1^2}{\bar{x}_1}-2r^2\sqrt{\frac{S_1^2}{\bar{x}_1}\frac{S_2^2}{\bar{x}_2}}+\frac{S_2^2}{\bar{x}_2}\right),
}

where \eqn{n} is the sample size, \eqn{r} is the sample Pearson correlation coefficient,
\eqn{S_1^2} and \eqn{S_2^2} are the two sample variances and \eqn{\bar{x}_1} and \eqn{\bar{x}_2}
are the two sample means. Under the null hypothesis the \eqn{I_{B^*}} follows asymptotically a
\eqn{\chi^2} with \eqn{2n-3} degrees of freedom. However, I did some simulations and I saw that
it does not perform very well in terms of the type I error. If you see the simulations in their book
(page 132) you will see this. For this reason, the function calculates the p-value of the
\eqn{I_{B^*}} using Monte Carlo (or parametric bootstrap).

The second function, bp.gof2(), is a vectorised version of the first, and much faster. I put both of
them here to show how one can vectorize a function and make it faster.
}

\value{
A list including:
\item{runtime}{
The duration of the algorithm.
}
\item{tab}{
The contingency table of the two variables.
}
\item{pvalue}{
The Monte-Carlo based estimated p-value.
}
}

\references{
Kocherlakota S. and Kocherlakota K. (1992). Bivariate discrete distributions. CRC Press.

Loukas S. and Kemp C. (1986). The index of dispersion test for the bivariate Poisson distribution.
Biometrics, 42(4): 941--948.

Rayner J. C., Thas O. and Best D. J. (2009). Smooth Tests of Goodness of Fit: Using R. John
Wiley & Sons.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}
\seealso{
\code{\link{bp.mle} }
}

\examples{
x <- rbp( 300, c(3, 5, 2) )
bp.gof(x)
bp.gof2(x)
}

