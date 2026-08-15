\name{Bell Touchard}
\alias{dbellt}
\alias{pbellt}
\alias{qbellt}
\alias{rbellt}
\alias{mle_bt}
\title{PMF, CDF, QF, random generation and parameters estimation based on the Bell Touchard distribution}
\description{Evaluates the PMF, CDF, QF, random generation and MLE based on the Bell Touchard distribution.
The PMF of the Bell Touchard distribution is as follows:
\deqn{
f\left(X=x\mid\lambda,\theta\right)=\exp\left\{ \theta\left[1-e^{\lambda}\right]\right\} \frac{\lambda^{x}\,T_{x}\left(\theta\right)}{x!};\qquad x=0,1,2,\dots,
}
where  \eqn{\lambda>0} and \eqn{\theta>0} are the two parameters and \eqn{T_{x}} are the Touchard polynomials, it is given by
\deqn{T_{n}=\frac{1}{e}\sum_{k=0}^{\infty}\frac{k^{n}}{k!}.}
It is important to note that when the parameter \eqn{\theta=1}, the Bell Touchard distribution reduces to Bell distribution.
}
\usage{
dbellt (x, lambda, theta, log = FALSE)
pbellt (q, lambda, theta, lower.tail = TRUE, log.p = FALSE)
qbellt (p, lambda, theta, lower.tail = TRUE, log.p = FALSE)
rbellt (n, lambda, theta)
mle_bt (x, lambda, theta)
}
\arguments{
\item{x}{
A vector of (non-negative integer) discrete values.
}
\item{lambda}{
A vector of (non-negative integer) values, \eqn{\lambda>0}.
}
\item{theta}{
A vector of (non-negative integer) values, \eqn{\theta>0}.
}
\item{q}{
A probablity or a vector of probablities.
}
\item{p}{
A probablity or a vector of probablities.
}
\item{n}{
A randomly generated values.
}
\item{lower.tail}{
logical; if TRUE (default), probabilities will be P [X <= x], otherwise, P [X> x].
}
\item{log}{
logical; if TRUE, probabilities p are given as log(p).
}
\item{log.p}{
logical; if TRUE, probabilities p are given as log(p).
}
}





\details{
Recently Castellares et al. (2020) proposed a two parameter discrete Bell Touchard
distribution that overcomes the issue of over-dispersion (variance larger than mean). Often, we experience the over-dispersed data in practice, where the Poisson model may not be suitable
because its variance is restricted to be equal to mean. The PMF of the Bell Touchard
distribution is simple and tractible in order to find probablities as well as several properties of the distribution.
}

\value{
dbellt gives the (log) probability function.
pbellt gives the (log) distribution function.
qbellt gives the quantile function.
rbellt generates random values.
mle_bt gives the maximum likelihood estimates with standard error of
the estimates and model selection measure, the Akaike information criterion (AIC).
}

\references{
Castellares, F., Lemonte, A. J., & Moreno–Arenas, G. (2020). On the two-parameter Bell–Touchard discrete distribution. Communications in Statistics-Theory and Methods, 49(19), 4834-4852.

Castellares, F., Ferrari, S. L., & Lemonte, A. J. (2018). On the Bell distribution and its associated regression model for count data. Applied Mathematical Modelling, 56, 172-185.
}

\author{
Muhammad Imran and M.H. Tahir.

R implementation and documentation: Muhammad Imran \email{imranshakoor84@yahoo.com} and M.H. Tahir <mht@iub.edu.pk>.
}


%\note{
%%  ~~further notes~~
%}
\seealso{
\code{\link{mle_borel}, \link{mle_poisson}, \link{bell_mle}
}
}

\examples{
x<-2
dbellt (x, 2, 2)
pbellt (0.2, 2, 2)
pbellt (0.2, 2, 2)
rbellt (10, 2, 1)
x <- data_sbirth
mle_bt (x, 0.12, 0.2)

# Plot of PMF of the Bell Touchard distribution for selected parameteric values.
# The Bell Touchard distribution reduces to the Bell distribution for fixing theta=1.
x <- 0:50
lambda <- 2
theta <- 1.5

plot(dbellt(x, lambda, theta), type = "h", col="red",lwd = 3,
     main = "Bell Touchard distribution",
     ylab = "P(X = x)", xlab = "Number of events")

legend("topright",c(expression(lambda==2~~theta==1.5)),lty=1, col="red", lwd=2,cex=1.0)

}



