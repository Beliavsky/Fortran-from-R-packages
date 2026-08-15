\name{Zero inflated Bell Touchard}
\alias{dzibellt}
\alias{pzibellt}
\alias{qzibellt}
\alias{rzibellt}
\alias{mle_zibellt}
\title{MLE of the zero inflated Bell Touchard distribution}

\description{Evaluates the maximum likelihood estimate of the zero inflated Bell Touchard (ZIBELLT) distribtion. The PMF of the ZIBELLT distribution is as follows:
\deqn{
f\left(X=x\mid p_{i},\,\lambda,\theta\right)=\begin{cases}
p_{i}+\left(1-p_{i}\right)\exp\left\{ \theta\left[1-e^{\lambda}\right]\right\} , & x=0\\
\left(1-p_{i}\right)\exp\left\{ \theta\left[1-e^{\lambda}\right]\right\} \frac{\lambda^{x}\,T_{x}\left(\theta\right)}{x!}, & x=1,2\cdots,
\end{cases}
}
where \eqn{pi\in(0,1)}, \eqn{\lambda>0} and \eqn{\theta>0}  \eqn{T_{x}} are the Touchard polynomials, it is given by
\deqn{T_{n}=\frac{1}{e}\sum_{k=0}^{\infty}\frac{k^{n}}{k!}.}
It is important to note that when the parameter \eqn{\theta=1}, the ZIBELLT distribution reduces to ZIBELL distribution. On the other side, when the parameter \eqn{\theta=1} and pi=0, the ZIBELLT distribution reduces to BELL distribution. So therefore, we can evaluate the PMF, CDF, QF and random numbers of the Bell and ZIBELL distribution by using the following functions.
}

\usage{
dzibellt (x, lambda, theta, pi, log = FALSE)
pzibellt (q, lambda, theta, pi, lower.tail = TRUE, log.p = FALSE)
qzibellt (p, lambda, theta, pi, lower.tail = TRUE, log.p = FALSE)
rzibellt (n, lambda, theta, pi)
mle_zibellt (x, lambda, theta, pi)
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
\item{n}{
The number of random values generated under zero inflated Bell Touchard distribution.
}
\item{pi}{
A vector of (non-negative integer) values, \eqn{pi\in(0,1)}.
}
\item{q}{
A vector of (non-negative integer) probablities.
}
\item{p}{
A vector of (non-negative integer) probablities.
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
distribution is simple and tractible in order to find probablities as well as several properties of the distribution. We extend it to the ZIBELLT distribution and evaluated by the following functions.
}

\value{
dzibellt gives the (log) probability function.
pzibellt gives the (log) distribution function.
qzibellt gives the quantile function.
rzibellt generates random values.
mle_zibell gives the maximum likelihood estimates with standard error of
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
\code{\link{mle_zibell}
}
}


\examples{
dzibellt (2, 0.12, 0.2,0.2)

x <- data_sbirth
mle_zibellt (x, 0.15,1.8,1.05)

# Plot of PMF of the ZIBELLT distribution for selected parameteric values.
# The ZIBELLT distribution reduces to the ZIBELL distribution for fixing theta=1.

x <- 0:30
lambda <- 1.2
theta <- 1.2
pi <- 0.01
plot(dzibellt(x, lambda, theta, pi), type = "h", col="red",lwd = 3,
     main = "Zero inflated Bell Touchard distribution",
     ylab = "P(X = x)", xlab = "Number of events")

legend("topright",c(expression(lambda==2~~theta==1.5~~pi==0.01)),
lty=1, col="red", lwd=2,cex=1.0)
}



