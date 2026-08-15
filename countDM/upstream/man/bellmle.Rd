\name{Bell}
\alias{bell_mle}
\alias{mle.bell}
\title{MLE of the Bell distribtion}
\description{Evaluates the maximum likelihood estimate of the Bell distribtion.
The PMF of the Bell distribution is as follows:
\deqn{
f(X=x\mid\theta)=\frac{\theta^{x}e^{e^{\theta}+1}B_{x}}{x!};\qquad x=0,1,2,\,\dots,
}
where \eqn{\theta>0} denotes the Bell parameter and \eqn{B_{x}} is the Bell number and it is given by
\deqn{B_{n}=\frac{1}{e}\sum_{k=0}^{\infty}\frac{k^{n}}{k!}.}
The Bell number \eqn{B_{n}} in the above equation is the nth moment of the Poisson distribution with parameter equal to 1.
}
\usage{
bell_mle (x)
mle.bell (x, theta)
}
\arguments{
\item{x}{
A vector of (non-negative integer) discrete values.
}
\item{theta}{
A vector of (non-negative integer) values.
}
}

\details{
The function allows to estimate the unknown parameter of the Bell distribution with loglik value using a Newton-Raphson
algorithm.
}

\value{
bell_mle gives the maximum liklihood estimate of parameter theta.
loglik gives value of the maximised log-likelihood. The mle.bell gives the maximum liklihood estimate with standard error and AIC,
}

\references{
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
\code{\link{mle_borel}, \link{mle_poisson}
}
}

\examples{
x <- data_sbirth
bell_mle (x)
mle.bell (x, 1.2)
}



