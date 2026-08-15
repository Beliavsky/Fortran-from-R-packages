\name{Poisson}
\alias{mle_poisson}
\title{MLE of the Poisson distribution}
\description{Evaluates the MLE of the Poisson distribtion.
The PMF of the Poisson distribution is as follows:
\deqn{
f\left(X=x\mid\,\theta\right)=\frac{\theta^{x}e^{-\theta}}{x!};\qquad x=0,1,2,\dots,
}
where parameter \eqn{\theta>0} and it is equal to the expected or mean value of X and also to its variance.
}
\usage{
mle_poisson (x, theta)
}
\arguments{
\item{x}{
A vector of (non-negative integer) discrete values.
}
\item{theta}{
A vector of (non-negative integer) values, \eqn{\theta>0}.
}
}

\details{
The function allows to estimate the unknown parameter of Poisson distribution with standard error of
the estimate and model selection measure, the Akaike information criterion (AIC).
}

\value{
mle_poisson gives the MLE along with standard error of the estimate and model selection measure AIC.
}

\references{
Frank A. Haight (1967). Handbook of the Poisson Distribution. New York: John Wiley & Sons.
}

\author{
Muhammad Imran and M.H. Tahir.

R implementation and documentation: Muhammad Imran \email{imranshakoor84@yahoo.com} and M.H. Tahir <mht@iub.edu.pk>.
}


%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{mle.bell}, \link{mle_borel}
}
}

\examples{
x <- data_sbirth
mle_poisson (x, 0.2)
}



