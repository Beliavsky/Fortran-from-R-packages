\name{Borel}
\alias{mle_borel}
\title{MLE of the Borel distribution}
\description{Evaluates the MLE of the Borel distribution. It is defined by the following PMF:
\deqn{
f(X=x\mid \alpha)=\frac{\left(\alpha x\right)^{x-1}e^{-\alpha x}}{x!};\qquad x=1,2,\dots,
}
where the parameter \eqn{\alpha\in(0,1)}.
}
\usage{
mle_borel (x, alpha)
}
\arguments{
\item{x}{
A vector of (non-negative integer) discrete values.
}
\item{alpha}{
A vector of (non-negative integer) values, \eqn{\alpha\in(0,1)}.
}
}

\details{
The function allows to estimate the unknown parameter of the Borel distribution with standard error of the estimate and model selection measure, the Akaike information criterion (AIC).
}

\value{
mle_borel gives the MLE along with standard error of the estimate and model selction measure AIC.

}

\references{
Tanner, J. C. (1961). A derivation of the Borel distribution. Biometrika, 48(1/2), 222-224.
}

\author{
Muhammad Imran and M.H. Tahir.

R implementation and documentation: Muhammad Imran \email{imranshakoor84@yahoo.com} and M.H. Tahir <mht@iub.edu.pk>.
}



%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{mle.bell}, \link{mle_poisson}
}
}

\examples{
x <- c(rep(1,48), rep(2,20), rep(3,7), rep(4,5),rep(5,2),rep(6,6))
mle_borel (x, 0.8)
}



