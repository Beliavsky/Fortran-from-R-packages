\name{Zero inflated Poisson}
\alias{mle_zip}
\title{MLE of the zero inflated Poisson distribtion}
\description{Evaluates the MLE of the zero inflated Poisson (ZIP) distribtion.
The PMF of the ZIP is as follows:
\deqn{
f\left(X=x\mid\alpha,\,\theta\right)=\begin{cases}
\alpha+\left(1-\alpha\right)e^{-\theta}, & x=0\\
\left(1-\alpha\right)\frac{\theta^{x}e^{-\theta}}{x!}, & x=1,2\cdots,
\end{cases}
}
where \eqn{\alpha\in(0,1)} denotes the probability of extra zeros and \eqn{\theta>0} is a Poisson parameter, which is also its mean and variance.
}
\usage{
mle_zip (x, alpha, theta)
}
\arguments{
\item{x}{
A vector of (non-negative integer) discrete values.
}
\item{theta}{
A vector of (non-negative integer) values, \eqn{\theta>0}.
}
\item{alpha}{
A vector of (non-negative integer) values, \eqn{\alpha\in(0,1)}.
}
}

\details{
The function allows to estimate the unknown parameter of the ZIP distribution with standard error of the estimate and model selection measure, the Akaike information criterion (AIC).
}

\value{
mle_zip gives the MLE along with standard error of the estimate and model selction measure AIC.
}

\references{
Beckett, S., Jee, J., Ncube, T., Pompilus, S., Washington, Q., Singh, A., & Pal, N. (2014). Zero-inflated Poisson (ZIP) distribution: Parameter estimation and applications to model data from natural calamities. Involve, a Journal of Mathematics, 7(6), 751-767.
}

\author{
Muhammad Imran and M.H. Tahir.

R implementation and documentation: Muhammad Imran \email{imranshakoor84@yahoo.com} and M.H. Tahir <mht@iub.edu.pk>.
}



%\note{
%%  ~~further notes~~
%}


\seealso{
\code{\link{mle_zibell}, \link{mle_zibellt}
}
}

\examples{
x <- data_sbirth
mle_zip (x, 0.2, 1.5)
}



