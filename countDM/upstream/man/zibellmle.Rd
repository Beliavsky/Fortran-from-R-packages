\name{Zero inflated Bell}
\alias{mle_zibell}
\title{MLE of the zero inflated Bell distribtion}
\description{Evaluates the MLE of the zero inflated Bell (ZIBELL) distribtion. The PMF of the ZIBELL distribution is as follows:
\deqn{
f\left(X=x\mid\alpha,\,\lambda\right)=\begin{cases}
\alpha+\left(1-\alpha\right)\exp\left\{ \theta\left[1-e^{\lambda}\right]\right\} , & x=0\\
\left(1-\alpha\right)\exp\left\{ \theta\left[1-e^{\lambda}\right]\right\} \frac{\lambda^{x}\,B_{x}}{x!}, & x=1,2\cdots,
\end{cases}
}
where \eqn{\alpha\in(0,1)},  \eqn{\lambda>0} and \eqn{B_{x}} are the Bell numbers and it is given by
\deqn{B_{n}=\frac{1}{e}\sum_{k=0}^{\infty}\frac{k^{n}}{k!}.}
}
\usage{
mle_zibell (x, alpha, lambda)
}
\arguments{
\item{x}{
A vector of (non-negative integer) values, discrete values.
}
\item{lambda}{
A vector of (non-negative integer) values, \eqn{\lambda>0}.
}
\item{alpha}{
A vector of (non-negative integer) values, \eqn{\alpha\in(0,1)}.
}
}

\details{
The function allows to estimate the unknown parameter of the ZIBELL distribution with standard error of the estimate and model selection measure, the Akaike information criterion (AIC).
}

\value{
mle_zibell gives the MLE along with standard error of the estimate and model selction measure AIC.
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
\code{\link{mle_zibellt}
}
}


\examples{
x <- data_sbirth
mle_zibell (x, 0.2, 1.5)
}



