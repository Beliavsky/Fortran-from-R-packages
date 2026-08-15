\name{Zero one inflated Poisson}
\alias{mle_zoip}
\title{MLE of the zero one inflated Poisson distrubution}
\description{Evaluates the MLE of the zero one inflated Poisson (ZOIP)distribtion.
\deqn{
f\left(X=x\mid\alpha,\,\beta,\,\theta\right)=\begin{cases}
\alpha+\left(1-\alpha-\beta\right)e^{-\theta}, & x=0\\
\beta+\left(1-\alpha-\beta\right)\theta e^{-\theta}, & x=1\\
\left(1-\alpha-\beta\right)e^{-\theta}\frac{\theta^{x}\,}{x!}, & x=2,3,\dots
\end{cases}
}
where \eqn{\theta>0}, the two parameters \eqn{\alpha\in(0,1)} and \eqn{\beta\in(0,1)}, respectively denotes the unknown proportion for incorporating extra zeros and extra ones than those allowed by the traditional Poisson distribution.
}

\usage{
mle_zoip (x, alpha, beta, theta)
}
\arguments{
\item{x}{
A vector of (non-negative integer) values, discrete values.
}
\item{alpha}{
A vector of (non-negative integer) values, \eqn{\alpha\in(0,1)}.
}
\item{beta}{
A vector of (non-negative integer) values, \eqn{\beta\in(0,1)}.
}
\item{theta}{
A vector of (non-negative integer) values, \eqn{\theta>0}.
}
}

\details{
The function allows to estimate the unknown parameter of the ZOIP distribution with standard error of the estimate and model selection measure, the Akaike information criterion (AIC).
}

\value{
mle_zoip gives the MLE along with standard error of the estimate and model selction measure AIC.
}

\references{
Zhang, C., Tian, G. L., & Ng, K. W. (2016). Properties of the zero-and-one inflated Poisson distribution and likelihood-based inference methods. Statistics and its interface, 9(1), 11-32.

Tang, Y., Liu, W., & Xu, A. (2017). Statistical inference for zero-and-one-inflated Poisson models. Statistical Theory and Related Fields, 1(2), 216-226.

Alshkaki, R. S. A. (2016). On the zero-one inflated Poisson distribution. Int J Stat Distrib Appl, 2(4), 42-8.
}

\author{
Muhammad Imran and M.H. Tahir.

R implementation and documentation: Muhammad Imran \email{imranshakoor84@yahoo.com} and M.H. Tahir <mht@iub.edu.pk>.
}



%\note{
%%  ~~further notes~~
%}


\seealso{
\code{\link{mle_zoibell}
}
}

\examples{
x <- data_sbirth
mle_zoip (x, 0.2,0.1, 0.5)
}



