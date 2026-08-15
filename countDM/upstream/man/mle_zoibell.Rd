\name{Zero one inflated Bell}
\alias{mle_zoibell}
\title{MLE of the zero one inflated Bell distrubution}
\description{Evaluates the MLE of the zero one inflated Bell (ZOIBELL)distribtion.
\deqn{
f\left(X=x\mid\alpha,\,\beta,\,\theta\right)=\begin{cases}
\alpha+\left(1-\alpha-\beta\right)\exp\left(1-e^{\theta}\right), & x=0\\
\beta+\left(1-\alpha-\beta\right)\theta\exp\left(1-e^{\theta}\right), & x=1\\
\left(1-\alpha-\beta\right)\exp\left(1-e^{\theta}\right)\frac{\theta^{x}\,B_{x}}{x!}, & x=2,3,\dots,
\end{cases}
}
where \eqn{\theta>0}, the two parameters \eqn{\alpha\in(0,1)} and \eqn{\beta\in(0,1)}, respectively, denotes the unknown proportion for incorporating extra zeros and extra ones than those allowed by the traditional Bell distribution.
}
\usage{
mle_zoibell (x, alpha, beta, theta)
}
\arguments{
\item{x}{
A vector of (non-negative integer) values, discrete values.
}
\item{theta}{
A vector of (non-negative integer) values, \eqn{\theta>0}.
}
\item{alpha}{
A vector of (non-negative integer) values, \eqn{\alpha\in(0,1).}
}
\item{beta}{
A vector of (non-negative integer) values, \eqn{\beta\in(0,1).
}
}
}

\details{
Inflated models have become quite popular in the recent applied statistical literature.
In many scientific studies, we often experience situations, the data consists of a large proportion of zeros and ones.
Castellares et al. proposed a single-parameter discrete Bell distribution. We extended the Bell distribution into zero one inflated Bell distribution.
}

\value{
mle_zoibell gives the MLE along with standard error of the estimate and model selction measure AIC.

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
\code{\link{mle_zoip}
}
}


\examples{
x <- data_sbirth
mle_zoibell (x, 0.1,0.2,0.2)
}



