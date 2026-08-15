\name{Touchard polynomials}
\alias{TP}
\title{Touchard polynomials}
\description{The function allows to compuate the Touchard polynomial. It is mathematically defined by
\deqn{
T_{x}\left(\theta\right)=\frac{1}{e^{\theta}}\sum_{k=0}^{\infty}\frac{k^{x}}{k!}\theta^{k}.
}
The first few Touchard polynomials are as follows:
\deqn{\begin{cases}
\begin{array}{ccccccccc}
T_{0}\left(\theta\right) & = & 1\\
T_{1}\left(\theta\right) & = & \theta\\
T_{2}\left(\theta\right) & = & \theta^{2} & + & \theta\\
T_{3}\left(\theta\right) & = & \theta^{3} & + & 3\theta^{2} & + & \theta\\
T_{4}\left(\theta\right) & = & \theta^{4} & + & 6\theta^{3} & + & 7\theta^{2} & + & \theta.
\end{array}\end{cases}.
}
}
\usage{
TP (x, theta)
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
The function allows to provide the Touchard polynomials.
}

\value{
TP gives the Touchard polynomials after specifying parameteric value.
}

\references{
Castellares, F., Lemonte, A. J., & Moreno–Arenas, G. (2020). On the two-parameter Bell–Touchard discrete distribution. Communications in Statistics-Theory and Methods, 49(19), 4834-4852.
}

\author{
Muhammad Imran and M.H. Tahir.

R implementation and documentation: Muhammad Imran \email{imranshakoor84@yahoo.com} and M.H. Tahir <mht@iub.edu.pk>.
}


%\note{
%%  ~~further notes~~
%}



\examples{
TP (2,3)
}



