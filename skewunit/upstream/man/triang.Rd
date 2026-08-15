\name{triang}
\alias{dtriang}
\alias{ptriang}
\alias{rtriang}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The triangular distribution
}
\description{
Density, distribution function and random generation for the triangular distribution.
}
\usage{
dtriang(x, log=FALSE)
ptriang(q, lower.tail=TRUE, log.p=FALSE)
rtriang(n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, q}{vector of quantiles.
}
  \item{n}{
number of observations. If length(n) > 1, the length is taken to be the number required.
}
\item{log, log.p}{
logical; if TRUE, probabilities p are given as log(p).
}
\item{lower.tail}{
logical; if TRUE (default), probabilities are \eqn{P[X\leq x]}, otherwise, \eqn{P[X>x]}.
}
}
\details{
The triangular distribution has density
\deqn{
f(x) = \left\{
\begin{array}{lr}
4x, & 0\leq x\leq 1/2,\\
4(1-x), & 1/2<x\leq 1,
\end{array}
\right.
}
and cumulative distribution function
\deqn{
F(x) = \left\{
\begin{array}{lr}
2x^2, & 0\leq x\leq 1/2,\\
2x^2-(2x-1)^2, & 1/2<x\leq 1,
\end{array}
\right.
}
}
\value{
dtriang gives the density, ptriang gives the distribution function, and rtriang generates random deviates.
The length of the result is determined by n for rasin, and is the maximum of the lengths of the numerical arguments for the other functions.
The numerical arguments other than n are recycled to the length of the result. Only the first elements of the logical arguments are used.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Diego Gallardo
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
dtriang(0.5)
ptriang(0.5)
rtriang(5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
