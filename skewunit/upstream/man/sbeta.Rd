\name{sbeta}
\alias{dsbeta}
\alias{psbeta}
\alias{rsbeta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The symmetrical beta distribution.
}
\description{
Density, distribution function and random generation for the symmetrical beta distribution.
}
\usage{
dsbeta(x, delta=1, log=FALSE)
psbeta(q, delta=1, lower.tail=TRUE, log.p=FALSE)
rsbeta(n, delta=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, q}{vector of quantiles.
}
  \item{n}{
number of observations. If length(n) > 1, the length is taken to be the number required.
}
  \item{delta}{
shape parameter (by default is 1).
}
\item{log, log.p}{
logical; if TRUE, probabilities p are given as log(p).
}
\item{lower.tail}{
logical; if TRUE (default), probabilities are \eqn{P[X\leq x]}, otherwise, \eqn{P[X>x]}.
}
}
\details{
The symmetrical beta distribution has density
\deqn{
f(x)=\frac{1}{B(\delta,\delta)}x^{\delta-1}(1-x)^{\delta-1}, \quad x \in (0,1), \delta>0,
}
where \eqn{B(a,b)} denotes the beta function. Its cumulative distribution function is
\deqn{
F(x)=I_x(\delta,\delta), \quad x \in (0,1).
}
}
\value{
dsbeta gives the density, psbeta gives the distribution function, and rsbeta generates random deviates.
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
dsbeta(0.5, 1.2)
psbeta(0.5, 0.5)
rsbeta(5, 1.5)
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
