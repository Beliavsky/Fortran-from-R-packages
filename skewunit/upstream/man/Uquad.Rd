\name{Uquad}
\alias{dUquad}
\alias{pUquad}
\alias{rUquad}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The U-quadratic distribution
}
\description{
Density, distribution function and random generation for the U-quadratic distribution.
}
\usage{
dUquad(x, a=0, b=1, log=FALSE)
pUquad(q, a=0, b=1, lower.tail=TRUE, log.p=FALSE)
rUquad(n, a=0, b=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, q}{vector of quantiles.
}
  \item{n}{
number of observations. If length(n) > 1, the length is taken to be the number required.
}
  \item{a, b}{
range of variable x. (\eqn{a<b}).
}
\item{log, log.p}{
logical; if TRUE, probabilities p are given as log(p).
}
\item{lower.tail}{
logical; if TRUE (default), probabilities are \eqn{P[X\leq x]}, otherwise, \eqn{P[X>x]}.
}
}
\details{
The U-quadratic distribution has density
\deqn{
f(x) = \alpha (x-\beta)^2, \quad x\in (a,b), a\leq x \leq b,
}
where \eqn{\alpha=12/(b-a)^3} and \eqn{\beta=(a+b)/2}. Its cumulative distribution function is
\deqn{
F(x) = \frac{\alpha}{3}[(x-\beta)^3+(\beta-a)^3], \quad x\in (a,b).
}
}
\value{
dUquad gives the density, pUquad gives the distribution function, and rUquad generates random deviates.
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
dUquad(0.5)
pUquad(0.5)
rUquad(5)
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
