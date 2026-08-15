\name{JSB}
\alias{dJSB}
\alias{pJSB}
\alias{rJSB}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The Johnson \eqn{S_B} distribution.
}
\description{
Density, distribution function and random generation for the Johnson \eqn{S_B} distribution.
}
\usage{
dJSB(x, delta=1, log=FALSE)
pJSB(q, delta=1, lower.tail=TRUE, log.p=FALSE)
rJSB(n, delta=1)
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
The Johnson \eqn{S_B} distribution has density
\deqn{
f(x)=\frac{\delta}{x(1-x)}\phi\left(\delta \eta(x)\right), \quad x \in (0,1),
}
where \eqn{\eta(x)=\log(\frac{x}{1-x})}, \eqn{\phi(\cdot)} denotes the density of the
standard normal distribution and \eqn{\delta>0}. Its cumulative distribution function is
\deqn{
F(x)=\Phi\left(\delta \eta(x)\right), \quad x \in (0,1),
}
where \eqn{\Phi(\cdot)} is the cumulative distribution function of the
standard normal distribution.
}
\value{
dJSB gives the density, pJSB gives the distribution function, and rJSB generates random deviates.
The length of the result is determined by n for rasin, and is the maximum of the lengths of the numerical arguments for the other functions.
The numerical arguments other than n are recycled to the length of the result. Only the first elements of the logical arguments are used.
}
\references{
Kotz, S., van Dorp, J.R. (2004). Beyond Beta. Other Continuous Families of Distributions with
Bounded Support and Applications. World Scientific.
}
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
dJSB(0.5, 1.2)
pJSB(0.5, 0.5)
rJSB(5, 1.5)
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
