\name{estimate.skewunit}
\alias{estimate.skewunit}
\alias{print.skewunit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimation for a Family of Skew Distributions with Bounded Support
}
\description{
Perform parameter estimation for a family of skew distributions with bounded support.
}
\usage{
estimate.skewunit(x, family1 = "asin", family2 = "asin", est.var = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
data in \eqn{(0,1)} interval.
}
  \item{family1}{
first family of distributions related to \eqn{f} (asin by default). See details Section.
}
  \item{family2}{
second family of distributions related to \eqn{G} (asin by default, NULL also can be specified). See details Section.
}
  \item{est.var}{
logical; if TRUE, estimate the standard errors of the estimators.
}
}
\details{
The Family of Skew Distributions with Bounded Support is defined by its density function given by
\deqn{
f(x)=2 G(\lambda(x-0.5)+0.5), \quad x \in (0,1), \lambda \in (-1,1),
}
where \eqn{f} is symmetric around 0.5, i.e., \eqn{f(x-0.5)=f(x+0.5)}.
The available options for family1 and family2 are asin, Uquad, triang, JSB and sbeta.
When family2 is NULL, f corresponds to the model specified in family1.
%
}
\value{
an object of class "skewunit" is returned. The object returned for this functions is a list containing the following
components:
\item{x}{x}
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Diego Gallardo, Emilio Gomez-Deniz, Osvaldo Venegas and Hector W. Gomez
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
set.seed(2100)
x=rskewunit(100, lambda=-0.5, delta=1.2, family1="asin", family2="JSB")
estimate.skewunit(x, family1="asin", family2="JSB")
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
