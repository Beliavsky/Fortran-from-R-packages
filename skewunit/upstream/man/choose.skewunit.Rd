\name{choose.skewunit}
\alias{choose.skewunit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Choose a Distribution in a Family of Skew Distributions with Bounded Support
}
\description{
choose.skewunit selects a combination of \eqn{f} and \eqn{G} in a Family of Skew Distributions with Bounded Support
based on the Akaike information criterion (AIC) or Bayesian information criterion (BIC).
}
\usage{
choose.skewunit(x, criteria="AIC")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
data in \eqn{(0,1)} interval.
}
  \item{criteria}{
criteria to choose a model: AIC (default) or BIC.
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
Therefore, also are included the baseline models asin, Uquad, triang, JSB and sbeta in the 
model selection.
%
}
\value{
an object of class "skewunit" is returned. The object returned for this functions is a list containing the following
components:
\item{summary}{a data frame with all the fitted models, including the specified criteria and number of parameters.}
\item{coefficients}{estimates and standard errores for the model with the lowest criteria.}
\item{logLik, AIC, BIC}{the corresponding criteria for the model with the lowest criteria.}
\item{family1, family2}{the specification for family1 and family2 for the model with the lowest criteria.}
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
x=rskewunit(100, lambda=-0.5, delta=1.2, family1="asin", family2="triang")
aux=choose.skewunit(x, criteria="AIC")
aux
aux$summary
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
