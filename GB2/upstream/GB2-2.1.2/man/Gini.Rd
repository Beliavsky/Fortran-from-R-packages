\name{Gini}
\Rdversion{2.1}
\alias{Gini}
\alias{gini.gb2}
\alias{gini.b2}
\alias{gini.dag}
\alias{gini.sm}
\title{
Computation of the Gini Coefficient for the GB2 Distribution and its Particular Cases.
}
\description{
Computes the Gini coefficient for the GB2 distribution using the function \code{\link{gb2.gini}}. 
Computes the Gini coefficient for the Beta Distribution of the Second Kind, Dagum and Singh-Maddala distributions.
}
\usage{
gini.gb2(shape1, shape2, shape3)
gini.b2(shape2, shape3)
gini.dag(shape1, shape2)
gini.sm(shape1, shape3)
}
\arguments{
  \item{shape1}{numeric; positive parameter.}
  \item{shape2, shape3}{numeric; positive parameters of the Beta distribution.}
}
\value{
The Gini coefficient.
}
\references{
 Kleiber, C. and Kotz, S. (2003) 
 \emph{Statistical Size Distributions in Economics and Actuarial Sciences}, chapter 6.
 Wiley, Ney York.
 
 McDonald, J. B. (1984)
 Some generalized functions for the size distribution of income. 
 \emph{Econometrica}, \bold{52}, 647--663.
}
\author{
Monique Graf
}
\seealso{
\code{\link{gb2.gini}}
}
%\examples{
%}
\keyword{distribution }

