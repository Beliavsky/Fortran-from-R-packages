\encoding{utf8} 
\name{dhermite}
\alias{dhermite}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Probability mass function for the generalized Hermite distribution
}
\description{
  Probability mass function for the generalized Hermite distribution with 
  parameters \code{a}, \code{b} and \code{m}.
}
\usage{
  dhermite(x, a, b, m=2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  vector of non-negative integer quantiles. 
}
  \item{a}{
first parameter for the Hermite distribution.
}
  \item{b}{
second parameter for the Hermite distribution.
}
  \item{m}{
degree of the generalized Hermite distribution. Its default value is \code{2}, 
corresponding to the standard Hermite distribution.
}
}

\value{
  Probability for a generalized Hermite random varible with parameters \code{a}, 
  \code{b} and \code{m} of taking \code{x} counts.
  
}

\references{
Kemp C D, Kemp A W. Some Properties of the Hermite Distribution. Biometrika 
1965;\bold{52} (3-4):381–394.

McKendrick A G Applications of Mathematics to Medical Problems. Proceedings of 
the Edinburgh Mathematical Society 1926;\bold{44}:98–130. 

Kemp A W, Kemp C D. An alternative derivation of the Hermite distribution. 
Biometrika 1966;\bold{53} (3-4):627–628.

Patel Y C. Even Point Estimation and Moment Estimation in Hermite Distribution. 
Biometrics 1976;\bold{32} (4):865–873.

Gupta R P, Jain G C. A Generalized Hermite distribution and Its Properties. 
SIAM Journal on Applied Mathematics 1974;\bold{27}:359–363.

Bekelis, D. Convolutions of the Poisson laws in number theory. In Analytic & 
Probabilistic Methods in Number Theory: Proceedings of the 2nd International 
Conference in Honour of J. Kubilius, Lithuania 1996;\bold{4}:283–296.

Zhang J, Huang H. On Nonnegative Integer-Valued Lévy Processes and Applications 
in Probabilistic Number Theory and Inventory Policies. American Journal of 
Theoretical and Applied Statistics 2013;\bold{2}:110–121.

Kotz S. Encyclopedia of statistical sciences. John Wiley 1982-1989.

Kotz S. Univariate discrete distributions. Norman L. Johnson 2005.

Puig P. (2003). Characterizing Additively Closed Discrete Models by a Property 
of Their Maximum Likelihood Estimators, with an Application to Generalized 
Hermite Distributions. Journal of the American Statistical Association 2003;
\bold{98}:687–692.

}

\author{
David Moriña, Manuel Higueras, Pedro Puig and María Oliveira


}

\seealso{
\code{\link{Distributions}} for some other distributions, 
\code{\link{qhermite}}, \code{\link{phermite}}, \code{\link{rhermite}}, 
\code{\link{hermite-package}}, \code{\link{glm.hermite}}
}

\examples{
d <- dhermite(3, 0.8, 0.3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hermite }
