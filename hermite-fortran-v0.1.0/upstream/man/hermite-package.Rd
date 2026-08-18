\encoding{utf8} 
\name{hermite-package}
\alias{hermite-package}
\alias{hermite}
\docType{package}
\title{
Generalized Hermite distribution
}
\description{
Probability mass, distribution and quantile functions; random generation; 
and regression models for the generalized Hermite distribution.

}
\details{
\tabular{ll}{
Package: \tab hermite\cr
Type: \tab Package\cr
Version: \tab 1.1.2\cr
Date: \tab 2017-05-17\cr
License: \tab GPL version 2 or newer\cr
LazyLoad: \tab yes\cr
}
The package implements probability mass function \code{\link{dhermite}}, 
distribution function \code{\link{phermite}}, quantile function 
\code{\link{qhermite}} and random generation \code{\link{rhermite}} for the 
generalized Hermite distribution. The probability mass function is usually 
parametrized in terms of the mean \eqn{\mu} and the index of dispersion 
\eqn{d = \frac{\sigma^2}{\mu}}:

\eqn{P(X=x) = P(X=0) \frac{\mu^x (m-d)^x}{(m-1)^x} \sum_{j=0}^{[x/m]} 
\frac{(d-1)^j (m-1)^{(m-1)j}}{m^j \mu^{(m-1)j} (m-d)^{mj} (x-mj)!j!}}
where \eqn{P(X=0) = exp(\mu (-1+ \frac{d-1}{m}))}, \code{m} is the degree of 
the generalized Poisson distribution and \eqn{[x/m]} is the integer part of
\eqn{x/m}.

The package is able to fit Hermite regression models as well, by means of the 
function  \code{\link{glm.hermite}}, also in the presence of covariates.
}
\author{
David Moriña, Manuel Higueras, Pedro Puig and María Oliveira

Mantainer: David Moriña Soler <david.morina@uab.cat>
}

\seealso{
\code{\link{Distributions}} for some other distributions, 
\code{\link{qhermite}}, \code{\link{phermite}}, \code{\link{rhermite}}, 
\code{\link{hermite-package}}, \code{\link{glm.hermite}}
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
\keyword{ hermite }

