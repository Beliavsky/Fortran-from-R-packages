\encoding{utf8} 
\name{glm.hermite}
\alias{glm.hermite}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Maximum likelihood estimation and Hermite regression
}
\description{
  \code{glm.hermite} is used to fit generalized linear models with count 
  responses following a Hermite distribution, specified by giving a symbolic 
  description of the linear predictor. A \code{summary} method providing the 
  most meaningful information on the fitted model is available for objects of 
  class \code{glm.hermite}.
}
\usage{
  glm.hermite(formula, data, link="log", start=NULL, m = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
symbolic description of the model. A typical predictor has the form 
\code{response ~ terms} where response is the (numeric) response vector and 
terms is a series of terms which specifies a linear predictor for response.
}
\item{data}{
an optional data frame containing the variables in the model.
}
\item{link}{
character specification of link function: "log" or "identity". By default 
\code{link="log"}.
}
\item{start}{
a vector containing the starting values for the parameters of the specified 
model. Its default value is \code{NULL}.
}
\item{m}{
value for parameter \code{m}. Its default value is \code{NULL}, and in that 
case it will be estimated inside the function.
}
}

\value{
  \code{glm.hermite} returns an object of class \code{glm.hermite}, which is a 
  list including the following components:

  \item{coefs}{
the vector of coefficients.
}
  \item{data}{
an optional data frame containing the variables in the model.
}
  \item{loglik}{
log-likelihood of the fitted model.
}
  \item{vcov}{
covariance matrix of all coefficients in the model (derived from the Hessian of 
the \code{maxLik} output).
}
  \item{hess}{Hessian matrix, returned by the \code{maxLik} output.
}
  \item{fitted.values}{
the fitted mean values, obtained by transforming the linear predictors by the 
inverse of the link
function.
}
  \item{w}{Likelihood ratio test statistic.
}
  \item{pval}{Likelihood ratio test p-value.
}
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
María Oliveira, Manuel Higueras, David Moriña and Pere Puig

}

\seealso{
\code{\link{Distributions}} for some other distributions, 
\code{\link{qhermite}}, \code{\link{phermite}}, \code{\link{rhermite}}, 
\code{\link{hermite-package}}
}

\examples{
data <- c(rep(0,122), rep(1,40), rep(2,14), rep(3,16), rep(4,6), rep(5,2))
mle1 <- glm.hermite(data~1, link="log", start=NULL, m=3)
mle1
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hermite }
