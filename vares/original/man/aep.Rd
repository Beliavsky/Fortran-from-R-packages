\name{aep}
\alias{daep}
\alias{paep}
\alias{varaep}
\alias{esaep}
\title{Asymmetric exponential power distribution}
\description{Computes the pdf, cdf, value at risk and expected shortfall for the asymmetric exponential power distribution due to Zhu and Zinde-Walsh (2009) given by
\deqn{\begin{array}{ll}
&\displaystyle
f (x) =  \left\{ \begin{array}{ll}
\displaystyle
\frac {\alpha}{\alpha^{*}} K \left( q_1 \right) \exp \left[ -\frac {1}{q_1} \left | \frac {x}{2 \alpha^{*}} \right |^{q_1} \right], & \mbox{if $x \leq 0$,}
\\
\\
\displaystyle
\frac {1 - \alpha}{1 - \alpha^{*}} K \left( q_2 \right) \exp \left[ -\frac {1}{q_2} \left | \frac {x}{2 - 2 \alpha^{*}} \right |^{q_2} \right], & \mbox{if $x > 0$,}
\end{array}
\right.
\\
&\displaystyle
F (x) = \left\{ \begin{array}{ll}
\displaystyle
\alpha Q \left( \frac {1}{q_1} \left( \frac {\mid x \mid}{2 \alpha^{*}} \right)^{q_1}, \frac {1}{q_1} \right), & \mbox{if $x \leq 0$,}
\\
\\
\displaystyle
1 - (1 - \alpha) Q \left( \frac {1}{q_2} \left( \frac {\mid x \mid}{2 - 2 \alpha^{*}} \right)^{q_2}, \frac {1}{q_2} \right), & \mbox{if $x > 0$,}
\end{array}
\right.
\\
&\displaystyle
{\rm VaR}_p (X) = \left\{ \begin{array}{ll}
\displaystyle
-2 \alpha^{*} \left[ q_1 Q^{-1} \left( \frac {p}{\alpha}, \frac {1}{q_1} \right) \right]^{\frac {1}{q_1}}, & \mbox{if $p \leq \alpha$,}
\\
\\
\displaystyle
2 \left(1 - \alpha^{*}\right) \left[ q_2 Q^{-1} \left( \frac {1 - p}{1 - \alpha}, \frac {1}{q_2} \right) \right]^{\frac {1}{q_2}}, & \mbox{if $p > \alpha$,}
\end{array}
\right.
\\
&\displaystyle
{\rm ES}_p (X) = \left\{ \begin{array}{ll}
\displaystyle
-\frac {2 \alpha^{*}}{p} \int_0^p \left[ q_1 Q^{-1} \left( \frac {v}{\alpha}, \frac {1}{q_1} \right) \right]^{\frac {1}{q_1}} dv, & \mbox{if $p \leq \alpha$,}
\\
\\
\displaystyle
-\frac {2 \alpha^{*}}{p} \int_0^\alpha \left[ q_1 Q^{-1} \left( \frac {v}{\alpha}, \frac {1}{q_1} \right) \right]^{\frac {1}{q_1}} dv & \
\\
\quad
\displaystyle
+\frac {2 \left(1 - \alpha^{*}\right)}{p} \int_\alpha^p \left[ q_2 Q^{-1} \left( \frac {1 - v}{1 - \alpha}, \frac {1}{q_2} \right) \right]^{\frac {1}{q_2}} dv, & \mbox{if $p > \alpha$}
\end{array}
\right.
\end{array}}
for \eqn{-\infty < x < \infty}, \eqn{0 < p < 1}, \eqn{0 < \alpha < 1}, the scale parameter, \eqn{q_1 > 0}, the first shape parameter, and \eqn{q_2 > 0}, the second shape parameter, where \eqn{\alpha^{*} = \alpha K \left( q_1 \right) / \left\{ \alpha K \left( q_1 \right) + (1 - \alpha) K \left( q_2 \right) \right\}}, \eqn{K (q) = \frac {1}{2 q^{1/q} \Gamma (1 + 1/q)}}, \eqn{Q (a, x) = \int_x^\infty t^{a - 1} \exp \left( -t \right) dt / \Gamma (a)} denotes the regularized complementary incomplete gamma function, \eqn{\Gamma (a) = \int_0^\infty t^{a - 1} \exp \left( -t \right) dt} denotes the gamma function, and \eqn{Q^{-1} (a, x)} denotes the inverse of \eqn{Q (a, x)}.}
\usage{daep(x, q1=1, q2=1, alpha=0.5, log=FALSE)
paep(x, q1=1, q2=1, alpha=0.5, log.p=FALSE, lower.tail=TRUE)
varaep(p, q1=1, q2=1, alpha=0.5, log.p=FALSE, lower.tail=TRUE)
esaep(p, q1=1, q2=1, alpha=0.5)}
\arguments{
  \item{x}{scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{p}{scaler or vector of values at which the value at risk or expected shortfall needs to be computed}
  \item{alpha}{the value of the scale parameter, must be in the unit interval, the default is 0.5}
  \item{q1}{the value of the first shape parameter, must be positive, the default is 1}
  \item{q2}{the value of the second shape parameter, must be positive, the default is 1}
  \item{log}{if TRUE then log(pdf) are returned}
  \item{log.p}{if TRUE then log(cdf) are returned and quantiles are computed for exp(p)}
  \item{lower.tail}{if FALSE then 1-cdf are returned and quantiles are computed for 1-p}}
\value{An object of the same length as \code{x}, giving the pdf or cdf values computed at \code{x} or an object of the same length as \code{p}, giving the values at risk or expected shortfall computed at \code{p}.}
\references{Stephen Chan, Saralees Nadarajah & Emmanuel Afuecheta (2016). An R Package for Value at Risk and Expected Shortfall, Communications in Statistics - Simulation and Computation, 45:9, 3416-3434, \doi{10.1080/03610918.2014.944658}}
\author{Saralees Nadarajah}
\examples{x=runif(10,min=0,max=1)
daep(x)
paep(x)
varaep(x)
esaep(x)}
 
