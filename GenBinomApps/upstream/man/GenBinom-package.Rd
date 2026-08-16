\name{GenBinomApps-package}
\alias{GenBinomApps-package}
\alias{GenBinomApps}
\docType{package}
\title{
Clopper-Pearson Confidence Interval and Generalized Binomial Distribution
}
\description{
Density, distribution function, quantile function, and random generation for the Generalized Binomial Distribution. Also included are functions to compute the Clopper-Pearson confidence interval limits for the standard case, for an enhanced model, and the required sample size for a given target probability for both models.
}
\details{
This package originates from semiconductor manufacturing but can also be used for other purposes. The functions are based on the paper \emph{Decision-Theoretical Model for Failures which are Tackled by Countermeasures}, Kurz et al. (2014).


The generalized binomial distribution is defined as the sum of independent, not identically binomial distributed random variables. That means, they have different success probabilities, and they can have different sample sizes.


\emph{Example:} A person has to drive \eqn{3} routes at each working day. The probabilities for a radar control on  these routes are  \eqn{0.1\%}{0.1\%} for the first route,  \eqn{0.5\%}{0.5\%} for the second route and  \eqn{1\%}{1\%} for the third route. The person has to drive route \eqn{1} and route \eqn{2} one time per day and route \eqn{3} two times per day. What are the probabilities to have \eqn{0,1,2,} more than \eqn{2} controls at \eqn{100} working days? \cr
Knowing that the number of controls is binomially distributed for each route:\cr
\eqn{R_{1}}{R1} ~ \eqn{binom(100,0.001)},  \eqn{R_{2}}{R2} ~ \eqn{binom(100,0.005)}, \eqn{R_{3}}{R3} ~ \eqn{binom(200,0.01)}

Thus the sum of these binomially distributed random variables has a generalized binomial distribution with parameters \eqn{n_{1}=100,n_{2}=100,n_{3}=200,p_{1}=0.001,p_{2}=0.005,p_{3}=0.01}{n1=100,n2=100,n3=200,p1=0.001,p2=0.005,p3=0.01}.\cr
\eqn{R=R_{1}+R_{2}+R_{3}}{R=R1+R2+R3}, \verb{} \eqn{R} ~ \eqn{gbinom(100,100,200,0.001,0.005,0.01)}

In this example the probabilities \eqn{P(R=0)}, \eqn{P(R=1)}, \eqn{P(R=2)}, \eqn{P(R > 2)} can be computed straightforward. \cr
See the examples for the results.


Consider now a burn-in study in which \eqn{k} failures are observed. The number of failures is binomially distributed. Thus, the Clopper-Pearson confidence interval limits can be used to obtain a confidence interval for the failure probability. 
If failures occur, countermeasures should be implemented with a type specific effectivity. Consider the case of different failure types. That leads to more than one countermeasure. Each countermeasure can have a different effectivity. The probability for solving a certain number of failures can be computed with the generalized binomial distribution. It gives the likelihoods for various possible outcome scenarios, if the countermeasures would have been introduced from the beginning on. Based on the model in Kurz et al. (2014), confidence intervals can be computed.
}


\author{
Horst Lewitschnig, David Lenzi.

Maintainer: 
Horst Lewitschnig <Horst.Lewitschnig@infineon.com>
}


\note{
The generalized binomial distribution described here is also known as Poisson-binomial distribution. 
}


\references{
D.Kurz, H.Lewitschnig, J.Pilz, \emph{Decision-Theoretical Model for Failures which are Tackled by Countermeasures}, IEEE Transactions on Reliability, Vol. 63, No. 2, June 2014.

K.J. Klauer, \emph{Kriteriumsorientierte Tests}, Verlag fuer Psychologie, Hogrefe, 1987, Goettingen, p. 208 ff.

M.Fisz, \emph{Wahrscheinlichkeitsrechnung und mathematische Statistik}, VEB Deutscher Verlag der Wissenschaften, 1973, p. 164 ff.

C.J.Clopper and E.S. Pearson, \emph{The use of confidence or fiducial limits illustrated in the case of the binomial}, Biometrika, vol. 26, 404-413, 1934.
}




\examples{
## n1=100, n2=100, n3=200, p1=0.001, p2=0.005, p3=0.01 
dgbinom(c(0:2),size=c(100,100,200),prob=c(0.001,0.005,0.01))
#  0.07343377 0.19260317 0.25173556
pgbinom(2,size=c(100,100,200),prob=c(0.001,0.005,0.01),lower.tail=FALSE)
# 0.4822275


## n=110000 tested devices, 2 failures divided in 2 failure types k1=1, k2=1.
## 2 countermeasures with effectivities p1=0.5, p2=0.8
cm.clopper.pearson.ci(110000,size=c(1,1), cm.effect=c(0.5,0.8))
# Confidence.Interval = upper
# Lower.limit = 0
# Upper.limit = 3.32087e-05
# alpha = 0.1

## target failure probability p=0.00001, 2 failures divided in 2 failure types k1=1, k2=1.
## 2 countermeasures with effectivities p1=0.5, p2=0.8
cm.n.clopper.pearson(0.00001,size=c(1,1), cm.effect=c(0.5,0.8))
# 365299
}
