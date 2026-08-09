\name{wrapnlsr}
\alias{wrapnlsr}
\title{wrapnlsr}
\description{
Provides class nls solution to a nonlinear least squares solution 
using the Nash Marquardt tools.
}
\usage{
wrapnlsr(formula = NULL,  data = NULL,  start = NULL,  control = NULL,
    trace = FALSE, subset = NULL, lower = -Inf, upper = Inf,  weights = NULL,
    ...)                           
}
\arguments{
\item{formula}{The modeling formula. Looks like 'y~b1/(1+b2*exp(-b3*T))'}

\item{data}{a data frame containing data for variables
used in the formula that are NOT the parameters. This
data may also be defined in the parent frame i.e.,
'global' to this function}

\item{start}{MUST be a named vector with all elements present
e.g., start=c(b1=200, b2=50, b3=0.3)}

\item{control}{a list of control parameters. See nlsr.control().}

\item{trace}{TRUE for console output during execution (default FALSE)}

\item{subset}{an optional vector specifying a subset of observations 
to be used in the fitting process. NOT used currently by nlxb()
or nlfb() and will throw an error if present and not NULL.}

\item{lower}{a vector of lower bounds on the parameters. 
If a single number, this will be applied to all parameters
Default \code{-Inf}.}

\item{upper}{a vector of upper bounds on the parameters. If a single number, 
this will be applied to all parameters. Default \code{Inf}.}

\item{weights}{A vector of (usually fixed) weights. The objective function that will be 
minimized is the sum of squares where each residual is multiplied by the 
square root of the corresponding weight. Default \code{NULL} implies 
unit weights.}

\item{...}{additional data needed to evaluate the modeling functions}
}
\value{
A solution object of type \code{nls}
}
\examples{
library(nlsr)
cat("kvanderpoel.R test of wrapnlsr\n")
x<-c(1,3,5,7)
y<-c(37.98,11.68,3.65,3.93)
pks28<-data.frame(x=x,y=y)
fit0<-try(nls(y~(a+b*exp(1)^(-c*x)), data=pks28, start=c(a=0,b=1,c=1), 
              trace=TRUE))
print(fit0)
fit1<-nlxb(y~(a+b*exp(-c*x)), data=pks28, start=c(a=0,b=1,c=1), trace = TRUE)
print(fit1) 
cat("\n\n or better\n")
fit2<-wrapnlsr(y~(a+b*exp(-c*x)), data=pks28, start=c(a=0,b=1,c=1), 
               lower=-Inf, upper=Inf, trace = TRUE)
fit2

weed <- c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443,
          38.558, 50.156, 62.948, 75.995, 91.972)
tt <- 1:12
weeddf <- data.frame(tt, weed)
hobbsu <- weed ~ b1/(1+b2*exp(-b3*tt))
st2 <- c(b1=200, b2=50, b3=0.3)
wts <- 0.5^tt # a straight scaling comes via wts <- rep(0.01, 12)
lo <- c(200, 0, 0)
up <- c(1000, 1000, 1000)
whuw2 <-  try(wrapnlsr(start=st2, formula=hobbsu, data=weeddf, subset=2:11,
                  weights=wts, trace=TRUE, lower=lo, upper=up))
summary(whuw2)
deviance(whuw2)
whuw2a <-  try(nlsr(start=st2, formula=hobbsu, data=weeddf, subset=2:11,
                  weights=wts, trace=TRUE, lower=lo, upper=up))
summary(whuw2a)
deviance(whuw2a)
}
