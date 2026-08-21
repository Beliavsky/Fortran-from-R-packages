\name{Bivariate Poisson lognormal}
\alias{dbipoilog}
\alias{rbipoilog}

\title{ Bivariate Poisson Lognormal Distribution }
\description{
  Density and random generation for the for the bivariate Poisson 
  lognormal distribution with parameters \code{mu1}, \code{mu2}, 
  \code{sig1}, \code{sig2} and \code{rho}.
}
\usage{
dbipoilog(n1, n2, mu1, mu2, sig1, sig2, rho)
rbipoilog(S, mu1, mu2, sig1, sig2, rho, nu1=1, nu2=1, 
          condS=FALSE, keep0=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n1}{ vector of observed individuals for each species in sample 1}
  \item{n2}{ vector of observed individuals for each species in sample 2 \cr(in the same order as in sample 1) }
  \item{mu1}{ mean of lognormal distribution in sample 1 }
  \item{mu2}{ mean of lognormal distribution in sample 1 }
  \item{sig1}{ standard deviation of lognormal distribution in sample 1 }
  \item{sig2}{ standard deviation of lognormal distribution in sample 2 }
  \item{rho}{ correlation among samples }
  \item{S}{ Total number of species in both communities }
  \item{nu1}{ sampling intensity sample 1 }
  \item{nu2}{ sampling intensity sample 2 }
  \item{condS}{ logical; if TRUE random deviates are conditional on S }
  \item{keep0}{ logical; if TRUE species with count 0 in both communities are included in the random deviates }
  
  
}
\details{
 The following is written from the perspective of using the Poisson lognormal distribution 
 to describe community structure (the distribution of species when sampling individuals 
 from a community of several species).
 
 The following assumes knowledge of the \bold{Details} section of \code{\link{dpoilog}}.
 
 Consider two communities jointly and assume that the log abundances among species have the binormal distribution
 with parameters (\code{mu1},\code{sig1},\code{mu2},\code{sig2},\code{rho}). If sampling intensities are 
 \code{nu1} = \code{nu2} = 1, samples from the communities will have the bivariate Poisson lognormal distribution
 
 \deqn{ P(N_1 = \code{n1}, N_2 = \code{n2}; \;\code{mu1},\code{sig1},\code{mu2},\code{sig2},\code{rho}) = }{
        P(N1 = \code{n1}, N2 = \code{n2}; \code{mu1},\code{sig1},\code{mu2},\code{sig2},\code{rho}) = }
 \deqn{ q(\code{n1},\code{n2};\;\code{mu1},\code{sig1},\code{mu2},\code{sig2},\code{rho}) = }{ 
        q(\code{n1},\code{n2}; \code{mu1},\code{sig1},\code{mu2},\code{sig2},\code{rho}) = }
 \deqn{\int\limits_{-\infty}^{\infty}\int\limits_{-\infty}^{\infty} g_{\code{n1}}(\code{mu1,sig1},u) g_{\code{n2}}(\code{mu2,sig2},v)\,\phi(u,v;\code{rho})\,du dv,}{
       \int_-infty^infty \int_-infty^infty g_n1(\code{mu1},\code{sig1},u) g_n2(\code{mu2},\code{sig2},v) phi(u,v; \code{rho}) du dv,}
 
 where \eqn{\phi(u,v; \code{rho})} here denotes the binormal distribution with zero means, unit variances and correlation \code{rho}.
 In the general case with sampling intensities \code{nu1} and \code{nu2}, \code{mu1} and \code{mu2} should be replaced by 
 \code{mu1} + ln \code{nu1} and \code{mu2} + ln \code{nu2} respectively. In this case, some species will be missing from both samples.
 The number of individuals for observed species will then have the truncated distribution
 
 \deqn{\frac{q(\code{n1},\code{n2};\;\code{mu1},\code{sig1},\code{mu2},\code{sig2},\code{rho})}{1-q(0,0;\;\code{mu1},\code{sig1},\code{mu2},\code{sig2},\code{rho})}}{
       q(\code{n1},\code{n2}; \code{mu1},\code{sig1},\code{mu2},\code{sig2},\code{rho}) / (1 - q(0,0; \code{mu1},\code{sig1},\code{mu2},\code{sig2},\code{rho}))}
       }
\value{
 \command{dbipoilog} returns the density\cr
 \command{rbipoilog} returns random deviates
}
\references{ Engen, S., R. Lande, T. Walla and P. J. DeVries. 2002. Analyzing spatial structure of communities using the two-dimensional 
            Poisson lognormal species abundance model. American Naturalist \bold{160}, 60-73. }
\author{ Vidar Grotan \email{vidar.grotan@ntnu.no} and Steinar Engen   }

\seealso{ \code{\link{bipoilogMLE}} for maximum likelihood estimation }
\examples{

### change in density of n2 for two different values of rho (given n1=10)   
barplot(rbind(dbipoilog(n1=rep(10,21),n2=0:20,mu1=0,mu2=0,sig=1,sig2=1,rho=0.0),
              dbipoilog(n1=rep(10,21),n2=0:20,mu1=0,mu2=0,sig=1,sig2=1,rho=0.8)),
              beside=TRUE,space=c(0,0.2),names.arg=0:20,xlab="n2",col=1:2)
legend(35,0.0012,c("rho=0","rho=0.8"),fill=1:2)


### draw random deviates from a community of 50 species 
rbipoilog(S=50,mu1=0,mu2=0,sig1=1,sig2=2,rho=0.7)

### draw random deviates including zeros
rbipoilog(S=50,mu1=0,mu2=0,sig1=1,sig2=2,rho=0.7,keep0=TRUE)

### draw random deviates with sampling intensities nu1=0.5 and nu2=0.7 
rbipoilog(S=50,mu1=0,mu2=0,sig1=1,sig2=2,rho=0.7,nu1=0.5,nu2=0.7)

### draw random deviates conditioned on a certain number of species 
rbipoilog(S=50,mu1=0,mu2=0,sig1=1,sig2=2,rho=0.7,nu1=0.5,nu2=0.7,condS=TRUE)


### how many species are likely to be observed in at least one of the samples
### (given S,mu1,mu2,sig1,sig2,rho)? 
hist(replicate(1000,nrow(rbipoilog(S=50,mu1=0,mu2=0,sig1=1,sig2=2,rho=0.7))),
     main="", xlab = "Number of species observed in at least one of the samples")

### how many individuals are likely to be observed 
### (given S,mu1,mu2,sig1,sig2,rho)? 
hist(replicate(1000,sum(rbipoilog(S=50,mu1=0,mu2=0,sig1=1,sig2=2,rho=0.7))),
     main="", xlab="sum nr of individuals in both samples")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ distribution }
