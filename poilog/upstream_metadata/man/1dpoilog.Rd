\name{Poisson lognormal}
\alias{dpoilog}
\alias{rpoilog}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Poisson lognormal distribution }
\description{
  Density and random generation for the Poisson lognormal distribution with parameters \code{mu} and \code{sig}.
}
\usage{
dpoilog(n, mu, sig)
rpoilog(S, mu, sig, nu=1, condS=FALSE, keep0=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{ vector of observed individuals for each species }
  \item{S}{ number of species in the community }
  \item{mu}{ mean of lognormal distribution  }
  \item{sig}{ standard deviation of lognormal distribution }
  \item{nu}{ sampling intensity, defaults to 1 }
  \item{condS}{ logical; if TRUE random deviates are conditional on S}
  \item{keep0}{ logical; if TRUE species with count 0 are included in the random deviates }
  
}
\details{

  The following is written from the perspective of using the Poisson lognormal distribution 
  to describe community structure (the distribution of species when sampling individuals 
  from a community of several species). 
  
  Under the assumption of random sampling, the number of individuals sampled from a given 
  species with abundance \emph{y}, say \emph{N}, is Poisson distributed with mean \eqn{\code{nu}\,y}{\code{nu}\emph{y}} 
  where the parameter \code{nu} expresses the sampling intensity. If ln \emph{y} is normally distributed 
  with mean \code{mu} and standard deviation \code{sig} among species, then the vector of individuals sampled 
  from all \code{S} species then constitutes a sample from the Poisson lognormal distribution
  with parameters (\code{mu} + ln \code{nu}, \code{sig}), where \code{mu} and \code{sig} 
  are the mean and standard deviation of the log abundances. For \code{nu} = 1, this is the Poisson 
  lognormal distribution with parameters (\code{mu},\code{sig}) which may be written in the form
    
  \deqn{P(N=\code{n};\code{mu},\code{sig}) = q(\code{n};\code{mu},\code{sig}) = \int\limits_{-\infty}^{\infty} g_\code{n}(\code{mu},\code{sig},u)\phi(u)\;du,}{
        P(N=\code{n}; \code{mu},\code{sig}) = q(\code{n}; \code{mu},\code{sig}) = \int_-infty^infty g_\code{n}(\code{mu},\code{sig},u) phi(u) du,}
  
  where \eqn{\phi(u)} is the standard normal distribution and
    
  \deqn{g_\code{n}(\code{mu},\code{sig},u) = \frac{\exp(u\,\code{sig}\,\code{n} + \code{mu}\,\code{n} - \exp(u\,\code{sig} + \code{mu}))}{\code{n}!}}{
        g_n(\code{mu},\code{sig},u) = exp(\code{mu} \code{sig} \code{n} + \code{mu} \code{n} + exp(-\code{mu} \code{sig} + \code{mu})) / \code{n}!} 
       
   
  Since \code{S} is usually unknown, we only consider the observed number of individuals for the observed species.
  With a general sampling intensity \code{nu}, the distribution of the number of individuals then follows the
  zero-truncated Poisson lognormal distribution
    
   \deqn{\frac{q(\code{n};\code{mu},\code{sig})}{1-q(0;\code{mu},\code{sig})}}{
          q(\code{n}; \code{mu},\code{sig})/(1 - q(0; \code{mu},\code{sig}))}
}
\value{
  \command{dpoilog} returns the density\cr    
  \command{rpoilog} returns random deviates
}
\references{ Engen, S., R. Lande, T. Walla & P. J. DeVries. 2002. Analyzing spatial structure of communities using the two-dimensional 
            Poisson lognormal species abundance model. American Naturalist 160: 60-73.
 }
\author{ Vidar Grotan \email{vidar.grotan@ntnu.no} and Steinar Engen   }

\seealso{\code{\link{poilogMLE}} for ML estimation}
\examples{

### plot density for given parameters 
barplot(dpoilog(n=0:20,mu=2,sig=1),names.arg=0:20)

### draw random deviates from a community of 50 species 
rpoilog(S=50,mu=2,sig=1)

### draw random deviates including zeros 
rpoilog(S=50,mu=2,sig=1,keep0=TRUE)

### draw random deviates with sampling intensity = 0.5 
rpoilog(S=50,mu=2,sig=1,nu=0.5)

### how many species are likely to be observed 
### (given S,mu,sig2 and nu)? 
hist(replicate(1000,length(rpoilog(S=30,mu=0,sig=3,nu=0.7))))

### how many individuals are likely to be observed
### (given S,mu,sig2 and nu)? 
hist(replicate(1000,sum(rpoilog(S=30,mu=0,sig=3,nu=0.7))))


}
\keyword{ distribution }
