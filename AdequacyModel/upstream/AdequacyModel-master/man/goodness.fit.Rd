\name{goodness.fit}
\alias{goodness.fit}
\encoding{UTF-8}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Adequacy of models
}
\description{
This function provides some useful statistics to assess the quality of fit of probabilistic models, including the statistics Cramér-von Mises and Anderson-Darling. These statistics are often used to compare models not fitted. You can also calculate other goodness of fit such as AIC, CAIC, BIC, HQIC and Kolmogorov-Smirnov test.
}
\usage{
goodness.fit(pdf, cdf, starts, data, method = "PSO", domain = c(0,Inf),
             mle = NULL,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pdf}{
        Probability density function;
}
  \item{cdf}{
        Cumulative distribution function;
}
  \item{starts}{
        Initial parameters to maximize the likelihood function;
}
  \item{data}{
        Data vector;
}
  \item{method}{
        Method used for minimization of the function \code{-log(likelihood)}. The methods supported are: \code{PSO} (default), \code{BFGS}, \code{Nelder-Mead}, \code{SANN}, \code{CG}. Can also be transmitted only the first letter of the methodology, i.e., \code{P}, \code{B}, \code{N}, \code{S} or \code{C} respectively;
}
  \item{domain}{
        Domain of probability density function. By default the domain of probability density function is the open interval 0 to infinity.This option must be an vector with two values;  
}
  \item{mle}{
        Vector with the estimation maximum likelihood. This option should be used if you already have knowledge of the maximum likelihood estimates. The default is \code{NULL}, ie, the function will try to obtain the estimates of maximum likelihoods;   
}
  \item{...}{If \code{method = "PSO"} or \code{method = "P"}, inform the arguments of the \code{\link{pso}} function. Get details about these arguments into \code{\link{pso}}. Basically the arguments that should be provided are the vectors \code{lim_inf} and \code{lim_sup}. The other parameters of the \code{\link{pso}} function can be informed in the desired configuration. However, may be omitted if the default configuration is sufficient.}
}
\details{
The function \code{goodness.fit} returns statistics KS (Kolmogorov-Smirnov), A (Anderson-Darling), W (Cramér-von Misses). Are also calculated other measures of goodness of fit. These functions are: AIC (Akaike Information Criterion), CAIC (Consistent Akaikes Information Criterion), BIC (Bayesian Information Criterion) and HQIC (Hannan-Quinn information criterion).

The Kolmogorov-Smirnov test may return \code{NA} with a certain frequency. The return \code{NA} informs that the statistical \code{KS} is not reliable for the data set used. More details about this issue can be obtained from \code{\link{ks.test}}.

By default, the function calculates the maximum likelihood estimates. The errors of the estimates are also calculated. In cases that the function can not obtain the maximum likelihood estimates, the change of the values initial, in some cases, resolve the problem. You can also enter with the maximum likelihood estimation if there is already prior knowledge. 
}

\value{
\item{W}{Statistic Cramér-von Misses;}
\item{A}{Statistic Anderson Darling;}
\item{KS}{Kolmogorov Smirnov test;}
\item{mle}{Maximum likelihood estimates;}
\item{AIC}{Akaike Information Criterion;}
\item{CAIC}{Consistent Akaikes Information Criterion;}
\item{BIC}{Bayesian Information Criterion;}
\item{HQIC}{Hannan-Quinn information criterion;}
\item{Erro}{Standard errors of the maximum likelihood estimates;}
\item{Value}{Minimum value of the function -log(likelihood);}
\item{Convergence}{0 indicates successful completion and 1 indicates that the iteration limit maxit had been reached. More details at \code{\link{optim}}.}
}

\seealso{For details about the optimization methodologies may view the functions \code{\link{pso}} and \code{\link{optim}}.}

\references{
Chen, G., Balakrishnan, N. (1995). A general purpose approximate goodness-of-fit test. Journal of Quality Technology, 27, 154-161.

Hannan, E. J. and Quinn, B. G. (1979). The Determination of the Order of an Autoregression. Journal of the Royal Statistical Society, Series B, 41, 190-195.

Nocedal, J. and Wright, S. J. (1999) Numerical Optimization. Springer.

Sakamoto, Y., Ishiguro, M. and Kitagawa G. (1986). Akaike Information Criterion Statistics. D. Reidel Publishing Company.
}

\note{
It is not necessary to define the likelihood function or log-likelihood. You only need to define the probability density function and distribution function. 
}

\examples{

# Example 1:

data(carbone)

# Exponentiated Weibull - Probability density function.
pdf_expweibull <- function(par,x){
  beta = par[1]
  c = par[2]
  a = par[3]
  a * beta * c * exp(-(beta*x)^c) * (beta*x)^(c-1) * (1 - exp(-(beta*x)^c))^(a-1)
}

# Exponentiated Weibull - Cumulative distribution function.
cdf_expweibull <- function(par,x){
  beta = par[1]
  c = par[2]
  a = par[3]
  (1 - exp(-(beta*x)^c))^a
}

set.seed(0)
result_1 = goodness.fit(pdf = pdf_expweibull, cdf = cdf_expweibull, 
                        starts = c(1,1,1), data = carbone, method = "PSO",
                        domain = c(0,Inf),mle = NULL, lim_inf = c(0,0,0),
                        lim_sup = c(2,2,2), S = 250, prop=0.1, N=50)
             
x = seq(0, 6, length.out = 500)
hist(carbone, probability = TRUE)
lines(x, pdf_expweibull(x, par = result_1$mle), col = "blue")

# Example 2:

pdf_weibull <- function(par,x){
  a = par[1]
  b = par[2]
  dweibull(x, shape = a, scale = b)
}

cdf_weibull <- function(par,x){
  a = par[1]
  b = par[2]
  pweibull(x, shape = a, scale = b)
}

set.seed(0)
random_data2 = rweibull(250,2,2)
result_2 = goodness.fit(pdf = pdf_weibull, cdf = cdf_weibull, starts = c(1,1), data = random_data2,
             method = "PSO", domain = c(0,Inf), mle = NULL, lim_inf = c(0,0), lim_sup = c(10,10),
             N = 100, S = 250)

x = seq(0,ceiling(max(random_data2)), length.out = 500)
hist(random_data2, probability = TRUE)
lines(x, pdf_weibull(par = result_2$mle, x), col = "blue")

# TO RUN THE CODE BELOW, UNCOMMENT THE CODES.

# Example 3:

# Kumaraswamy Beta - Probability density function.
#pdf_kwbeta <- function(par,x){
#  beta = par[1]
#  a = par[2]
#  alpha = par[3]
#  b = par[4]
#  (a*b*x^(alpha-1)*(1-x)^(beta-1)*(pbeta(x,alpha,beta))^(a-1)*
#  (1-pbeta(x,alpha,beta)^a)^(b-1))/beta(alpha,beta) 
#}

# Kumaraswamy Beta - Cumulative distribution function.
#cdf_kwbeta <- function(par,x){
#  beta = par[1]
#  a = par[2]
#  alpha = par[3]
#  b = par[4]
#  1 - (1 - pbeta(x,alpha,beta)^a)^b
#}

#set.seed(0)
#random_data3 = rbeta(150,2,2.2)

#system.time(goodness.fit(pdf = pdf_kwbeta, cdf = cdf_kwbeta, starts = c(1,1,1,1),
#              data = random_data3, method = "PSO", domain = c(0,1), lim_inf = c(0,0,0,0),
#              lim_sup = c(10,10,10,10), S = 200, prop = 0.1, N = 40))


}

\author{
Pedro Rafael Diniz Marinho \email{pedro.rafael.marinho@gmail.com}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{AdequacyModel}
\keyword{pso}
\keyword{AIC}
\keyword{BIC}
\keyword{CAIC}
\keyword{distribution}
\keyword{survival}
