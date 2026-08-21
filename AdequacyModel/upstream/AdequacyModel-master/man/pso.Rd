\name{pso}
\alias{pso}
\encoding{UTF-8}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Adequacy of models
}
\description{
In computer science, the PSO is a computational method for optimization of parametric and multiparametric functions. The PSO algorithm is a meta-heuristic method, which has been providing good solutions for problems of global optimization functions with box-constrained. As in most heuristic methods that are inspired by biological phenomena, the PSO method is inspired by the behavior of flying birds. The philosophical idea of the PSO algorithm is based on the collective behavior of birds (particle) in search of food (point of global optimal).

The \code{pso} function is an efficient function for global minimization, wherein it is not necessary to provide Initial kicks. This is the function for general purpose optimization. 
}
\usage{
pso(func, S = 350, lim_inf, lim_sup, e = 0.0001, data = NULL, N = 500, prop = 0.2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{func}{
        Objective function, i.e, function to be minimized;
}
  \item{S}{
        Particle number considered. By default, S = 350;
}
  \item{lim_inf}{
        Vector with the lower limit of the search for the parameters of the objective function that will be minimized;
}
  \item{lim_sup}{
        Vector with the upper limits of search for the parameters of the objective function that will be minimized;
}
  \item{e}{
        Stop value of the algorithm, i.e., if the variance of the last 20 minimum values is less than or equal to \code{e}, the algorithm will converge to the global minimum. By default, \code{e = 0.0001};
}
  \item{data}{
        Vector of data provided in the event of function to be minimized (passed as an argument for \code{func}) involve some data set. An example of a function that you should inform a data set is when we want to minimize the log-likelihood function multiplied by -1 (-log-likelihood). By defatul, \code{data = NULL};
}
  \item{N}{
        Minimum number of iterations. By default, \code{N = 500};   
}

  \item{prop}{
        Proportion of last minimum value that is calculated variance used as a stopping criterion. That is, if the number of iterations is greater or equal to the minimum number of iterations \code{N}, calculate the variance of the last values of minimum obtained, wherein 0 <= prop <= 1}.   


}
\details{
The PSO optimizes a problem by having a population of candidate solutions and moving these particles around in the search-space according to simple mathematical formulae over the particle’s position and velocity. The movement of the particles in the search space is randomized. Each iteration of the PSO algorithm, there is a leader particle, which is the particle that minimizes the objective function in the corresponding iteration. The remaining particles arranged in the search region will follow the leader particle randomly and sweep the area around this leading particle. In this local search process, another particle may become the new leader particle and the other particles will follow the new leader randomly. Each particle
arranged in the search region has a velocity vector and position vector and its movement in the search region is given by changes in these vectors.

As a stopping criterion is considered the variance of the last 20 minimum values estimated by the algorithm. If this variance is less or equal the \code{e} the algorithm will stop providing the global minimum value. This is a conditional criterion, which will only be evaluated if the number of iterations is greater than or equal to the minimum number of iterations set to  \code{N}.

The amount of minimum values considered in the calculation of the variance is given by the proportion of minimum values established by the argument \code{prop} which by default is \code{prop = 0.2}. That is, if the last 20\% (\code{prop = 0.2}) of the minimum values has less variance than or equal to \code{e}, the algorithm will stop global search, indicating convergence according to the established criteria. This indicates that there was no significant improvements in this proportion of last iterations. 
}



\value{
\item{par_pso}{Global minimum point;}
\item{f_pso}{Global minimum value.}
}

%\seealso{For details about the optimization methodologies may view the functions \code{\link{optim}}, \code{\link{ks.test}}, %\code{\link{nlminb}}.}

\references{
Beni G, Wang J (1993). Swarm intelligence in cellular robotic systems. pp. 703-712.

Eberhart RC, Kennedy J (1995). A new optimizer using particle swarm theory. In Proceedings of the sixth international symposium on micro machine and human science, volume 1, pp. 39-43. New York, NY.

Kennedy J, Kennedy JF, Eberhart RC, Shi Y (2001). Swarm intelligence. Morgan Kaufmann.

Shi Y, Eberhart R (1998). A modified particle swarm optimizer. In Evolutionary Computation Proceedings, 1998. IEEE World Congress on Computational Intelligence., The 1998 IEEE International Conference on, pp. 69-73. IEEE.
}

\note{
In other versions of the package, the paper with more details that complement the documentation of this function will be provided in the above references and this note will be undone.
}

\examples{

# The objective functions below are rather difficult to be optimized. 
# However, the function pso has great results.

# Example 1 (Easom function): 

easom_function <- function(par,x){
  x1 = par[1]
  x2 = par[2]
  -cos(x1)*cos(x2)*exp(-((x1-pi)^2 + (x2-pi)^2))
}

set.seed(0)
result_1 = pso(func = easom_function, S = 500, lim_inf = c(-10,-10), lim_sup = c(10,10),
               e = 0.00001)
result_1$par

# Example 2 (Holder table function):

holder <- function(par,x){
  x1 = par[1]
  x2 = par[2]
  -abs(sin(x1)*cos(x2) * exp(abs(1 - sqrt(x1^2+x2^2)/pi)))
}

set.seed(0)
result_2 = pso(func = holder, S = 700, lim_inf = c(-10,-10), lim_sup = c(10,10),
               e = 0.00001, N=500)
result_2$par

# Example 3:

f_pso <- function(par,x){
  theta = par[1]
  -(6 + theta^2 * sin(14*theta))
}

set.seed(0)
result_3 <- pso(func = f_pso, S = 500, lim_inf = c(-2.5), lim_sup = c(2.5), e = 0.0001)
result_3$par

# TO RUN THE CODE BELOW, UNCOMMENT THE CODES.

# Example 4 (maximizing a function of the log-likelihood function):

# pdf_exp <- function(par,x){
#  lambda = par[1]
#  lambda*exp(-lambda*x)
#}

# -log-likelihood function of the exponential distribution.
#likelihood <- function(par,x){
#  lambda = par[1]
#  -sum(log(pdf_exp(par,x)))
#}

#set.seed(0)
#random_data1 = rexp(500,1)
#result_1 = pso(func = likelihood, S = 250, lim_inf = c(0), lim_sup = c(100), e = 0.0001,
#    data = random_data1, N = 50, prop = 0.2)

#x = seq(0,ceiling(max(random_data1)), length.out = 500)
#hist(random_data1, probability = TRUE)
#lines(x, pdf_exp(par = result_1$par, x), col = "blue")

# Example 5 (maximizing a function of the log-likelihood function):

# Probability density function (Weibull) 
#pdf_weibull <- function(par,x){
#  a = par[1]
#  b = par[2]
#  dweibull(x,shape=a,scale=b)
#}

# -log-likelihood function of the Weibull distribution.
#likelihood <- function(par,x){
#  -sum(log(pdf_weibull(par,x)))
#}

#set.seed(0)
#random_data2 = rweibull(250,2,2)
#result_2 = pso(func = likelihood, S = 250, lim_inf = c(0,0), lim_sup = c(10,10), e = 0.0001,
#               data = random_data2, N = 50, prop = 0.2)
    
#x = seq(0,ceiling(max(random_data2)), length.out = 500)
#hist(random_data2, probability = TRUE, ylim = c(0,0.5))
#lines(x, pdf_weibull(par = result_2$par, x), col = "blue")
    
}

\author{
Pedro Rafael Diniz Marinho \email{pedro.rafael.marinho@gmail.com}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{AdequacyModel}
\keyword{optimization}
\keyword{PSO}
\keyword{pso}
\keyword{distribution}
\keyword{survival}