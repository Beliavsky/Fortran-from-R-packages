\name{psoptim}
\alias{psoptim}

\title{Particle Swarm OPTIMization}

\description{
 Particle swarm optimization. The maximum is searched.
}	
	
\usage{
psoptim(FUN, n=100, max.loop=100, w=0.9, c1=0.2, c2=0.2,
        xmin, xmax, vmax=c(4,4), seed=10, anim=TRUE)
}

\arguments{
	  \item{FUN}{
	the optimized function with a vector as parameter
	}
	  \item{n}{
	number of particles
	}
	  \item{max.loop}{
	maximal number of iterations
	}
	  \item{w}{
	inertia weight
	}
	  \item{c1}{
	coefficient of the self-recognition component
	}
	  \item{c2}{
	coefficient of the social component
	}
	  \item{xmin}{
	vector of position constraints - minimal values
	}
	  \item{xmax}{
	vector of position constraints - maximal values
	}
	\item{vmax}{
	vector of velocity constraints in each direction
	}
	\item{seed}{
	seed for random values
	}
	\item{anim}{
	logical; if \code{TRUE} (dafault), animation of the optimization process is shown
	}
}

\details{
The \code{i}-th particle velocity \code{v} in \code{j}-th direction is calculated in \code{t} iteration according to:

\code{
v[ij](t+1) = w*v[ij](t) + c1*r1*(xP[ij](t) - x[ij](t)) + c2*r2*(xS[j](t) - x[ij](t)).
}

where: \code{r1} and \code{r2} are random values, \code{w} is inertia weight, \code{c1} is a coefficient of the self-recognition component and \code{c2} is a coefficient of the social component. \code{xP} denotes so far best position of the particle and \code{xS} - the best position among the swarm.

The new position (coordinates) is calculated as:

\code{
    x[ij](t+1) = x[ij](t) + v[ij](t+1).
}
In the current version of the package, the function works without checking the correctness of the given arguments.
}
\value{
A list with the two components:

\code{sol}	solution, i.e. the best set of parameters found.

\code{val}	the best fitness function found.
}
\references{
Abraham A, Guo H, Liu H. (2006) \emph{Swarm Intelligence: Foundations, Perspectives and Applications} in 
Nedjah N, Mourelle L. (eds.): "Swarm Intelligent Systems", Springer, Berlin Heidelberg, pp. 3-25.

Banks A, Vincent J, Anyakoha C. (2007) \emph{A review of particle swarm optimization. Part I: background and development}. Natural Computing, vol. 6, No. 4, pp. 467-484.

Dorigo M, Stutzle T. (2004) \emph{Ant Colony Optimization}, MIT Press.

Eberhart R, Yuhui S. (2001) \emph{Particle swarm optimization: developments, applications and resources}, Congress on Evolutionary Computation. Seoul, 
Korea.
}



\author{
Krzysztof Ciupke, \code{<krzysztof.ciupke at polsl.pl>}
}

\examples{
n <- 50
m.l <- 50
w <- 0.95
c1 <- 0.2
c2 <- 0.2
xmin <- c(-5.12, -5.12)
xmax <- c(5.12, 5.12)
vmax <- c(4, 4)

g <- function(x){  
  -(20 + x[,1]^2 + x[,2]^2 - 10*(cos(2*pi*x[,1]) + cos(2*pi*x[,2])))
}

psoptim(FUN=g, n=n, max.loop=m.l, w=w, c1=c1, c2=c2,
        xmin=xmin, xmax=xmax, vmax=vmax, seed=5, anim=FALSE)
}
