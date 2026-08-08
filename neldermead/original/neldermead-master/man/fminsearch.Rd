\name{fminsearch}
\alias{fminsearch}
\title{
  Computation of the unconstrained minimum of given function with the
  Nelder-Mead algorithm.
}
\description{
  This function searches for the unconstrained minimum of a given cost function.
  The provided algorithm is a direct search algorithm, i.e. an algorithm which
  does not use the derivative of the cost function. It is based on the update of
  a simplex, which is a set of k>=n+1 vertices, where each vertex is associated
  with one point and one function value. This algorithm is the Nelder-Mead
  algorithm. This function is based on a specialized use of the more general
  \code{neldermead} function bundle. Users who want to have a more flexible
  solution based on direct search algorithms should consider using the
  \code{neldermead} functions instead of the \code{fminsearch} function.
}
\usage{
  fminsearch(fun = NULL, x0 = NULL, options = NULL, verbose=FALSE)
}
\arguments{
  \item{fun}{A cost function return a numeric scalar.}
  \item{x0}{A numerical vector of initial guesses (length n).}
  \item{options}{A list of optimization options, which drives the behaviour of
    \code{fminsearch}. These options must be set with the \code{optimset}
    function (see \code{?optimset}) which returns a list with the following
    elements: \describe{
      \item{MaxIter}{The maximum number of iterations. The default is 200 * n.}
      \item{MaxFunEvals}{The maximum number of evaluations of the cost function.
        The default is 200 * n.}
      \item{TolFun}{The absolute tolerance on function value. The default value
        is 1.e-4.}
      \item{TolX}{The absolute tolerance on simplex size. The default value is
        1.e-4.}
      \item{Display}{The verbose level.}
      \item{OutputFcn}{The output function, or a list of output functions
        called at the end of each iteration. The default value is NULL.}
      \item{PlotFcns}{The plot function, or a list of plotput functions called
        at the end of each iteration. The default value is empty.}
    }
  }
  \item{verbose}{The verbose option, controlling the amount of messages.}
}
\details{
  \strong{Termination criteria}

  In this section, we describe the termination criteria used by fminsearch.
  The criteria is based on the following variables: \describe{
    \item{ssize}{the current simplex size,}
    \item{shiftfv}{the absolute value of the difference of function value
      between the highest and lowest vertices.}
  }

  If both \code{ssize < options$TolX} and \code{shiftfv < options$TolFun}
  conditions are true, then the iterations stop. The size of the simplex is
  computed using the 'sigmaplus' method of the \pkg{optimsimplex} package. The
  'sigmamplus' size is the maximum length of the vector from each vertex to the
  first vertex. It requires one loop over the vertices of the simplex.

  \strong{The initial simplex}
  
  The \code{fminsearch} algorithm uses a special initial simplex, which is an
  heuristic depending on the initial guess. The strategy chosen by
  \code{fminsearch} corresponds to the content of \code{simplex0method}
  element of the neldermead object (set to 'pfeffer').  It is applied using 
  the content of the \code{simplex0deltausual} (0.05) and
  \code{simplex0deltazero} (0.0075) elements. Pfeffer's method is an heuristic
  which is presented in 'Global Optimization Of Lennard-Jones Atomic Clusters'
  by Ellen Fan. It is due to L. Pfeffer at Stanford. See in the help of
  optimsimplex for more details.

  \strong{The number of iterations}
  
  In this section, we present the default values for the number of iterations in
  \code{fminsearch}.

  The \code{options} input argument is an optional list which can
  contain the \code{MaxIter} field, which stores the maximum number of
  iterations. The default value is 200n, where n is the number of variables.
  The factor 200 has not been chosen by chance, but is the result of experiments
  performed against quadratic functions with increasing space dimension.
  This result is presented in 'Effect of dimensionality on the Nelder-mead
  simplex method' by Lixing Han and Michael Neumann. This paper is based on
  Lixing Han's  PhD, 'Algorithms in Unconstrained Optimization'. The study is
  based on numerical experiments with a quadratic function where the number of
  terms depends on the dimension of the space (i.e. the number of variables).
  Their study showed that the number of iterations required to reach the
  tolerance criteria is roughly 100n. Most iterations are based on inside
  contractions. Since each step of the Nelder-Mead algorithm only require one
  or two function evaluations, the number of required function evaluations in
  this experiment is also roughly 100n.

  \strong{Output and plot functions}
  
  The \code{optimset} function can be used to configure one or more output and
  plot functions.
  The output or plot function is expected to have the following definition:

  \code{myfun <- function(x , optimValues , state)}

  The input arguments \code{x}, \code{optimValues} and \code{state} are
  described in detail in the \code{optimset} help page. The
  \code{optimValues$procedure} field represents the type of step performed at
  the current iteration and can be equal to one of the following strings:
  \itemize{
    \item '' (the empty string),
    \item 'initial simplex',
    \item 'expand',
    \item 'reflect',
    \item 'contract inside',
    \item 'contract outside'.
  }
}
\value{
  Return a object of class neldermead. Use the \code{neldermead.get} to extract
  the following element from the returned object:\describe{
    \item{xopt}{The vector of n numeric values, minimizing the cost function.}
    \item{fopt}{The minimum value of the cost function.}
    \item{exitflag}{The flag associated with exist status of the algorithm.
      The following values are available:\describe{
        \item{-1}{The maximum number of iterations has been reached.}
        \item{0}{The maximum number of function evaluations has been reached.}
        \item{1}{The tolerance on the simplex size and function value delta has
          been reached. This signifies that the algorithm has converged,
          probably to a solution of the problem.}
      }
    }
    \item{output}{A list which stores detailed information about the exit of the
      algorithm. This list contains the following fields:\describe{
        \item{algorithm}{A string containing the definition of the algorithm
          used, i.e. 'Nelder-Mead simplex direct search'.}
        \item{funcCount}{The number of function evaluations.}
        \item{iterations}{The number of iterations.}
        \item{message}{A string containing a termination message.}
      }
    }
  }
}
\author{
Author of Scilab neldermead module: Michael Baudin (INRIA - Digiteo)

Author of R adaptation: Sebastien Bihorel (\email{sb.pmlab@gmail.com})
}
\references{
'Sequential Application of Simplex Designs in Optimisation and Evolutionary
Operation', Spendley, W. and Hext, G. R. and Himsworth, F. R., American
Statistical Association and American Society for Quality, 1962

'A Simplex Method for Function Minimization', Nelder, J. A. and Mead, R., The
Computer Journal, 1965

'Iterative Methods for Optimization', C. T. Kelley, SIAM Frontiers in Applied
Mathematics, 1999

'Algorithm AS47 - Function minimization using a simplex procedure', O'Neill, R.,
Applied Statistics, 1971

'Effect of dimensionality on the nelder-mead simplex method', Lixing Han and
Michael Neumann, Optimization Methods and Software, 21, 1, 1--16, 2006.

'Algorithms in Unconstrained Optimization', Lixing Han, Ph.D., The University of
Connecticut, 2000.

'Global Optimization Of Lennard-Jones Atomic Clusters' Ellen Fan, Thesis,
February 26, 2002, McMaster University
}
\seealso{
  \code{\link[neldermead]{optimset}}
  \code{\link[neldermead]{neldermead}}
}
\examples{
#In the following example, we use the fminsearch function to compute the minimum
#of the Rosenbrock function. We first define the function 'banana', and then use
#the fminsearch function to search the minimum, starting with the initial guess
#(-1.2, 1.0). In this particular case, 85 iterations are performed with 159
#function evaluations
  banana <- function(x){
    y <- 100*(x[2]-x[1]^2)^2 + (1-x[1])^2
  }
  sol <- fminsearch(banana, c(-1.2,1))
  sol

#In the following example, we configure the absolute tolerance on the size of
#the simplex to a larger value, so that the algorithm performs less iterations.
#Since the default value of 'TolX' for the fminsearch function is 1.e-4, we
#decide to use 1.e-2. The optimset function is used to create an optimization
#option list and the field 'TolX' is set to 1.e-2. The options list is then
#passed to the fminsearch function as the third input argument. In this
#particular case, the number of iterations is 70 with 130 function evaluations.

  opt <- optimset(TolX=1.e-2)
  sol <- fminsearch(banana, c(-1.2,1), opt)
  sol
  
#In the following example, we want to produce intermediate outputs of the
#algorithm. We define the outfun function, which takes the current point x as
#input argument. The function plots the current point into the current graphic
#window with the plot function. We use the 'OutputFcn' feature of the optimset
#function and set it to the output function. Then the option list is passed 
#to the fminsearch function. At each iteration, the output function is called 
#back, which creates and update a plot. While this example creates a 2D plot,
#the user may customized the output function so that it writes a message in 
#the console, write some data into a data file, etc... The user can distinguish 
#between the output function (associated with the 'OutputFcn' option) and the
#plot function (associated with the 'PlotFcns' option). See the optimset for
#more details on this feature.

  outfun <- function(x, optimValues, state){
    plot(x[1],x[2],xlim=c(-1.5,1.5),ylim=c(-1.5,1.5))
    par(new=TRUE)
  }
  opt <- optimset(OutputFcn=outfun)
  sol <- fminsearch(banana, c(-1.2,1), opt)
  sol

#The 'Display' option allows to get some input about the intermediate steps of
#the algorithm as well as to be warned in case of a convergence problem.
#In the following example, we present what happens in case of a convergence
#problem. We set the number of iterations to 10, instead of the default 400
#iterations. We know that 85 iterations are required to reach the convergence
#criteria. Therefore, the convergence criteria is not met and the maximum number
#of iterations is reached.

  opt <- optimset(MaxIter=10)
  sol <- fminsearch(banana, c(-1.2,1), opt)

#Since the default value of the 'Display' option is 'notify', a message is
#generated, which warns the user about a possible convergence problem. The
#previous script produces the following output.
# Exiting: Maximum number of iterations has been exceeded
#          - increase MaxIter option.
#          Current function value: 4.1355598

#In the following example, we present how to display intermediate steps used by
#the algorithm. We simply set the 'Display' option to the 'iter' value. This 
#option allows to see the number of function evaluations, the minimum function
#value and which type of simplex step is used for the iteration.
  opt <- optimset(Display='iter')
  sol <- fminsearch(banana, c(-1.2,1), opt)
  sol
}
\keyword{ method }
