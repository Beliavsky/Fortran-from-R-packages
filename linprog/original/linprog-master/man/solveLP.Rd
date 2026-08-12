\name{solveLP}
\alias{solveLP}

\title{Solve Linear Programming / Optimization Problems}

\description{
  Minimizes (or maximizes) \eqn{c'x}, subject to \eqn{A x <= b} and \eqn{x >= 0}.

Note that the inequality signs \code{<=} of the individual linear constraints
in \eqn{A x <= b} can be changed with argument \code{const.dir}.
}

\usage{
solveLP( cvec, bvec, Amat, maximum = FALSE,
   const.dir = rep( "<=", length( bvec ) ),
   maxiter = 1000, zero = 1e-9, tol = 1e-6, dualtol = tol,
   lpSolve = FALSE, solve.dual = FALSE, verbose = 0 )
}

\arguments{
   \item{cvec}{vector \eqn{c} (containing \eqn{n} elements).}
   \item{bvec}{vector \eqn{b} (containing \eqn{m} elements).}
   \item{Amat}{matrix A (of dimension \eqn{m \times n}).}
   \item{maximum}{logical. Should we maximize or minimize (the default)?}
   \item{const.dir}{vector of character strings giving the directions
      of the constraints: each value should be one of "<," "<=," "=," "==,"
      ">," or ">=". (In each pair the two values are identical.)}
   \item{maxiter}{maximum number of iterations.}
   \item{zero}{numbers smaller than this value (in absolute terms) are
      set to zero.}
   \item{tol}{if the constraints are violated by more than this number,
      the returned component \code{status} is set to \code{3}.}
   \item{dualtol}{if the constraints in the dual problem are violated by more
      than this number, the returned status is non-zero.}
   \item{lpSolve}{logical. Should the package 'lpSolve' be used to solve
      the LP problem?}
   \item{solve.dual}{logical value indicating if the dual problem should
      also be solved.}
   \item{verbose}{an optional integer variable to indicate how many
      intermediate results should be printed
      (0 = no output; 4 = maximum output).}
}

\details{
   This function uses the Simplex algorithm of George B. Dantzig (1947)
   and provides detailed results (e.g. dual prices, sensitivity analysis
   and stability analysis).\cr
   If the solution \eqn{x=0} is not feasible, a 2-phase procedure is
   applied.\cr
   Values of the simplex tableau that are actually zero might get small
   (positive or negative) numbers due to rounding errors, which might
   lead to artificial restrictions. Therefore, all values that are smaller
   (in absolute terms) than the value of \code{zero} (default is 1e-10) are
   set to 0.\cr
   Solving the Linear Programming problem by the package \code{lpSolve}
   (of course) requires the installation of this package, which is available
   on CRAN (\url{https://cran.r-project.org/package=lpSolve}).
   Since the \code{lpSolve} package uses C-code and this (\code{linprog})
   package is not optimized for speed, the former is much faster.
   However, this package provides more detailed results (e.g. dual values,
   stability and sensitivity analysis).\cr
   This function has not been tested extensively and might not solve all
   feasible problems (or might even lead to wrong results). However, you can
   export your LP to a standard MPS file via \code{\link{writeMps}} and check
   it with other software (e.g. \code{lp_solve}, see
   \url{https://lpsolve.sourceforge.net/5.5/}).\cr
   Equality constraints are not implemented yet.
}

\value{
  \code{solveLP} returns a list of the class \code{solveLP}
  containing following objects:

  \item{opt}{optimal value (minimum or maximum) of the objective function.}
  \item{solution}{vector of optimal values of the variables.}
  \item{iter1}{iterations of Simplex algorithm in phase 1.}
  \item{iter2}{iterations of Simplex algorithm in phase 2.}
  \item{basvar}{vector of basic (=non-zero) variables (at optimum).}
  \item{con}{matrix of results regarding the constraints:\cr
            1st column = maximum values (=vector \eqn{b});\cr
            2nd column = actual values;\cr
            3rd column = differences between maximum and actual values;\cr
            4th column = dual prices (shadow prices);\cr
            5th column = valid region for dual prices.}
  \item{allvar}{matrix of results regarding all variables (including slack variables):\cr
            1st column = optimal values;\cr
            2nd column = values of vector \eqn{c};\cr
            3rd column = minimum of vector \eqn{c} that does \emph{not} change the solution;\cr
            4th column = maximum of vector \eqn{c} that does \emph{not} change the solution;\cr
            5th column = derivatives to the objective function;\cr
            6th column = valid region for these derivatives.}
   \item{status}{numeric. Indicates if the optimization did succeed:\cr
      0 = success; 1 = lpSolve did not succeed;
      2 = solving the dual problem did not succeed;
      3 = constraints are violated at the solution
         (internal error or large rounding errors);
      4 = simplex algorithm phase 1 did not find a solution within
         the number of iterations specified by argument \code{maxiter};
      5 = simplex algorithm phase 2 did not find the optimal solution within
         the number of iterations specified by argument \code{maxiter}.}
   \item{lpStatus}{numeric. Return code of \code{\link[lpSolve]{lp}}
      (only if argument \code{lpSolve} is \code{TRUE}).}
   \item{dualStatus}{numeric. Return code from solving the dual problem
      (only if argument \code{solve.dual} is \code{TRUE}).}
   \item{maximum}{logical. Indicates whether the objective function
            was maximized or minimized.}
   \item{Tab}{final 'Tableau' of the Simplex algorith.}
   \item{lpSolve}{logical. Has the package 'lpSolve' been used to solve
      the LP problem.}
   \item{solve.dual}{logical. Argument \code{solve.dual}.}
   \item{maxiter}{numeric. Argument \code{maxiter}.}
}

\references{

  Dantzig, George B. (1951),
  \emph{Maximization of a linear function of variables subject to linear inequalities},
  in Koopmans, T.C. (ed.), Activity analysis of production and allocation,
  John Wiley & Sons, New York, p. 339-347.

  Steinhauser, Hugo; Cay Langbehn and Uwe Peters (1992),
  Einfuehrung in die landwirtschaftliche Betriebslehre. Allgemeiner Teil,
  5th ed., Ulmer, Stuttgart.

  Witte, Thomas; Joerg-Frieder Deppe and Axel Born (1975),
  Lineare Programmierung. Einfuehrung fuer Wirtschaftswissenschaftler,
  Gabler-Verlag, Wiesbaden.

}

\author{
   Arne Henningsen
}

\seealso{
   \code{\link{readMps}} and \code{\link{writeMps}}
}

\examples{

## example of Steinhauser, Langbehn and Peters (1992)
## Production activities
cvec <- c(1800, 600, 600)  # gross margins
names(cvec) <- c("Cows","Bulls","Pigs")

## Constraints (quasi-fix factors)
bvec <- c(40, 90, 2500)  # endowment
names(bvec) <- c("Land","Stable","Labor")

## Needs of Production activities
Amat <- rbind( c(  0.7,   0.35,   0 ),
               c(  1.5,   1,      3 ),
               c( 50,    12.5,   20 ) )

## Maximize the gross margin
solveLP( cvec, bvec, Amat, TRUE )


## example 1.1.3 of Witte, Deppe and Born (1975)
## Two types of Feed
cvec <- c(2.5, 2 )  # prices of feed
names(cvec) <- c("Feed1","Feed2")

## Constraints (minimum (<0) and maximum (>0) contents)
bvec <- c(-10, -1.5, 12)
names(bvec) <- c("Protein","Fat","Fibre")

## Matrix A
Amat <- rbind( c( -1.6,  -2.4 ),
               c( -0.5,  -0.2 ),
               c(  2.0,   2.0 ) )

## Minimize the cost
solveLP( cvec, bvec, Amat )

# the same optimisation using argument const.dir
solveLP( cvec, abs( bvec ), abs( Amat ), const.dir = c( ">=", ">=", "<=" ) )


## There are also several other ways to put the data into the arrays, e.g.:
bvec <- c( Protein = -10.0,
           Fat     =  -1.5,
           Fibre   =  12.0 )
cvec <- c( Feed1 = 2.5,
           Feed2 = 2.0 )
Amat <- matrix( 0, length(bvec), length(cvec) )
rownames(Amat) <- names(bvec)
colnames(Amat) <- names(cvec)
Amat[ "Protein", "Feed1" ] <-  -1.6
Amat[ "Fat",     "Feed1" ] <-  -0.5
Amat[ "Fibre",   "Feed1" ] <-   2.0
Amat[ "Protein", "Feed2" ] <-  -2.4
Amat[ "Fat",     "Feed2" ] <-  -0.2
Amat[ "Fibre",   "Feed2" ] <-   2.0
solveLP( cvec, bvec, Amat )
}

\keyword{optimize}
