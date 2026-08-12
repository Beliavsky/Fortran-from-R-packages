\name{print.solveLP}
\alias{print.solveLP}
\title{Print Objects of Class solveLP}

\description{
  This method prints the results of the Linear Programming algorithm.
}

\usage{
   \method{print}{solveLP}( x, digits=6, ...)
}

\arguments{
  \item{x}{an object returned by \code{\link{solveLP}}.}
  \item{digits}{number of digits to print.}
  \item{...}{currently ignored.}
}

\value{
  \code{print.solveLP} invisibly returns the object given in argument \code{x}.
}

\author{
   Arne Henningsen
}

\seealso{
   \code{\link{solveLP}}, \code{\link{summary.solveLP}},
   \code{\link{readMps}}, \code{\link{writeMps}}
}


\examples{

## example of Steinhauser, Langbehn and Peters (1992)
\dontrun{library( linprog )}

## Production activities
cvec <- c(1800, 600, 600)  # gross margins
names(cvec) <- c("Milk","Bulls","Pigs")

## Constraints (quasi-fix factors)
bvec <- c(40, 90, 2500)  # endowment
names(bvec) <- c("Land","Stable","Labor")

## Needs of Production activities
Amat <- rbind( c(  0.7,   0.35,   0 ),
               c(  1.5,   1,      3 ),
               c( 50,    12.5,   20 ) )

## Maximize the gross margin
res <- solveLP( cvec, bvec, Amat, TRUE )

## print the results
print( res )
}

\keyword{ optimize }
