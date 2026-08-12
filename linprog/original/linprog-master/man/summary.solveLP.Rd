\name{summary.solveLP}
\alias{summary.solveLP}
\alias{print.summary.solveLP}
\title{Summary Results for Objects of Class solveLP}

\description{
  These methods prepare and print summary results of the Linear Programming
  algorithm.
}

\usage{
   \method{summary}{solveLP}(object,...)
   \method{print}{summary.solveLP}(x,...)
}

\arguments{
  \item{object}{an object returned by \code{\link{solveLP}}.}
  \item{x}{an object returned by \code{summary.solveLP}.}
  \item{...}{currently ignored.}
}

\value{
  \code{summary.solveLP} returns an object of class \code{summary.solveLP}.
  \code{print.summary.solveLP} invisibly returns the object given
  in argument \code{x}.
}

\author{
   Arne Henningsen
}

\seealso{
   \code{\link{solveLP}}, \code{\link{print.solveLP}},
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

## prepare and print the summary results
summary( res )
}

\keyword{ optimize }



