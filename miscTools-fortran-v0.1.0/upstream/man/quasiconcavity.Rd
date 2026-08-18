\name{quasiconcavity}
\alias{quasiconcavity}
\alias{quasiconvexity}
\title{Test for quasiconcavity / quasiconvexity}

\description{
   Test wether a function is quasiconcave or quasiconvex.
   The bordered Hessian of this function is checked by
   \code{quasiconcavity}() or \code{quasiconvexity}().
}

\usage{
   quasiconcavity( m, tol = .Machine$double.eps )
   quasiconvexity( m, tol = .Machine$double.eps )
}

\arguments{
   \item{m}{a bordered Hessian matrix or a list containing bordered Hessian matrices}
   \item{tol}{tolerance level (values between \code{-tol} and \code{tol}
      are considered to be zero).}
}

\value{
   locigal or a logical vector (if \code{m} is a list).
}

\references{
   Chiang, A.C. (1984)
   \emph{Fundamental Methods of Mathematical Economics},
   3rd ed., McGraw-Hill.
}

\author{Arne Henningsen}

\examples{
   quasiconcavity( matrix( 0, 3, 3 ) )

   quasiconvexity( matrix( 0, 3, 3 ) )

   m <- list()
   m[[1]] <- matrix( c( 0,-1,-1, -1,-2,3, -1,3,5 ), 3, 3 )
   m[[2]] <- matrix( c( 0,1,-1, 1,-2,3, -1,3,5 ), 3, 3 )

   quasiconcavity( m )

   quasiconvexity( m )
}

\keyword{array}
