\name{dex}
\alias{dex}
\title{dex}
\usage{
dex(x, do_substitute = NA, verbose = FALSE)
}
\arguments{
\item{x}{An expression represented in a variety of ways.
See Details.}

\item{do_substitute}{Whether to use the expression passed as \code{x}, or
to evaluate it and use its value.}

\item{verbose}{Print messages describing the process.}
}
\description{
Calculate expression for derivative calculations.
Mainly for internal use, \code{dex()} converts 
its input to an expression suitable for use
in \code{\link{nlsDeriv}} and related functions.
}
\section{Details}{

   If \code{do_substitute} is \code{NA}, the following 
   rules are used:
     
     \itemize{
        \item An attempt is made to evaluate \code{x}.  If that fails,
         the expression is used.
     
        \item If the evaluation succeeds and the value is a character
          vector, it is parsed.
     
        \item If the value is not a character vector and the
         expression is a single name, the value is used.
     }
     Otherwise, the expression is used.
     
     Once the expression is determined it may be simplified,
     by extracting the language object from a length-one
     expression vector, or the right-hand-side from a 
     formula.
     
     Normally a warning will be issued if \code{x} is a formula
     containing a left-hand side.  To suppress this, 
     wrap the formula in \code{expression()}, or pass it 
     as a character string to be parsed.
}

\section{value}{

     An expression or language object suitable as input
     to \code{\link{nlsDeriv}} and related functions.
}

\examples{
   aa <- dex(~ x^2)
   aa
   str(aa)
   bb <- dex(expression(x^2))
   bb
   str(bb)
   cc <- dex("x^2")
   cc
   str(cc)
   
}
\author{
Duncan Murdoch
}
\keyword{internal}
