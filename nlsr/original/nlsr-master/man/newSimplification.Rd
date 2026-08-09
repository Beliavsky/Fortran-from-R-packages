\name{newSimplification}
\alias{newSimplification}
\title{newSimplification}
\usage{
newSimplification(expr, test, simplification, do_eval = FALSE, 
                       simpEnv = sysSimplifications)
}
\arguments{
\item{expr}{An expression to simplify.}

\item{test}{An expression to evaluate whether the simplification
should be applied to \code{expr}.}

\item{simplification}{An equivalent but simpler version of the expression.}

\item{do_eval}{Whether the simplification should be evaluated
or stored as-is.}

\item{simpEnv}{The environment in which the simplification is stored.}
}
\description{
Define a new simplification expression
}
\value{
If \code{expr} is missing, list all the functions with known 
simplifications.

If \code{test} is missing, list all the known simplifications for
\code{expr}.

Otherwise, add the new simplification to the list of possible simplifications.
}
\examples{
# The unary + operation can always be safely removed:
newSimplification(+a, TRUE, a)

# The unary - operation can be absorbed into numeric values:
newSimplification(-a, is.numeric(a), -a, do_eval = TRUE)

# Adding zero to anything can be skipped:
newSimplification(a + b, isZERO(b), a)
}
\keyword{internal}