    \name{romano_ex}
     \docType{data}
     \alias{romano_ex}
     \title{Simulated values of a Romano & Siegel example}
     \description{Data frame containing 100 values of x,y,z generated as follows: x0 = rnorm(100); y = rnorm(100); z = rnorm(100); x = abs(x0)*sign(y*z). All pairs are independent but the three series are not.}
     \usage{data(romano_ex)}
     \format{dataframe}
     \keyword{datasets}

\examples{
data(romano_ex)
plot(romano_ex$x)
}
\keyword{datasets}

