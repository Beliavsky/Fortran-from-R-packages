\name{dfoptim}
\alias{dfoptim-package}
\docType{package}
\title{
Derivative-Free Optimization
}
\description{
Derivative-Free optimization algorithms. These algorithms do not require gradient information.  
More importantly, they can be used to solve non-smooth optimization problems. They can also handle box constraints on parameters.
}
\details{
\tabular{ll}{
Package: \tab dfoptim\cr
Type: \tab Package\cr
Version: \tab 2023.1.0\cr
Date: \tab 2023-08-21\cr
License: \tab GPL-2 or greater\cr
LazyLoad: \tab yes\cr
}
Derivative-Free optimization algorithms. These algorithms do not require gradient information.  
More importantly, they can be used to solve non-smooth optimization problems. 
These algorithms were translated from the Matlab code of Prof. C.T. Kelley, given in his book "Iterative methods for optimization".
However, there are some non-trivial modifications of the algorithm. \cr

Currently, the Nelder-Mead and Hooke-Jeeves algorithms is implemented.  In future, more derivative-free algorithms may be added.
}
\author{
Ravi Varadhan, Johns Hopkins University   \cr
URL:  http://www.jhsph.edu/agingandhealth/People/Faculty_personal_pages/Varadhan.html  \cr
Hans W. Borchers, ABB Corporate Research \cr
Maintainer:  Ravi Varadhan <ravi.varadhan@jhu.edu>
}
\references{
C.T. Kelley (1999), Iterative Methods for Optimization, SIAM.
}
\keyword{optimize}
