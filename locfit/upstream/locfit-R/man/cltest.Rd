\name{cltest}
\alias{cltest}
\title{Test dataset for classification}
\usage{data(cltest)}
\format{
Data Frame. Three variables x1, x2 and y. The latter indicates
class membership.
}
\description{
200 observations from a 2 population model. Under population 0,
\eqn{x_{1,i}} has a standard normal distribution, and
\eqn{x_{2,i} = (2-x_{1,i}^2+z_i)/3}, where \eqn{z_i} is also standard normal.
Under population 1, \eqn{x_{2,i} = -(2-x_{1,i}^2+z_i)/3}.
The optimal classification regions form a checkerboard pattern,
with horizontal boundary at \eqn{x_2=0}, vertical boundaries at
\eqn{x_1 = \pm \sqrt{2}}.

This is the same model as the cltrain dataset.
}
\keyword{datasets}
