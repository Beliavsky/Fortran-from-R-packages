\name{Wilks.test}
\alias{Wilks.test}
\alias{Wilks.test.default}
\alias{Wilks.test.data.frame}
\alias{Wilks.test.formula}
\alias{Wilks.test.matrix}
\alias{model.frame.Wilks.test}
\title{
Classical and Robust One-way MANOVA: Wilks Lambda
}
\description{
Classical and Robust One-way MANOVA: Wilks Lambda
}
\usage{

\method{Wilks.test}{formula}(formula, data, \dots, subset, na.action)

\method{Wilks.test}{default}(x, grouping, method=c("c", "mcd", "rank"), 
    approximation=c("Bartlett", "Rao", "empirical"),
    xd=NULL, xq=NULL, xfn = NULL, xwl=NULL, nrep=3000, trace=FALSE, \dots)

\method{Wilks.test}{data.frame}(x, \dots)

\method{Wilks.test}{matrix}(x, grouping, \dots, subset, na.action)
}
\arguments{
  \item{formula}{
    A formula of the form \code{groups ~ x1 + x2 + \dots}  That is, the
    response is the grouping factor and the right hand side specifies
    the (non-factor) variables.
  }
  \item{data}{
    Data frame from which variables specified in \code{formula} are
    to be taken.
  }
  \item{x}{
    (required if no formula is given as the principal argument.)
    a matrix or data frame or Matrix containing the explanatory variables.
  }
  \item{grouping}{
    grouping variable - a factor specifying the class for each 
    observation (required if no formula argument is given.)
  }
  \item{subset}{
    An index vector specifying the cases to be used. 
  }
  \item{na.action}{
    A function to specify the action to be taken if \code{NA}s are found.
    The default action is for the procedure to fail.  An alternative is
    \code{na.omit}, which leads to rejection of cases with missing values on
    any required variable. 
  }
  \item{method}{
    \code{"c"} for standard estimators of the mean and variance,
    \code{"mcd"} for MCD estimators of mean and variances and 
    \code{"rank"} for rank based wilks' lambda as proposed by Nath and Pavur (1985).
  }
  \item{approximation}{
    \code{"Bartlett"} for Bartlett approximation (default),
    \code{"Rao"} for rao approximation (only for method="c") and
    \code{"empirical"} for simulated empirical distribution.
  }
  \item{xd}{multiplication factor for the approximate distribution of 
    the robust Lambda statistic. If \code{xd=NULL} the factor will
    computed by simulation and will be returned in the value (see Details)
  }
  \item{xq}{the degrees of freedom for the approximate \eqn{\chi^2} distribution of 
    the robust Lambda statistic. If \code{xq=NULL} the degrees of freedom will
    computed by simulation and will be returned in the value (see Details)
  }
  \item{xfn}{the empirical distribution function. If \code{xfn=NULL} the empirical function  will
    be estimated by simulation and will be returned in the value (see Details)
  }
  \item{xwl}{the simulated values of the robust statistic. If \code{xwl=NULL} the simulation 
    will be performed and the calculated result will be returned in the value (see Details)
  }
  \item{nrep}{number of trials for the simulations for computing the 
    multiplication factor \code{xd} and the degrees of freedom \code{xq}. Default is \code{nrep=3000}.
  }
 \item{trace}{whether to print intermediate results. Default is \code{trace = FALSE}}  
 \item{\dots}{arguments passed to or from other methods.}
}
\value{
  A list with class \code{"htest"} containing the following components:
  \item{statistic}{the value of the Wilks' Lambda statistic.}
  \item{parameter}{The corresponding approximation of the Wilks' lambda statistic and the degrees of freedom.}
  \item{p.value}{the p-value for the test.}
  \item{estimate}{the estimated mean vectors.}
  \item{method}{a character string indicating what type of test was
    performed.}
  \item{data.name}{a character string giving the name of the data.}
  \item{xd}{multiplication factor for the approximate distribution of 
    the robust Lambda statistic. 
  }
  \item{xq}{the degrees of freedom for the approximate \eqn{\chi^2} distribution of 
    the robust Lambda statistic. 
  }
}
\details{
    The classical Wilks' Lambda statistic for testing the equality of 
    the group means of two or more groups is modified into a robust 
    one through substituting the classical estimates by the highly robust 
    and efficient reweighted MCD estimates, which can be computed efficiently 
    by the FAST-MCD algorithm - see \code{\link{CovMcd}}. 
    An approximation for the finite sample distribution of the
    Lambda statistic is obtained, based on matching the mean and 
    variance of a multiple of an \eqn{\chi^2} distribution which 
    are computed by simultaion.
}
\note{
This function may be called giving either a formula and
optional data frame, or a matrix and grouping factor as the first
two arguments.  All other arguments are optional.
}

\references{ 
 Todorov, V. and Filzmoser, P. (2007) Robust statistic for the one-way MANOVA,
 \emph{submetted to the Journal of Environmetrics}.

 Todorov, V. (2007) Robust selection of variables in linear discriminant 
 analysis, \emph{Statistical Methods and Applications}, \bold{15}, 
 395.407, doi:10.1007/s10260-006-0032-6.
 
 Nath, R. and Pavur, R. (1985) A new statistic in the one way 
 multivariate analysis of variance, \emph{Computatational 
 Statistics and Data Analysis}, \bold{2}, 297--315
}

\author{ Valentin Todorov \email{valentin.todorov@chello.at}
}

\seealso{ \code{\link{CovMcd}}, \code{\link{T2.test}}}

\examples{
library(MASS)
data(anorexia)
grp <- as.factor(anorexia[,1])
x <- as.matrix(anorexia[,2:3])
##  Using the default interface, classical test
Wilks.test(x, grouping=grp, method="c")

##  Using the default interface, rank based test
Wilks.test(x, grouping=grp, method="rank")

## For this data set: p=2, n=n1+n2+n3=29+26+17
## were computed the following multiplication factor xd and degrees of freedom xq
##  for the MCD estimates with alpha=0.5
xd <-  -0.02162666
xq <- 3.63971
Wilks.test(x, grouping=grp, method="mcd", xd=xd, xq=xq)

## Now the same with the formula interface
Wilks.test(Treat~Prewt+Postwt, data=anorexia, method="mcd", xd=xd, xq=xq)

##Iris data with formula interface
data(iris)
Wilks.test(Species~., data=iris, method="c")

## and with default interface
Wilks.test(iris[,1:4],grouping=iris[,5], method="c")

# hemophilia data - classical, rank and MCD test
data(hemophilia)
hemophilia$gr <- as.factor(hemophilia$gr)

Wilks.test(gr~., data=hemophilia, method="c")
Wilks.test(gr~., data=hemophilia, method="rank")
## already simulated parameters for MCD with alpha=0.5
xd <- -0.01805436
xq <- 1.950301
Wilks.test(gr~., data=hemophilia, xd=xd, xq=xq, method="mcd")

}
\keyword{multivariate}
\keyword{robust}
