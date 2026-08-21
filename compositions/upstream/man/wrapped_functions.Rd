\name{wrapped_functions}
\alias{anova}
%\alias{summary.manova}
\title{Standard R functions wrapped for compatibility}
\description{
  These functions of a standard R distribution (from package "base" or "stats") are wrapped by naive functions in "compositions" with the goal to ensure their normal behavior with compositional data objects. 
}
\usage{
    anova(\dots)
   % \method{summary}{manova}( object, \dots )
}

\arguments{
  %\item{object}{argument passed to the original `summary` method}
  \item{\dots}{arguments passed to the original function (all!)}
}
\value{
  The same as the original function from package base (i.e. search for it with `?base::anova`). 
  %or `?stats::summary.manova`).
}
\details{
  The functions documented in this page are just wrapers around base functions 
  from R that, due to a variety of reasons, need pre- or post-processing when 
  "compositions" is loaded. Pre-processing are e.g., converting "rmult" class 
  objects to plain "matrix" objects, or removing sticky class behaviour (see 
  \code{\link{getStickyClassOption}})
}

\author{Raimon Tolosana-Delgado \url{http://www.stat.boogaart.de}}
\seealso{
  \code{anova} in package "base" . %,   \code{summary.manova} in package "stats"

}


\examples{
# anova:
data("Hydrochem") # load data
Z = acomp(Hydrochem[,7:19]) # select composition
Hydrochem$compo = Z # attach to dataset
md = lm(alr(compo)~log(H), data=Hydrochem) # fit model
anova(md)   # anova test

}
\keyword{multivariate}
