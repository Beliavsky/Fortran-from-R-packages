\name{trees513}
\alias{trees513}
\encoding{utf-8}
\docType{data}
\title{ Frankonian Tree Damage Data }
\description{
  Damages on young trees caused by deer browsing.
}
\usage{data("trees513")}
\format{
  A data frame with 2700 observations on the following 4 variables.
  \describe{
    \item{\code{damage}}{a factor with levels \code{yes} and \code{no}
        indicating whether or not the trees has been damaged by game animals,
        mostly roe deer.}
    \item{\code{species}}{a factor with levels \code{spruce}, \code{fir}, 
         \code{pine}, \code{softwood (other)}, \code{beech}, \code{oak}, \code{ash/maple/elm/lime}, 
         and \code{hardwood (other)}.}
    \item{\code{lattice}}{a factor with levels \code{1}, ..., \code{53}, 
          essentially a number indicating the position of the sampled area.}
    \item{\code{plot}}{a factor with levels \code{x_1}, ..., \code{x_5} where
          \code{x} is the lattice. \code{plot}
          is nested within \code{lattice} and is a replication for each lattice point.}
  }
}
\details{

In most parts of Germany, the natural or artificial
regeneration of forests is difficult due to a high browsing
intensity. Young trees suffer from browsing damage, mostly by roe
and red deer. In order to estimate the browsing intensity for
several tree species, the Bavarian State Ministry of Agriculture
and Foresty conducts a survey every three years. Based on
the estimated percentage of damaged trees, suggestions for
the implementation or modification of deer management plans are made.
The survey takes place in all 756 game management districts
(`Hegegemeinschaften') in Bavaria. The data given here are from 
the game management district number 513 `Unterer Aischgrund'
(located in Frankonia between Erlangen and Höchstadt) in 2006. The data
of 2700 trees include the species and a binary
variable indicating whether or not the tree suffers from damage caused
by deer browsing.

Data and analysis are presented by \bibcitet{hothorn2008b}.

}
\references{\bibshow{*}}
\examples{

  summary(trees513)

}
\keyword{datasets}
