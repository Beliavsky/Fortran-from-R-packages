\name{subsetting}
\alias{[.acomp}
\alias{[.rcomp}
\alias{[.ccomp}
\alias{[.aplus}
\alias{[.rplus}
\alias{[.rmult}
\alias{$.acomp}
\alias{$.rcomp}
\alias{$.ccomp}
\alias{$.aplus}
\alias{$.rplus}
\alias{$.rmult}
\alias{getStickyClassOption}
\alias{setStickyClassOption}
\title{Subsetting of compositions}
\description{
  Extract subsets (rows) or subsompositions (columns) of a compositional data set
}
\usage{
  getStickyClassOption()
  setStickyClassOption(value)
\method{[}{acomp}(x, i, j, drop=gsi.LengthOne(j))
\method{[}{rcomp}(x, i, j, drop=gsi.LengthOne(j))
\method{[}{aplus}(x, i, j, drop=gsi.LengthOne(j))
\method{[}{rplus}(x, i, j, drop=gsi.LengthOne(j))
\method{[}{ccomp}(x, i, j, drop=gsi.LengthOne(j))
\method{[}{rmult}(x, i, j, drop=gsi.LengthOne(j))
\method{$}{acomp}(x, name) 
\method{$}{rcomp}(x, name)
\method{$}{aplus}(x, name) 
\method{$}{rplus}(x, name)
\method{$}{ccomp}(x, name) 
\method{$}{rmult}(x, name)
}
\arguments{
  \item{x}{vector or dataset of a compositions class}
  \item{i}{row indices/names to select/exclude, resp. boolean of fitting length (recyling appplied if length(i)<nrow(x)); if x is a compositional vector, this gives the elements (equivalent to variables) to be extracted/selected}
  \item{j}{column indices/names to select/exclude, resp. boolean of fitting length (recyling appplied if length(i)<ncol(x))}
  \item{drop}{boolean, should matrices be simplified to vectors? defaults to FALSE (a difference with standard R). If set to TRUE, it has the extra effect of removing the compositional class}
  \item{name}{column name of the variable to be extracted OR name of a scaling function to be applied. It accepts \code{unclass} (or \code{raw}), \code{clo}, \code{pwlr}, \code{alr}, \code{clr}, \code{ilr}, \code{apt}, \code{cpt},
  \code{ipt}, \code{iit}, \code{ilt}, \code{cdt}, \code{idt} and their inverses \code{*Inv}}
  \item{value}{logical, controlling the global options for sticky classes}
}
\value{
  For \code{[} a vector or matrix with the relevant elements selected. 
  When selecting rows, this object is of the same class than \code{x}, i.e. the class is sticky.
  When selecting columns, the class depends on the number of columns selected and the value of \code{drop}. With \code{drop=T}, output is always a matrix or a vector. The same happens if 
  \code{gsi.LengthOne(j)==TRUE}, which happens if and only if \code{j} is a non-null vector of
  length one (i.e. if you only want one single column).
  
  If you want to get rid of sticky classes and return to the behaviour of "compositions" v1.xx, call \code{setStickyClassOption(FALSE)}. This may be a good idea if you run old scripts written for that versions of "compositions". You can recover the default behaviour from "compositions" v2 with \code{setStickyClassOption(TRUE)}, and check which sticky class status is currently defined in the global options with \code{getStickyClassOption()}.
  
  For \code{$} the output is either a transformed data set of the appropriate class, or the selected column as a class-less vector. The transformation ability is particularly useful if you have put a whole compositional class into one column of a data set, in which case you can confortably use the transformations in formula interfaces (see example below). This is NEVER sticky.
}
\author{R. Tolosana-Delgado}
\seealso{
  \code{\link{rmult}}, \code{\link{acomp}}, \code{\link{rcomp}},
  \code{\link{aplus}},
  \code{\link{rplus}}, \code{\link{ccomp}},
}

\examples{
data(Hydrochem)
xc = acomp(Hydrochem[,6:10])
xc[1:3,]
xc[-(10:nrow(xc)),]
xc[1:3,1:3]
xc[1:3,1:3, drop=TRUE]

xc[1:3,1]
class(xc[1:4,1])
class(xc[1:4,1, drop=TRUE])

data("Hydrochem")
xc = acomp(Hydrochem[, 9:14])
Hydrochem$compo = xc
lm(compo$clr~River, data=Hydrochem)
}
\keyword{multivariate}
