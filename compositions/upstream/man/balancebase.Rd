\name{balance}
\alias{balanceBase}
\alias{balance01}
\alias{balance01.acomp}
\alias{balance01.rcomp}
\alias{balance}
\alias{balance.acomp}
\alias{balance.rcomp}
\alias{balance.aplus}
\alias{balance.rplus}
\alias{balanceBase.acomp}
\alias{balanceBase.rcomp}
\alias{balanceBase.aplus}
\alias{balanceBase.rplus}
\title{Compute balances for a compositional dataset.}
\description{
  Compute balances in a compositional dataset.
}
\usage{
balance(X,...)
\method{balance}{acomp}(X,expr,...)
\method{balance}{rcomp}(X,expr,...)
\method{balance}{aplus}(X,expr,...)
\method{balance}{rplus}(X,expr,...)
balance01(X,...)
\method{balance01}{acomp}(X,expr,...)
\method{balance01}{rcomp}(X,expr,...)
balanceBase(X,...)
\method{balanceBase}{acomp}(X,expr,...)
\method{balanceBase}{rcomp}(X,expr,...)
\method{balanceBase}{acomp}(X,expr,...)
\method{balanceBase}{rcomp}(X,expr,...)

}
\arguments{
  \item{X}{ compositional dataset (or optionally just its column names
    for balanceBase)}
  \item{expr}{ a \code{~} formula using the column names of X as
variables and
    seperating them by \code{/} and organize by paranthesis
    \code{()}. \code{:} and \code{*} can be used instead of \code{/}
    when the corresponding balance should not be created. \code{-} can
    be used as an synonym to \code{/} in the real geometries. \code{1}
    can be used in the unclosed geometries to level against a constant. 
  }
  \item{\dots}{for future perposes}
}
\value{
  \item{balance}{a matrix (or vector) with the corresponding balances
    of the dataset.}
  \item{balance01}{a matrix (or vector) with the corresponding balances
    in the dataset transformed in the given geometry to a value between
    0 and 1.}
  \item{balanceBase}{a matrix (or vector) with column vectors giving the
    transform in the cdt-transform used to achiev the correponding
    balances. }
}
\details{
  For \code{acomp}-compositions balances are defined as orthogonal
  projections representing the log ratio of the geometric means of
  subsets of elements. Based on a recursive subdivision (provided by the
  \code{expr=}) this projections provide a (complete or incomplete) basis of
  the clr-plane. The basis is given by the \code{balanceBase}
  functions. The transform is given by the \code{balance} functions. The
  \code{balance01} functions are a backtransform of the balances to the
  amount of the first portion if this was the only balance in a 2
  element composition, providing an "interpretation" for the values of
  the balances.

  The package tries to give similar concepts for the other scales.
  For \code{rcomp} objects the concept is mainly unchanges but augmented
  by a virtual component 1, which always has portion 1. 


  For
  \code{rcomp} objects, we choose not a "orthogonal" transformation
  since such a concept anyway does not really exist in the given space,
  but merily use the difference of one subset to the other. The
  balance01 is than not really a transform of the balance but simply the
  portion of the first group of parts in all contrasted parts. 

  For \code{rplus} objects we just used an analog to generalisation from
  the \code{rcomp} defintion as \code{aplus} is generalized from
  \code{acomp}. However at this time we have no idea wether this has any
  usefull interpretation.
  
  
} \seealso{
  \code{\link{clr}},\code{\link{ilr}},\code{\link{ipt}}, \code{\link{ilrBase}}
}
\references{
  \url{https://ima.udg.edu/Activitats/CoDaWork08/} Papers of Boogaart and Tolosana
  \url{https://ima.udg.edu/Activitats/CoDaWork05/} Paper of Egozcue
  }
\examples{
X <- rnorm(100)
Y <- rnorm.acomp(100,acomp(c(A=1,B=1,C=1)),0.1*diag(3))+acomp(t(outer(c(0.2,0.3,0.4),X,"^")))
colnames(Y) <- c("A","B","C")

subComps <- function(X,...,all=list(...)) {
  X <- oneOrDataset(X)
  nams <- sapply(all,function(x) paste(x[[2]],x[[3]],sep=","))
  val  <- sapply(all,function(x){
              a = X[,match(as.character(x[[2]]),colnames(X)) ]
              b = X[,match(as.character(x[[2]]),colnames(X)) ]
              c = X[,match(as.character(x[[3]]),colnames(X)) ] 
              return(a/(b+c))
             }) 
  colnames(val)<-nams
  val
}
pairs(cbind(ilr(Y),X),panel=function(x,y,...) {points(x,y,...);abline(lm(y~x))})
pairs(cbind(balance(Y,~A/B/C),X),panel=function(x,y,...) {points(x,y,...);abline(lm(y~x))})

pairwisePlot(balance(Y,~A/B/C),X)
pairwisePlot(X,balance(Y,~A/B/C),panel=function(x,y,...) {plot(x,y,...);abline(lm(y~x))})
pairwisePlot(X,balance01(Y,~A/B/C))
pairwisePlot(X,subComps(Y,A~B,A~C,B~C))



balance(rcomp(Y),~A/B/C)
balance(aplus(Y),~A/B/C)
balance(rplus(Y),~A/B/C)



}
\keyword{multivariate}
