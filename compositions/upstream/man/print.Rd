\name{print.acomp}
\alias{print.acomp}
\alias{print.rcomp}
\alias{print.rplus}
\alias{print.aplus}
\title{Printing compositional data.}
\description{
  Prints compositional objects with appropriate missing encodings.
}
\usage{
\method{print}{acomp}(x,...,replace0=TRUE)
\method{print}{aplus}(x,...,replace0=TRUE)
\method{print}{rcomp}(x,...,replace0=FALSE)
\method{print}{rplus}(x,...,replace0=FALSE)
}
\arguments{
  \item{x}{a compositional object}
  \item{\ldots}{further arguments to \code{print.default}}
  \item{replace0}{logical: Shall 0 be treated as "Below detection Limit"
    with unkown limit.}
}
\value{
An invisible version of x.
}

\details{
  Missings are displayed with an appropriate encoding:
  \describe{
    \item{MAR}{Missing at random: The value is missing independently of
      its true value.}
    \item{MNAR}{Missing NOT at random: The value is missing dependently of
      its true value, but without a known systematic. Maybe a better name
      would be: Value dependen missingness.}
    \item{BDL}{below detection limit (with unspecified detection limit):
      The value is missing because it was below an unkown detection limit.}
    \item{<Detectionlimit}{below detection limit (with specified detection limit):
      The value is below the displayed detection limit. }
    \item{SZ}{Structural Zero: A true value is either bound to be zero or
      does not exist for structural
      nonrandom reasons. E.g. the
      portion of pregnant girls at a boys school.}
    \item{ERR}{Error: An illegal encoding value was found in the
      object. }
  }
}
\section{Missing Policy}{
  The policy of treatment of zeroes, missing values and values 
  below detecion limit is explained in depth in \link{compositions.missings}. 
  }

\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}, Raimon
Tolosana-Delgado}

\seealso{
  \code{\link{clr}},\code{\link{acomp}}, 
  \code{\link{plot.acomp}}, \code{\link{boxplot.acomp}},
  \code{\link{barplot.acomp}}, \code{\link{mean.acomp}},
  \code{\link{var.acomp}}, \code{\link{variation.acomp}},
  \code{\link{zeroreplace}}
}

\references{
  Boogaart, K.G. v.d., R. Tolosana-Delgado, M. Bren (2006) Concepts for 
  handling of zeros and missing values in compositional data, in: 
  E. Pirard (ed.) (2006)Proceedings of the IAMG'2006 Annual Conference on 
  "Quantitative Geology from multiple sources", September 2006, Liege, 
  Belgium,, S07-01, 4pages, ISBN 978-2-9600644-0-7

}
\examples{
data(SimulatedAmounts)
mydata <- simulateMissings(sa.groups5,dl=0.01,knownlimit=TRUE,
                     MAR=0.05,MNARprob=0.05,SZprob=0.05)
mydata[1,1]<-BDLvalue
print(aplus(mydata))
print(aplus(mydata),digits=3)
print(acomp(mydata))
print(rplus(mydata))
print(rcomp(mydata))

}
\keyword{classes}
