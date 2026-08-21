\name{binary}
\alias{binary}
\alias{bit}
\alias{bit.numeric}
\alias{bit.character}
\alias{bit<-}
\alias{bit<-.numeric}
\alias{bit<-.character}
\alias{maxBit}
\alias{maxBit.numeric}
\alias{maxBit.character}
\alias{bitCount}
\alias{gsi.orSum}
\alias{whichBits}
\alias{binary2logical}
\alias{unbinary}
\alias{gadic}
\title{Treating binary and g-adic numbers}
\description{
  Allows the access to individual digits in binary (and general
  g-adic) numbers. 
}
\usage{
binary(x,mb=max(maxBit(x,g)),g=2)
unbinary(x,g=2)
bit(x,b,g=2)                       
\method{bit}{numeric}(x,b=0:maxBit(x,g),g=2)
\method{bit}{character}(x,b=0:maxBit(x,g),g=2)
bit(x,b,g=2) <- value                      
\method{bit}{numeric}(x,b=0:maxBit(x,g),g=2) <- value
\method{bit}{character}(x,b=0:maxBit(x,g),g=2) <- value
maxBit(x,g=2)
\method{maxBit}{numeric}(x,g=2)
\method{maxBit}{character}(x,g=2)
bitCount(x,mb=max(maxBit(x,g)),g=2)
gsi.orSum(...,g=2)
whichBits(x,mb=max(maxBit(x,g)),g=2,values=c(TRUE))
binary2logical(x,mb=max(maxBit(x,g)),g=2,values=c(TRUE))
}
\arguments{
  \item{x}{a number either represented a g-adic character string or as a
    integeral numeric value}
  \item{b}{the indicies of the bits to be processes. The least
    significant bit has index 0.}
  \item{mb}{maximal bit. The index of the most significant bit to be treated}
  \item{g}{the base of the g-adic representation. 2 corresponds to
    binary numbers, 8 to octal numbers, 16 to hexadecimal numbers. g is
    limited by 36. }
  \item{value}{ a vector of bit values to be selected or setted. }
  \item{values}{ a vector of bit values that should be considered as TRUE. }
  \item{\dots}{some binary numbers}
}
\value{
\item{binary}{returns a standard binary (or g-adic) character representation of
  the number}
\item{unbinary}{returns a binary (or g-adic) representation of
  the number}
\item{bit}{returns the values of the requested bits. The values are
  returned as a logical vector for binary numbers an as numeric digit
  values for other g-adic numbers.}
\item{maxBit}{returns the most significant bit represented in the number. This is
the highest bit set in numeric numbers and the highest actually given
character in a character representation. }
\item{bitCount}{returns the g-adic crossfoot of the number. For a binary
  number this is the number of bits set}
\item{gsi.orSum}{ Only works for binary numbers and does a parallel or
  on each of the bits for a list of binary numbers.}
\item{whichBits}{returns the indices of the bits acutally set (or more
  precisely of the bits with value in \code{values})}
\item{binary2logical}{returns the a true false vector of the bits
acutally set (or more
  precisely of the bits with value in \code{values})}
}
\details{
  These routines are primerily intended to manipulate g-adic numbers for
  user interface purposes and condensed representation of
  information. They are not intended for a long number arithmetic. 
  
}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{outlierplot}}
}

\examples{
(x<-unbinary("10101010"))
(y<-binary(x))
bit(x,1:3)
bit(y,0:3)
maxBit(x)
maxBit(y)
whichBits(x)
whichBits(y)
binary2logical(y)
bit(x) %*% 2^(0:maxBits(x))
bit(y) %*% 2^(0:maxBits(y))
gsi.orSum(y,1)
bitCount(x)
bitCount(y)
bit(x,2)<-1
x
bit(y,2)<-1
y

}
\keyword{logic}
