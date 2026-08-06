%  file sn/man/barolo.Rd  
%  This file is a component of the package 'sn' for R
%  copyright (C) 2013 Adelchi Azzalini
%---------------------
\name{barolo}
\alias{barolo}
\docType{data}
\title{Price of Barolo wine}

\description{A data frame with prices of bottles of Barolo wine 
      and some auxiliary variables}

\usage{data(barolo)}
\format{A data frame with 307 observations on five variables, as follows:
\tabular{ll}{%
 \code{reseller}\tab reseller code (factor with levels \code{A, B, C, D}) \cr
 \code{vintage} \tab vintage year (numeric) \cr
 \code{volume}  \tab content volume in centilitres (numeric) \cr
 \code{price}   \tab price in Euro (numeric) \cr
 \code{age}     \tab age in 2010 (numeric)
}
For six items, \code{vintage} is \code{NA}'s and so also \code{age}.
Three items have a non-standard volume of 50 cl.
}
\details{The data have been obtained in July 2010 from the websites  
of four Italian wine resellers, selecting only quotations of Barolo 
wine, which is produced in the Piedmont region of Italy. 
The price does not include the delivery charge.

The data have been presented in Section 4.3.2 of the reference below,
where a subset of them has been used for illustrative purposes.
This subset refers to reseller \code{"A"} and bottles of 75cl.
}

\source{
 Azzalini, A. with the collaboration of Capitanio, A. (2014). 
 \emph{The Skew-Normal and Related Families}. 
 Cambridge University Press, IMS Monographs series.
}

\examples{
data(barolo)
attach(barolo)  
f <- cut(age, c(0, 5,  6, 8, 11,  30))
table(volume, f)
plot(volume, price, col=as.numeric(f), pch=as.character(reseller))
legend(400, 990, col=1:5, lty=1, title="age class",
       legend=c("4-5", "6", "7-8", "9-11", "12-30"))
#
A75 <- (reseller=="A" & volume==75)
hist(log(price[A75],10), col="gray85")
# see Figure 4.7 of the source
}
\keyword{datasets}
