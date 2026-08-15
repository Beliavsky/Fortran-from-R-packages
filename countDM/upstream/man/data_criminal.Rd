\name{Criminal acts}
\alias{data_criminal}
\title{The crime sociology consisting a sample of 4301 people with deviating behavior}
\description{The data set from crime sociology consisting a sample of 4301 people with deviating behavior.
}

\usage{
data_criminal
}
\arguments{
\item{data_criminal}{
A vector of (non-negative integer) count values.
}
}

\details{
The data set from crime sociology consisting a sample of 4301 people with deviating behavior. Recently, it is used by Zhang et al. (2016), fitted the zero one inflated Poisson distribution.
}

\value{
data_criminal gives the crime sociology consisting a sample of 4301 people with deviating behavior.
}

\references{
Zhang, C., Tian, G. L., & Ng, K. W. (2016). Properties of the zero-and-one inflated Poisson distribution and likelihood-based inference methods. Statistics and its interface, 9(1), 11-32.
}

\author{
Muhammad Imran and M.H. Tahir.

R implementation and documentation: Muhammad Imran \email{imranshakoor84@yahoo.com} and M.H. Tahir <mht@iub.edu.pk>.
}



%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{data_sbirth}
}
}


\examples{
x <- data_criminal
table (x)
}



