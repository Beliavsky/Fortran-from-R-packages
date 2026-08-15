\name{Stillbirths}
\alias{data_sbirth}
\title{Still births of New Zealand white rabbits}
\description{The number of
stillbirths of New Zealand white rabbits.
}
\usage{
data_sbirth
}
\arguments{
\item{data_sbirth}{
A vector of (non-negative integer) count values.
}
}

\details{
The data set consists of frequencies of still births in 402
litters of New Zealand white rabbits, originally used by
Morgan et al. (2007).
}

\value{
data_sbirth gives the frequencies of still births in 402
litters of New Zealand white rabbits.
}

\references{
Morgan, B. T., Palmer, K. J., & Ridout, M. S. (2007). Negative score test statistic. The American Statistician, 61(4), 285-288.

Alshkaki, R. S. A. (2016). On the zero-one inflated Poisson distribution. Int J Stat Distrib Appl, 2(4), 42-8.

Pudprommarat, C. (2020, March). Zero-one inflated negative binomial-Sushila distribution and its application. In International Academic Multidisciplinary Research Conference in Rome 2020 (pp. 20-28).
}

\author{
Muhammad Imran and M.H. Tahir.

R implementation and documentation: Muhammad Imran \email{imranshakoor84@yahoo.com} and M.H. Tahir <mht@iub.edu.pk>.
}



%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{data_criminal}
}
}


\examples{
x <- data_sbirth
table (x)
}



