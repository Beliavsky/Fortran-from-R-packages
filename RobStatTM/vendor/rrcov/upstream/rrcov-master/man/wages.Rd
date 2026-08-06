\name{wages}
\alias{wages}
\docType{data}
\title{Wages and Hours}
\description{
The data are from a national sample of 6000 households with a male head earning less than 
USD 15,000 annually in 1966. The data were clasified into 39 demographic groups for 
analysis. The study was undertaken in the context of proposals for a guaranteed 
annual wage (negative income tax). At issue was the response of labor supply 
(average hours) to increasing hourly wages. 
The study was undertaken to estimate this response from available data.
}
\usage{data(wages)}
\format{
    A data frame with 39 observations on the following 10 variables: 
 
    \describe{
        \item{\code{HRS}}{Average hours worked during the year}
        \item{\code{RATE}}{Average hourly wage (USD)}
        \item{\code{ERSP}}{Average yearly earnings of spouse (USD)}
        \item{\code{ERNO}}{Average yearly earnings of other family members (USD)}
        \item{\code{NEIN}}{Average yearly non-earned income}
        \item{\code{ASSET}}{Average family asset holdings (Bank account, etc.) (USD)}
        \item{\code{AGE}}{Average age of respondent}
        \item{\code{DEP}}{Average number of dependents}
        \item{\code{RACE}}{Percent of white respondents}
        \item{\code{SCHOOL}}{Average highest grade of school completed}
    }
}
%%\details{}
\source{
    DASL library
%%    \url{http://lib.stat.cmu.edu/DASL/Datafiles/wagesdat.html}
    'http://lib.stat.cmu.edu/DASL/Datafiles/wagesdat.html'
}
\references{
    D.H. Greenberg and M. Kosters, (1970). Income Guarantees and the Working Poor, The Rand Corporation. 
}
\examples{
    data(wages)
    names(wages)
    x <- as.matrix(wages)
    ok <- is.finite(x \%*\% rep(1, ncol(x)))
    wages <- wages[ok, , drop = FALSE]
    wages.lm <- lm(HRS~AGE, data=wages)
    plot(HRS ~ AGE, data = wages)
    abline(wages.lm)
    class(wages.lm)
    names(wages.lm)
    summary(wages.lm)
    
    wages.mm <- lmrob(HRS~AGE, data=wages)
    plot(HRS ~ AGE, data = wages)
    abline(wages.mm)
    class(wages.mm)
    names(wages.mm)
    summary(wages.mm)   
}
\keyword{datasets}
