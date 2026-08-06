\name{bushmiss}
\alias{bushmiss}
\docType{data}
\title{ Campbell Bushfire Data with added missing data items}
\description{
  This data set is based on the bushfire data set which was used by 
  Campbell (1984) to locate bushfire scars - see  \code{\link[robustbase]{bushfire}}
  in package \code{robustbase}. The original dataset contains satelite 
  measurements on five frequency bands, corresponding to each of 38 pixels. 
}
\usage{data(bushmiss)}
\format{
  A data frame with 190 observations on 6 variables.
  
  The original data set consists of 38 observations in 5 variables. 
  Based on it four new data sets are created in which some of the data 
  items are replaced by missing values with a simple "missing completely 
  at random " mechanism. For this purpose independent Bernoulli trials 
  are realized for each data item with a probability of success 0.1, 0.2, 0.3, 0.4,
  where success means that the corresponding item is set to missing.  The obtained five 
  data sets, including the original one (each with probability of a data 
  item to be missing equal to 0, 0.1, 0.2, 0.3 and 0.4 which is reflected 
  in the new variable \code{MPROB}) are merged. (See also Beguin and Hulliger (2004).)
  %\describe{}
}
%\Note{}
\source{
  Maronna, R.A. and Yohai, V.J. (1995)
  The Behavoiur of the Stahel-Donoho Robust Multivariate Estimator.
  \emph{Journal of the American Statistical Association} \bold{90}, 330--341.
  
  Beguin, C. and Hulliger, B. (2004) Multivariate outlier detection in incomplete
    survey data: the epidemic algorithm and transformed rank correlations.
    \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)} \bold{127}, 2, 275--294.
}
%\seealso{}
\examples{
## The following code will result in exactly the same output
##  as the one obtained from the original data set
data(bushmiss)
bf <- bushmiss[bushmiss$MPROB==0,1:5]
plot(bf)
covMcd(bf)


\dontrun{
##  This is the code with which the missing data were created:
##
##  Creates a data set with missing values (for testing purposes)
##  from a complete data set 'x'. The probability of
##  each item being missing is 'pr' (Bernoulli trials).
##
getmiss <- function(x, pr=0.1)
{
    n <- nrow(x)
    p <- ncol(x)
    done <- FALSE
    iter <- 0
    while(iter <= 50){
        bt <- rbinom(n*p, 1, pr)
        btmat <- matrix(bt, nrow=n)
        btmiss <- ifelse(btmat==1, NA, 0)
        y <- x+btmiss
        if(length(which(rowSums(nanmap(y)) == p)) == 0)
            return (y)
        iter <- iter + 1
    }
    y
}
}

}
\keyword{datasets}
