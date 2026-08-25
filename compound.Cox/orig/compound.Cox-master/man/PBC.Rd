\name{PBC}
\alias{PBC}
\docType{data}
\title{
Primary biliary cirrhosis (PBC) of the liver data
}
\description{
A subset of primary biliary cirrhosis (PBC) of the liver data in the book "Counting Process & Survival Analysis" by Fleming & Harrington (1991).
This subset is used in Tibshirani (1997).
}
\usage{data(PBC)}
\format{
  A data frame with 276 observations on the following 19 variables.
  \describe{
    \item{\code{T}}{Survival times (either time to death or censoring) in days}
    \item{\code{d}}{Censoring indicator, 1=death, 0=censoring}
    \item{\code{trt}}{Treatment indicator, 1=treatment by D-penicillamine, 0=placebo}
    \item{\code{age}}{Age in years (days divided by 365.25)}
    \item{\code{sex}}{Sex, 0=male, 1=female}
    \item{\code{asc}}{Presence of ascites, 0=no, 1=yes}
    \item{\code{hep}}{Presence of hepatomegaly, 0=no, 1=yes}
    \item{\code{spi}}{Presence of spiders, 0=no, 1=yes}
    \item{\code{ede}}{Presence of edema, 0=no edema, 0.5=edema resolved by therapy, 1=edema not resolved by therapy}
    \item{\code{bil}}{log(bililubin, mg/dl)}
    \item{\code{cho}}{log(cholesterol, mg/dl)}
    \item{\code{alb}}{log(albumin, gm/dl)}
    \item{\code{cop}}{log(urine copper, mg/day)}
    \item{\code{alk}}{log(alkarine, U/liter)}
    \item{\code{SGO}}{log(SGOT, in U/ml)}
    \item{\code{tri}}{log(triglycerides, in mg/dl)}
    \item{\code{pla}}{log(platelet count, [the number of platelets per-cubic-milliliter of blood]/1000)}
    \item{\code{pro}}{log(prothrombin time, in seconds)}
    \item{\code{gra}}{Histologic stage of disease, graded 1, 2, 3, or 4}
  }
}
\details{
Survival data consisting of 276 patients with 17 covariates.
Among them, 111 patients died (d=1) while others were censored (d=0).
The covariates consist of a treatment indicator (trt), age, sex, 5 categorical variables
 (ascites, hepatomegaly, spider, edema, and stage of disease) and
 9 log-transformed continuous variables (bilirubin, cholesterol, albumin, urine copper,
 alkarine, SGOT, triglycerides, platelet count, and prothrombine). 
}
\source{
Fleming & Harrignton (1991); Tibshirani (1997)
}
\references{
Tibshirani R (1997), The Lasso method for variable selection in the Cox model, Statistics in Medicine, 385-395.
}
\examples{
data(PBC)
PBC[1:5,] ### profiles for the first 5 patients ###
# See also Appendix D.1 of Fleming & Harrington, Counting Process & Survival Analysis (1991) # 
}
\keyword{datasets}
\keyword{PBC}
