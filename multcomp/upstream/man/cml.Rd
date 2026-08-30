\name{cml}
\alias{cml}
\docType{data}
\title{Chronic Myelogenous Leukemia survival data.}
\description{
  Survival in a randomised trial comparing three treatments for 
  Chronic Myelogeneous Leukemia (simulated data).
}
\usage{data("cml")}
\format{
  A data frame with 507 observations on the following 7 variables.
  \describe{
    \item{\code{center}}{a factor with 54 levels indicating the study center.}
    \item{\code{treatment}}{a factor with levels \code{trt1}, \code{trt2}, \code{trt3} indicating the treatment group.}
    \item{\code{sex}}{sex (0 = female, 1 = male)}
    \item{\code{age}}{age in years}
    \item{\code{riskgroup}}{risk group (0 = low, 1 = medium, 2 = high)}
    \item{\code{status}}{censoring status (FALSE = censored, TRUE = dead)}
    \item{\code{time}}{survival or censoring time in days.}
    }
}
\details{
The data are simulated according to structure of the data by the German CML 
Study Group used in \bibcitet{hehlmann1994}.
}
\references{\bibshow{*}}
\examples{

if (require("coxme")) {
    data("cml")
    ### one-sided simultaneous confidence intervals for many-to-one 
    ### comparisons of treatment effects concerning time of survival 
    ### modeled by a frailty Cox model with adjustment for further 
    ### covariates and center-specific random effect.
    cml_coxme <- coxme(Surv(time, status) ~ treatment + sex + age + riskgroup + (1|center), 
                       data = cml)
    glht_coxme <- glht(model = cml_coxme, linfct = mcp(treatment = "Dunnett"), 
                       alternative = "greater")
    ci_coxme <- confint(glht_coxme)
    exp(ci_coxme$confint)[1:2,]
}
}
\keyword{datasets}
