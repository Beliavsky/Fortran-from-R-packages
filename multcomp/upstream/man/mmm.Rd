\name{mmm}

\alias{mmm}
\alias{mlf}

\title{Simultaneous Inference for Multiple Marginal Models}

\description{Calculation of correlation between test statistics from multiple marginal models 
             using the score decomposition}

\usage{
mmm(...)
mlf(...)
}

\arguments{
  \item{\dots}{A names argument list containing fitted models (\code{mmm}) or
               definitions of linear functions (\code{mlf}). If
               only one linear function is defined for \code{mlf},
               it will be applied to all models in \code{mmm} by
               \code{\link{glht.mlf}}.}
}

\details{

  Estimated correlations of the estimated parameters of interest from the
  multiple marginal models are obtained using a stacked version of the
  i.i.d.  decomposition of parameter estimates by means of score components
  (first derivatives of the log likelihood).  The method is less
  conservative than the Bonferroni correction.  The details are provided by
  Pipper, Ritz and Bisgaard (2012).

  The implementation assumes that the model were fitted to the same data,
  i.e., the rows of the matrices returned by \code{estfun} belong to the
  same observations for each model.

  The reference distribution is always multivariate normal, if you want
  to use the multivariate t, please specify the corresponding degrees of
  freedom as an additional \code{df} argument to \code{\link{glht}}.

  Observations with missing values contribute zero to the score function.
  Models have to be fitted using \code{\link{na.exclude}} as \code{na.action}
  argument.

}

\value{

  The procedure implements multiple marginal models as introduced by
  \bibcitet{multcomp::pipper2011}.

  An object of class \code{mmm} or \code{mlf}, basically a named list of the
  arguments with a special method for \code{\link{glht}} being available for
  the latter.  \code{vcov}, \code{\link[sandwich]{estfun}}, and
  \code{\link[sandwich]{bread}} methods are available for objects of class
  \code{mmm}.
}
\references{\bibshow{*}}
\author{
  Code for the computation of the joint covariance and 
  sandwich matrices was contributed by Christian Ritz and 
  Christian B. Pipper.
}
\examples{

### replicate analysis of Hasler & Hothorn (2011), 
### A Dunnett-Type Procedure for Multiple Endpoints,
### The International Journal of Biostatistics: Vol. 7: Iss. 1, Article 3.
### DOI: 10.2202/1557-4679.1258

library("sandwich")

### see ?coagulation
if (require("SimComp")) {
    data("coagulation", package = "SimComp")

    ### level "S" is the standard, "H" and "B" are novel procedures
    coagulation$Group <- relevel(coagulation$Group, ref = "S")

    ### fit marginal models
    (m1 <- lm(Thromb.count ~ Group, data = coagulation))
    (m2 <- lm(ADP ~ Group, data = coagulation))
    (m3 <- lm(TRAP ~ Group, data = coagulation))

    ### set-up Dunnett comparisons for H - S and B - S 
    ### for all three models
    g <- glht(mmm(Thromb = m1, ADP = m2, TRAP = m3),
              mlf(mcp(Group = "Dunnett")), alternative = "greater")

    ### joint correlation
    cov2cor(vcov(g))

    ### simultaneous p-values adjusted by taking the correlation
    ### between the score contributions into account
    summary(g)
    ### simultaneous confidence intervals
    confint(g)

    ### compare with
    \dontrun{
        library("SimComp")
        SimCiDiff(data = coagulation, grp = "Group",
                  resp = c("Thromb.count","ADP","TRAP"), 
                  type = "Dunnett", alternative = "greater",
                  covar.equal = TRUE)
    }
 
    ### use sandwich variance matrix
    g <- glht(mmm(Thromb = m1, ADP = m2, TRAP = m3),
              mlf(mcp(Group = "Dunnett")), 
              alternative = "greater", vcov. = sandwich)
    summary(g)
    confint(g)
}

### attitude towards science data
data("mn6.9", package = "TH.data")

### one model for each item
mn6.9.y1 <- glm(y1 ~ group, family = binomial(), 
                na.action = na.omit, data = mn6.9)
mn6.9.y2 <- glm(y2 ~ group, family = binomial(), 
                na.action = na.omit, data = mn6.9)
mn6.9.y3 <- glm(y3 ~ group, family = binomial(), 
                na.action = na.omit, data = mn6.9)
mn6.9.y4 <- glm(y4 ~ group, family = binomial(), 
                na.action = na.omit, data = mn6.9)

### test all parameters simulaneously
summary(glht(mmm(mn6.9.y1, mn6.9.y2, mn6.9.y3, mn6.9.y4), 
             mlf(diag(2))))
### group differences
summary(glht(mmm(mn6.9.y1, mn6.9.y2, mn6.9.y3, mn6.9.y4), 
             mlf("group2 = 0")))

### alternative analysis of Klingenberg & Satopaa (2013),
### Simultaneous Confidence Intervals for Comparing Margins of
### Multivariate Binary Data, CSDA, 64, 87-98
### http://dx.doi.org/10.1016/j.csda.2013.02.016

### see supplementary material for data description
### NOTE: this is not the real data but only a subsample
influenza <- structure(list(
HEADACHE = c(1L, 0L, 0L, 1L, 0L, 0L, 1L, 1L, 1L,
0L, 0L, 1L, 0L, 1L, 0L, 1L, 1L, 1L, 1L, 1L, 1L, 0L, 0L, 0L, 0L,
1L, 1L), MALAISE = c(0L, 0L, 1L, 1L, 0L, 1L, 1L, 1L, 0L, 1L,
0L, 0L, 1L, 1L, 0L, 0L, 1L, 0L, 1L, 0L, 1L, 0L, 1L, 1L, 0L, 1L,
0L), PYREXIA = c(0L, 0L, 0L, 0L, 0L, 1L, 0L, 1L, 0L, 0L, 1L,
1L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 1L, 1L, 1L, 0L, 0L, 0L, 1L, 1L
), ARTHRALGIA = c(0L, 0L, 0L, 0L, 1L, 0L, 1L, 0L, 1L, 1L, 0L,
0L, 1L, 1L, 0L, 1L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L
), group = structure(c(2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L,
2L, 2L, 2L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
1L), levels = c("pla", "trt"), class = "factor"), Freq = c(32L,
165L, 10L, 23L, 3L, 1L, 4L, 2L, 4L, 2L, 1L, 1L, 1L, 1L, 167L,
1L, 11L, 37L, 7L, 7L, 5L, 3L, 3L, 1L, 2L, 4L, 2L)), names = c("HEADACHE",
"MALAISE", "PYREXIA", "ARTHRALGIA", "group", "Freq"), row.names = c(1L,
2L, 3L, 5L, 9L, 36L, 43L, 50L, 74L, 83L, 139L, 175L, 183L, 205L,
251L, 254L, 255L, 259L, 279L, 281L, 282L, 286L, 302L, 322L, 323L,
366L, 382L), class = "data.frame")
influenza <- influenza[rep(1:nrow(influenza), influenza$Freq), 1:5]

### Fitting marginal logistic regression models
(head_logreg <- glm(HEADACHE ~ group, data = influenza, 
                    family = binomial()))
(mala_logreg <- glm(MALAISE ~ group, data = influenza, 
                    family = binomial()))
(pyre_logreg <- glm(PYREXIA ~ group, data = influenza, 
                    family = binomial()))
(arth_logreg <- glm(ARTHRALGIA ~ group, data = influenza, 
                    family = binomial()))

### Simultaneous inference for log-odds
xy.sim <- glht(mmm(head = head_logreg,
                   mala = mala_logreg,
                   pyre = pyre_logreg,
                   arth = arth_logreg),
               mlf("grouptrt = 0"))
summary(xy.sim)
confint(xy.sim)

### Artificial examples
### Combining linear regression and logistic regression
set.seed(29)
y1 <- rnorm(100)
y2 <- factor(y1 + rnorm(100, sd = .1) > 0)
x1 <- gl(4, 25) 
x2 <- runif(100, 0, 10)

m1 <- lm(y1 ~ x1 + x2)
m2 <- glm(y2 ~ x1 + x2, family = binomial())
### Note that the same explanatory variables are considered in both models
### but the resulting parameter estimates are on 2 different scales 
### (original and log-odds scales)

### Simultaneous inference for the same parameter in the 2 model fits
summary(glht(mmm(m1 = m1, m2 = m2), mlf("x12 = 0")))

### Simultaneous inference for different parameters in the 2 model fits
summary(glht(mmm(m1 = m1, m2 = m2),
             mlf(m1 = "x12 = 0", m2 = "x13 = 0")))

### Simultaneous inference for different and identical parameters in the 2
### model fits
summary(glht(mmm(m1 = m1, m2 = m2),
             mlf(m1 = c("x12 = 0", "x13 = 0"), m2 = "x13 = 0")))

### Examples for binomial data
### Two independent outcomes
y1.1 <- rbinom(100, 1, 0.45)
y1.2 <- rbinom(100, 1, 0.55)
group <- factor(rep(c("A", "B"), 50))

m1 <- glm(y1.1 ~ group, family = binomial)
m2 <- glm(y1.2 ~ group, family = binomial)

summary(glht(mmm(m1 = m1, m2 = m2), 
             mlf("groupB = 0")))

### Two perfectly correlated outcomes
y2.1 <- rbinom(100, 1, 0.45)
y2.2 <- y2.1
group <- factor(rep(c("A", "B"), 50))

m1 <- glm(y2.1 ~ group, family = binomial)
m2 <- glm(y2.2 ~ group, family = binomial)

summary(glht(mmm(m1 = m1, m2 = m2), 
             mlf("groupB = 0")))

### use sandwich covariance matrix
summary(glht(mmm(m1 = m1, m2 = m2), 
             mlf("groupB = 0"), vcov. = sandwich))

}
\keyword{models}
