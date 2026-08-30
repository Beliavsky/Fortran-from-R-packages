\name{recovery}
\alias{recovery}
\docType{data}
\title{ Recovery Time Data Set }
\usage{data("recovery")}
\description{
  Recovery time after surgery.
}
\format{
  This data frame contains the following variables
  \describe{
    \item{blanket}{blanket type, a factor at four levels: \code{b0}, \code{b1},
                   \code{b2}, and \code{b3}.}
    \item{minutes}{response variable: recovery time after a surgical 
                   procedure.}
  }
}
\details{

  A company developed specialized heating blankets designed to help the
  body heat following a surgical procedure. Four types of blankets
  were tried on surgical patients with the aim of comparing the
  recovery time of patients. One of the blanket was a standard blanket
  that had been in use already in various hospitals. Data were published by
  \bibcitet{multcomp::westfall1999}, page 66.

}
\references{\bibshow{*}}
\examples{

  ### set up one-way ANOVA
  amod <- aov(minutes ~ blanket, data = recovery)

  ### set up multiple comparisons: one-sided Dunnett contrasts
  rht <- glht(amod, linfct = mcp(blanket = "Dunnett"), 
              alternative = "less")

  ### cf. Westfall et al. (1999, p. 80)
  confint(rht, level = 0.9)

  ### the same
  rht <- glht(amod, linfct = mcp(blanket = c("b1 - b0 >= 0", 
                                             "b2 - b0 >= 0", 
                                             "b3 - b0 >= 0")))
  confint(rht, level = 0.9)

}
\keyword{datasets}
