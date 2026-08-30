\name{cld}
\alias{cld}
\alias{cld.glht}
\alias{cld.summary.glht}
\alias{cld.confint.glht}
\title{Set up a compact letter display of all pair-wise comparisons}
\description{
  Extract information from \code{glht}, \code{summary.glht} or
  \code{confint.glht} objects which is required to create
  and plot compact letter displays of all pair-wise comparisons.
}
\usage{
cld(object, ...)
\method{cld}{summary.glht}(object, level = 0.05, decreasing = FALSE, ...)
\method{cld}{glht}(object, level = 0.05, decreasing = FALSE, ...)
\method{cld}{confint.glht}(object, decreasing = FALSE, ...)
}
\arguments{
  \item{object}{
    An object of class \code{glht}, \code{summary.glht} or \code{confint.glht}.
  }
  \item{level}{
    Significance-level to be used to term a specific pair-wise
    comparison significant.
  }
  \item{decreasing}{
  logical. Should the order of the letters be increasing or decreasing?
  }
  \item{...}{additional arguments.}
}
\details{
  This function extracts all the information from \code{glht},
  \code{summary.glht} or \code{confint.glht} objects that is required
  to create a compact letter display of all pair-wise comparisons
  \bibcitep{multcomp::piepho2004}.
  In case the contrast matrix is not of type \code{"Tukey"}, an error
  is issued. In case of \code{confint.glht} objects, a pair-wise comparison
  is termed significant whenever a particular confidence interval contains 0.
  Otherwise, p-values are compared to the value of \code{"level"}.
  Once, this information is extracted, plotting of all pair-wise
  comparisons can be carried out.
}
\value{
  An object of class \code{cld}, a list with items:
  \item{y}{
    Values of the response variable of the original model.
  }
  \item{yname}{
    Name of the response variable.
  }
  \item{x}{
    Values of the variable used to compute Tukey contrasts.
  }
  \item{weights}{
    Weights used in the fitting process.
  }
  \item{lp}{
    Predictions from the fitted model.
  }
  \item{covar}{
    A logical indicating whether the fitted model contained covariates.
  }
  \item{signif}{
    Vector of logicals indicating significant differences with
    hyphenated names that identify pair-wise comparisons.
  }
}
\references{\bibshow{*}}
\seealso{
  \code{\link{glht}}
  \code{\link{plot.cld}}
}
\examples{
  ### multiple comparison procedures
  ### set up a one-way ANOVA
  data(warpbreaks)
  amod <- aov(breaks ~ tension, data = warpbreaks)
  ### specify all pair-wise comparisons among levels of variable "tension"
  tuk <- glht(amod, linfct = mcp(tension = "Tukey"))
  ### extract information
  tuk.cld <- cld(tuk)
  ### use sufficiently large upper margin
  old.par <- par(mai=c(1,1,1.25,1), no.readonly = TRUE)
  ### plot
  plot(tuk.cld)
  par(old.par)
  
  ### now using covariates
  data(warpbreaks)
  amod2 <- aov(breaks ~ tension + wool, data = warpbreaks)
  ### specify all pair-wise comparisons among levels of variable "tension"
  tuk2 <- glht(amod2, linfct = mcp(tension = "Tukey"))
  ### extract information
  tuk.cld2 <- cld(tuk2)
  ### use sufficiently large upper margin
  old.par <- par(mai=c(1,1,1.25,1), no.readonly = TRUE)
  ### plot using different colors
  plot(tuk.cld2, col=c("black", "red", "blue"))
  par(old.par)

  ### set up all pair-wise comparisons for count data
  data(Titanic)
  mod <- glm(Survived ~ Class, data = as.data.frame(Titanic), weights = Freq, family = binomial())
  ### specify all pair-wise comparisons among levels of variable "Class"
  glht.mod <- glht(mod, mcp(Class = "Tukey"))
  ### extract information
  mod.cld <- cld(glht.mod)
  ### use sufficiently large upper margin
  old.par <- par(mai=c(1,1,1.5,1), no.readonly = TRUE)
  ### plot
  plot(mod.cld)
  par(old.par)
}
