\name{plot.cld}
\alias{plot.cld}
\title{Plot a cld object}
\description{
  Plot information of \code{glht}, \code{summary.glht} or \code{confint.glht}
  objects stored as \code{cld} objects together with a compact
  letter display of all pair-wise comparisons.
}
\usage{
\method{plot}{cld}(x, type = c("response", "lp"), ...)
}
\arguments{
  \item{x}{
    An object of class \code{cld}.
  }
  \item{type}{
    Should the response or the linear predictor (lp) be plotted.
    If there are any covariates, the lp is automatically used. To
    use the response variable, set \code{type="response"} and \code{covar=FALSE}
    of the \code{cld} object.
  }
  \item{...}{
    Other optional print parameters which are passed to the plotting functions.
  }
}
\details{
  This function plots the information stored in \code{glht}, \code{summary.glht} or
  \code{confint.glht} objects. Prior to plotting, these objects have to be converted to
  \code{cld} objects (see \code{\link{cld}} for details).
  All types of plots include a compact letter display (cld) of all pair-wise comparisons.
  Equal letters indicate no significant differences. Two levels are significantly
  different, in case they do not have any letters in common.
  If the fitted model contains any covariates, a boxplot of the linear predictor is
  generated with the cld within the upper margin. Otherwise, three different types
  of plots are used depending on the class of variable \code{y} of the \code{cld} object.
  In case of \code{class(y) == "numeric"}, a boxplot is generated using the response variable,
  classified according to the levels of the variable used for the Tukey contrast
  matrix. Is \code{class(y) == "factor"}, a mosaic plot is generated, and the cld is printed
  above. In case of \code{class(y) == "Surv"}, a plot of fitted survival functions is generated
  where the cld is plotted within the legend.
  The compact letter display is computed using the algorithm of Piepho (2004).
  Note: The user has to provide a sufficiently large upper margin which can be used to
        depict the compact letter display (see examples).
}
\references{
  Hans-Peter Piepho (2004), An Algorithm for a Letter-Based Representation of All-Pairwise Comparisons,
  \emph{Journal of Computational and Graphical Statistics}, \bold{13}(2), 456--466.
}
\seealso{
  \code{\link{glht}}
  \code{\link{cld}}
  \code{\link{cld.summary.glht}}
  \code{\link{cld.confint.glht}}
  \code{\link{cld.glht}}
  \code{\link{boxplot}}
  \code{\link{mosaicplot}}
  \code{\link[survival]{plot.survfit}}
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
  old.par <- par(mai=c(1,1,1.25,1), no.readonly=TRUE)
  ### plot
  plot(tuk.cld)
  par(old.par)

  ### now using covariates
  amod2 <- aov(breaks ~ tension + wool, data = warpbreaks)
  tuk2 <- glht(amod2, linfct = mcp(tension = "Tukey"))
  tuk.cld2 <- cld(tuk2)
  old.par <- par(mai=c(1,1,1.25,1), no.readonly=TRUE)
  ### use different colors for boxes
  plot(tuk.cld2, col=c("green", "red", "blue"))
  par(old.par)
  
  ### get confidence intervals
  ci.glht <- confint(tuk)
  ### plot them
  plot(ci.glht)
  old.par <- par(mai=c(1,1,1.25,1), no.readonly=TRUE)
  ### use 'confint.glht' object to plot all pair-wise comparisons
  plot(cld(ci.glht), col=c("white", "blue", "green"))
  par(old.par)
  
  ### set up all pair-wise comparisons for count data
  data(Titanic)
  mod <- glm(Survived ~ Class, data = as.data.frame(Titanic), 
             weights = Freq, family = binomial())
  ### specify all pair-wise comparisons among levels of variable "Class"
  glht.mod <- glht(mod, mcp(Class = "Tukey"))
  ### extract information
  mod.cld <- cld(glht.mod)
  ### use sufficiently large upper margin
  old.par <- par(mai=c(1,1,1.5,1), no.readonly=TRUE)
  ### plot
  plot(mod.cld)
  par(old.par)
  
  ### set up all pair-wise comparisons of a Cox-model
  if (require("survival") && require("MASS")) {
    ### construct 4 classes of age
    Melanoma$Cage <- factor(sapply(Melanoma$age, function(x){
                            if( x <= 25 ) return(1)
                            if( x > 25 & x <= 50 ) return(2)
                            if( x > 50 & x <= 75 ) return(3)
                            if( x > 75 & x <= 100) return(4) }
                           ))
    ### fit Cox-model
    cm <- coxph(Surv(time, status == 1) ~ Cage, data = Melanoma)
    ### specify all pair-wise comparisons among levels of "Cage"
    cm.glht <- glht(cm, mcp(Cage = "Tukey"))
    # extract information & plot
    old.par <- par(no.readonly=TRUE)
    ### use mono font family
    if (dev.interactive())
        old.par <- par(family = "mono")
    plot(cld(cm.glht), col=c("black", "red", "blue", "green"))
    par(old.par)
  }

  if (require("nlme") && require("lme4")) {
    data("ergoStool", package = "nlme")

    stool.lmer <- lmer(effort ~ Type + (1 | Subject),
                       data = ergoStool)
    glme41 <- glht(stool.lmer, mcp(Type = "Tukey"))

    old.par <- par(mai=c(1,1,1.5,1), no.readonly=TRUE)
    plot(cld(glme41))
    par(old.par)
  }
}
\keyword{hplot}
