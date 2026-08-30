
### fixed effects models
### methods and tests contributed by Grant McDermott (@grantmcdermott)

library("multcomp", quietly = TRUE)

fixestOK <- require("fixest", quietly = TRUE)
if (fixestOK) {
    lmod  <- lm(Sepal.Length ~ Sepal.Width + as.factor(Species), iris)
    fmod  <- feols(Sepal.Length ~ Sepal.Width + as.factor(Species), iris)  
    fmod2 <- feols(Sepal.Length ~ Sepal.Width | Species, iris, vcov = "iid") ## see next model too
    fmod3 <- feols(Sepal.Length ~ Sepal.Width | Species, iris) ## default vcov is clustered by Species

    glmod  <- glht(lmod, "Sepal.Width==0")
    gfmod  <- glht(fmod, "Sepal.Width==0")
    gfmod2 <- glht(fmod2, "Sepal.Width==0")
    gfmod3 <- glht(fmod3, "Sepal.Width==0", vcov = "iid")
    stopifnot(all.equal(confint(glmod)$confint, confint(gfmod)$confint))
    stopifnot(all.equal(confint(glmod)$confint, confint(gfmod2)$confint))
    stopifnot(all.equal(confint(glmod)$confint, confint(gfmod3)$confint))
}
