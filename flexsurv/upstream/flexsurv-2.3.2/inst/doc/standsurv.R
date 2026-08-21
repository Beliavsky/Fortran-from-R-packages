## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  fig.dim = c(8, 6),
  collapse = TRUE,
  comment = "#>"
)

## ----message = FALSE, warning = FALSE-----------------------------------------
library(flexsurv)
library(flexsurvcure)
library(ggplot2)
library(dplyr)
library(survminer)

## -----------------------------------------------------------------------------
data(bc)
set.seed(236236)
## Age at diagnosis is correlated with survival time. A longer survival time 
## gives a younger mean age
bc$age <- rnorm(dim(bc)[1], mean = 65 - scale(bc$recyrs, scale=F), sd = 5)
## Create age at diagnosis in days - used later for matching to expected rates
bc$agedays <- floor(bc$age * 365.25)
## Create some random diagnosis dates between 01/01/1984 and 31/12/1989
bc$diag <- as.Date(floor(runif(dim(bc)[1], as.Date("01/01/1984", "%d/%m/%Y"), 
                               as.Date("31/12/1989", "%d/%m/%Y"))), 
                   origin="1970-01-01")
## Create sex (assume all are female)
bc$sex <- factor("female")
## 2-level prognostic variable
bc$group2 <- ifelse(bc$group=="Good", "Good", "Medium/Poor")
head(bc)

## -----------------------------------------------------------------------------
km <- survfit(Surv(recyrs, censrec)~group2, data=bc)
kmsurvplot <- ggsurvplot(km) 
kmsurvplot + xlab("Time from diagnosis (years)")

## -----------------------------------------------------------------------------
model.weibull.sep <- flexsurvreg(Surv(recyrs, censrec)~group2, 
                                 anc = list(shape = ~ group2), 
                                 data=bc, dist="weibullPH")
model.weibull.sep

## -----------------------------------------------------------------------------
predictions <- summary(model.weibull.sep, type = "survival", tidy=T)
ggplot() + geom_line(aes(x=time, y=est, color = group2), data=predictions) + 
  geom_step(aes(x=time, y=surv, group=strata), data=kmsurvplot$data.survplot)

## -----------------------------------------------------------------------------
ss.weibull.sep.surv <- standsurv(model.weibull.sep,
                                             type = "survival",
                                             at = list(list(group2 = "Good"),
                                                       list(group2 = "Medium/Poor")))
ss.weibull.sep.surv

## -----------------------------------------------------------------------------
plot(ss.weibull.sep.surv)

## -----------------------------------------------------------------------------
plot(ss.weibull.sep.surv) + xlab("Time since diagnosis (years)") +
  geom_step(aes(x=time, y=surv, group=strata), data=kmsurvplot$data.survplot)

## -----------------------------------------------------------------------------
ss.weibull.sep.haz <- standsurv(model.weibull.sep,
                                            type = "hazard",
                                            at = list(list(group2 = "Good"),
                                                      list(group2 = "Medium/Poor")))
plot(ss.weibull.sep.haz) + xlab("Time since diagnosis (years)")

## -----------------------------------------------------------------------------
ss.weibull.sep.rmst <- standsurv(model.weibull.sep,
                                             type = "rmst",
                                             at = list(list(group2 = "Good"),
                                                       list(group2 = "Medium/Poor")))
plot(ss.weibull.sep.rmst) + xlab("Time since diagnosis (years)")

## -----------------------------------------------------------------------------
ss.weibull.sep.3 <- standsurv(model.weibull.sep,
                                          type = "survival",
                                                 at = list(list(group2 = "Good"),
                                                           list(group2 = "Medium/Poor")),
                                          contrast = "difference")
plot(ss.weibull.sep.3, contrast=TRUE) + xlab("Time since diagnosis (years)") +
  ylab("Difference in survival probabilities")  + geom_hline(yintercept = 0)

## -----------------------------------------------------------------------------
ss.weibull.sep.4 <- standsurv(model.weibull.sep,
                                          type = "hazard",
                                                 at = list(list(group2 = "Good"),
                                                           list(group2 = "Medium/Poor")),
                                          contrast = "ratio")
plot(ss.weibull.sep.4, contrast=TRUE) + xlab("Time since diagnosis (years)") +
  ylab("Hazard ratio") + geom_hline(yintercept = 1)

## -----------------------------------------------------------------------------
ss.weibull.sep.boot <- standsurv(model.weibull.sep,
                                          type = "survival",
                                                 at = list(list(group2 = "Good"),
                                                           list(group2 = "Medium/Poor")),
                                          t = seq(0,7,length=10),
                                          ci = TRUE,
                                          boot = TRUE,
                                          B = 100,
                                          seed = 2367)

ss.weibull.sep.deltam <- standsurv(model.weibull.sep,
                                          type = "survival",
                                                 at = list(list(group2 = "Good"),
                                                           list(group2 = "Medium/Poor")),
                                          t = seq(0,7,length=10),
                                          ci = TRUE,
                                          boot = FALSE)

plot(ss.weibull.sep.boot, ci = TRUE) +  
  geom_ribbon(aes(x=time, ymin=survival_lci, ymax=survival_uci, color=at, linetype = "Delta method"), fill=NA,
              data=attr(ss.weibull.sep.deltam,"standpred_at")) +
  scale_linetype_manual(values = c("Bootstrap" = "solid", "Delta method"= "dashed")) +
  ggtitle("Comparison of bootstrap and delta method confidence intervals")


## -----------------------------------------------------------------------------
model.weibull.age.sep <- flexsurvreg(Surv(recyrs, censrec)~group2 + age, 
                                 anc = list(shape = ~ group2 + age), 
                                 data=bc, dist="weibullPH")

## Marginal survival standardized to age distribution of study population
ss.weibull.age.sep.surv <- standsurv(model.weibull.age.sep,
                                             type = "survival",
                                             at = list(list(group2 = "Good"),
                                                       list(group2 = "Medium/Poor")),
                                             t = seq(0,7,length=50)
)

## Marginal survival standardized to an older population
# create a new prediction dataset as a copy of the bc data but whose ages are drawn from
# a normal distribution with mean age 75, sd 5.
newpred.data <- bc
set.seed(247)
newpred.data$age = rnorm(dim(bc)[1], 75, 5)
ss.weibull.age2.sep.surv <- standsurv(model.weibull.age.sep,
                                             type = "survival",
                                             at = list(list(group2 = "Good"),
                                                       list(group2 = "Medium/Poor")),
                                             t = seq(0,7,length=50),
                                             newdata=newpred.data)

## Overlay both marginal survival curves
plot(ss.weibull.age.sep.surv) + 
  geom_line(aes(x=time, y=survival, color=at, linetype = "Older population"),
            data = attr(ss.weibull.age2.sep.surv, "standpred_at") ) +
  scale_linetype_manual(values = c("Study" = "solid", "Older population"= "dashed"))


## -----------------------------------------------------------------------------
summary(survexp.us)

## -----------------------------------------------------------------------------
ss.weibull.sep.expected <- standsurv(model.weibull.sep,
                                                 type = "survival",
                                                 at = list(list(group2 = "Good"),
                                                           list(group2 = "Medium/Poor")),
                                                 t = seq(0,7,length=50),
                                                 rmap=list(sex = sex,
                                                           year = diag,
                                                           age = agedays
                                                 ),
                                                 ratetable = survexp.us,
                                                 scale.ratetable = 365.25,
                                                 newdata = bc
)
plot(ss.weibull.sep.expected, expected = T)

## -----------------------------------------------------------------------------
ss.weibull.sep.expectedh <- standsurv(model.weibull.sep,
                                                 type = "hazard",
                                                 at = list(list(group2 = "Good"),
                                                           list(group2 = "Medium/Poor")),
                                                 t = seq(0,7,length=50),
                                                 rmap=list(sex = sex,
                                                           year = diag,
                                                           age = agedays
                                                 ),
                                                 ratetable = survexp.us,
                                                 scale.ratetable = 365.25,
                                                 newdata = bc
)
plot(ss.weibull.sep.expectedh, expected = T)

## -----------------------------------------------------------------------------
## reshape US lifetable to be a tidy data.frame, and convert rates to per person-year as flexsurv regression is in years
survexp.us.df <- as.data.frame.table(survexp.us, responseName = "exprate") %>%
  mutate(exprate = 365.25 * exprate)
survexp.us.df$age <- as.numeric(as.character(survexp.us.df$age))
survexp.us.df$year <- as.numeric(as.character(survexp.us.df$year))

## Obtain attained age and attained calendar year in (whole) years
bc <- bc %>% mutate(attained.age.yr = floor(age + recyrs),
                    attained.year = lubridate::year(diag + rectime))

## merge in (left join) expected rates at event time
bc <- bc %>% left_join(survexp.us.df, by = c("attained.age.yr"="age", 
                                        "attained.year"="year", 
                                        "sex"="sex")) 

# A stratified relative survival mixture-cure model
model.weibull.sep.rs <- flexsurvcure(Surv(recyrs, censrec)~group2, 
                                 anc = list(shape = ~ group2,
                                            scale = ~ group2), 
                                 data=bc, dist="weibullPH",
                                 bhazard=exprate)

model.weibull.sep.rs

## -----------------------------------------------------------------------------
## All-cause survival
ss.weibull.sep.rs.surv <- standsurv(model.weibull.sep.rs,
                                                 type = "survival",
                                                 at = list(list(group2 = "Good"),
                                                           list(group2 = "Medium/Poor")),
                                                 t = seq(0,30,length=50),
                                                 rmap=list(sex = sex,
                                                           year = diag,
                                                           age = agedays
                                                 ),
                                                 ratetable = survexp.us,
                                                 scale.ratetable = 365.25,
                                                 newdata = bc
)
plot(ss.weibull.sep.rs.surv, expected = T)

# All-cause hazard
ss.weibull.sep.rs.haz <- standsurv(model.weibull.sep.rs,
                                                 type = "hazard",
                                                 at = list(list(group2 = "Good"),
                                                           list(group2 = "Medium/Poor")),
                                                 t = seq(0,30,length=50),
                                                 rmap=list(sex = sex,
                                                           year = diag,
                                                           age = agedays
                                                 ),
                                                 ratetable = survexp.us,
                                                 scale.ratetable = 365.25,
                                                 newdata = bc
)
plot(ss.weibull.sep.rs.haz, expected = T)

## -----------------------------------------------------------------------------
# Excess hazard
ss.weibull.sep.rs.excesshaz <- standsurv(model.weibull.sep.rs,
                                                 type = "excesshazard",
                                                 at = list(list(group2 = "Good"),
                                                           list(group2 = "Medium/Poor")),
                                                 t = seq(0,30,length=50),
                                                 rmap=list(sex = sex,
                                                           year = diag,
                                                           age = agedays
                                                 ),
                                                 ratetable = survexp.us,
                                                 scale.ratetable = 365.25,
                                                 newdata = bc
)
plot(ss.weibull.sep.rs.excesshaz)

