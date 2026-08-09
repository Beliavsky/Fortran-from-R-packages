## ----setup, echo=FALSE--------------------------------------------------------
rm(list=ls()) # clear workspace for each major section

## ----nlsbtsimple, echo=TRUE---------------------------------------------------
# Bounds Test nlsbtsimple.R (see BT.RES in Nash and Walker-Smith (1987))
rm(list=ls())
bt.res<-function(x){ x }
bt.jac<-function(x){nn <- length(x); JJ <- diag(nn); attr(JJ, "gradient") <- JJ;  JJ}
n <- 4
x<-rep(0,n) ; lower<-rep(NA,n); upper<-lower # to get arrays set
for (i in 1:n) { lower[i]<-1.0*(i-1)*(n-1)/n; upper[i]<-1.0*i*(n+1)/n }
x <-0.5*(lower+upper) # start on mean
xnames<-as.character(1:n) # name our parameters just in case
for (i in 1:n){ xnames[i]<-paste("p",xnames[i],sep='') }
names(x) <- xnames
require(minpack.lm)
require(nlsr)
bsnlf0<-nlfb(start=x, resfn=bt.res, jacfn=bt.jac) # unconstrained
bsnlf0
bsnlm0<-nls.lm(par=x, fn=bt.res, jac=bt.jac) # unconstrained
bsnlm0
bsnlf1<-nlfb(start=x, resfn=bt.res, jacfn=bt.jac, lower=lower, upper=upper)
bsnlf1
bsnlm1<-nls.lm(par=x, fn=bt.res, jac=bt.jac, lower=lower, upper=upper)
bsnlm1

# single value bounds
bsnlf2<-nlfb(start=x, resfn=bt.res, jacfn=bt.jac, lower=0.25, upper=4)
bsnlf2
# nls.lm will NOT expand single value bounds
bsnlm2<-try(nls.lm(par=x, fn=bt.res, jac=bt.jac, lower=0.25, upper=4))

## ----nlsLMoob, echo=TRUE------------------------------------------------------
x<-rep(0,n) # resetting to this puts x out of bounds
cat("lower:"); print(lower)
cat("upper:"); print(upper)
names(x) <- xnames # to ensure names set
obnlm<-nls.lm(par=x, fn=bt.res, jac=bt.jac, lower=lower, upper=upper)
obnlm

## ----Hobbssetup, echo=TRUE----------------------------------------------------
## Use the Hobbs Weed problem
weed <- c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443,
          38.558, 50.156, 62.948, 75.995, 91.972)
tt <- 1:12
weeddf <- data.frame(y=weed, tt=tt)
st <- c(b1=1, b2=1, b3=1) # a default starting vector (named!)
## Unscaled model
wmodu <- y ~ b1/(1+b2*exp(-b3*tt))
## Scaled model
wmods <- y ~ 100*b1/(1+10*b2*exp(-0.1*b3*tt))
## We can provide the residual and Jacobian as functions
# Unscaled Hobbs problem
hobbs.res  <-  function(x){ # scaled Hobbs weeds problem -- residual
  # This variant uses looping
  if(length(x) != 3) stop("hobbs.res -- parameter vector n!=3")
  y  <-  c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443, 
           38.558, 50.156, 62.948, 75.995, 91.972)
  tt  <-  1:12
  res  <-  x[1]/(1+x[2]*exp(-x[3]*tt)) - y
}

hobbs.jac  <-  function(x) { # scaled Hobbs weeds problem -- Jacobian
  jj  <-  matrix(0.0, 12, 3)
  tt  <-  1:12
  yy  <-  exp(-x[3]*tt)
  zz  <-  1.0/(1+x[2]*yy)
  jj[tt,1]   <-   zz
  jj[tt,2]   <-   -x[1]*zz*zz*yy
  jj[tt,3]   <-   x[1]*zz*zz*yy*x[2]*tt
  attr(jj, "gradient") <- jj
  jj
}
# Scaled Hobbs problem
shobbs.res  <-  function(x){ # scaled Hobbs weeds problem -- residual
  # This variant uses looping
  if(length(x) != 3) stop("shobbs.res -- parameter vector n!=3")
  y  <-  c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443, 
           38.558, 50.156, 62.948, 75.995, 91.972)
  tt  <-  1:12
  res  <-  100.0*x[1]/(1+x[2]*10.*exp(-0.1*x[3]*tt)) - y
}

shobbs.jac  <-  function(x) { # scaled Hobbs weeds problem -- Jacobian
  jj  <-  matrix(0.0, 12, 3)
  tt  <-  1:12
  yy  <-  exp(-0.1*x[3]*tt)
  zz  <-  100.0/(1+10.*x[2]*yy)
  jj[tt,1]   <-   zz
  jj[tt,2]   <-   -0.1*x[1]*zz*zz*yy
  jj[tt,3]   <-   0.01*x[1]*zz*zz*yy*x[2]*tt
  attr(jj, "gradient") <- jj
  jj
}

## ----exnlxb-------------------------------------------------------------------
require(nlsr)
# Use Hobbs scaled formula model
anlxb1  <-  try(nlxb(wmods, start=st, data=weeddf))
print(anlxb1)
# summary(anlxb1) # for different output
# Test without starting parameters
anlxb1nostart  <-  try(nlxb(wmodu, data=weeddf))
print(anlxb1nostart)

## ----nlxb1side1---------------------------------------------------------------
# One-sided unscaled Hobbs weed model formula
wmodu1  <-  ~ b1/(1+b2*exp(-b3*tt)) - y
anlxb11  <-  try(nlxb(wmodu1, start=st, data=weeddf))
print(anlxb11)

## ----c003---------------------------------------------------------------------
# st  <-  c(b1=1, b2=1, b3=1)
wrj <- model2rjfun(wmods, st, data=weeddf)
wrj
weedux <- nlxb(wmods, start=st, data=weeddf) 
print(weedux)

wss <- model2ssgrfun(wmods, st, data=weeddf)
print(wss)

# We can get expressions used to calculate these as follows:
wexpr.rj <- modelexpr(wrj) 
print(wexpr.rj)

wexpr.ss <- modelexpr(wss) 
print(wexpr.ss)

## ----c004---------------------------------------------------------------------
require(nlsr)

cat("try nlfb\n")
st  <-  c(b1=1, b2=1, b3=1)
ans1 <- nlfb(st, shobbs.res, shobbs.jac, trace=FALSE)
summary(ans1)

## No jacobian function -- use internal approximation
ans1n <- nlfb(st, shobbs.res, trace=FALSE, control=list(japprox="jafwd", watch=FALSE)) # NO jacfn -- tell it fwd approx
print(ans1n)
## difference
coef(ans1)-coef(ans1n)

## ----c001---------------------------------------------------------------------
coef(anlxb1) # this is solution of scaled problem from unit start

## ----c006---------------------------------------------------------------------
### From examples above
print(weedux)
print(ans1)

## ----c006a--------------------------------------------------------------------
### From examples above
summary(weedux)
summary(ans1)

## ----c009, echo=TRUE----------------------------------------------------------
## The following attempt at Hobbs unscaled with nls() fails!
rm(list=ls()) # clear before starting
weed <- c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443,
          38.558, 50.156, 62.948, 75.995, 91.972)
tt <- 1:12
weeddf <- data.frame(tt, weed)
st <- c(b1=1, b2=1, b3=1) # a default starting vector (named!)
wmodu<- weed ~ b1/(1 + b2 * exp(- b3 * tt))
anls1u  <-  try(nls(wmodu, start=st, trace=FALSE, data=weeddf)) # fails
## But we succeed by calling nlxb first.
library(nlsr)
anlxb1un  <-  try(nlxb(wmodu, start=st, trace=FALSE, data=weeddf))
print(anlxb1un)
st2 <- coef(anlxb1un) # We try to start nls from solution using nlxb
anls2u  <-  try(nls(wmodu, start=st2, trace=FALSE, data=weeddf))
print(anls2u)
## Or we can simply call wrapnlsr FROM THE ORIGINAL start
anls2a  <-  try(wrapnlsr(wmodu, start=st, trace=FALSE, data=weeddf))
summary(anls2a)
# # For comparison could run nlsLM
# require(minpack.lm)
# anlsLM1u  <-  try(nlsLM(wmodu, start=st, trace=FALSE, data=weeddf))
# print(anlsLM1u)

## ----c007, echo=TRUE----------------------------------------------------------
## Use shobbs example -- Scaled Hobbs problem
shobbs.res  <-  function(x){ # scaled Hobbs weeds problem -- residual
  if(length(x) != 3) stop("shobbs.res -- parameter vector n!=3")
  y  <-  c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443, 
           38.558, 50.156, 62.948, 75.995, 91.972)
  tt  <-  1:12
  res  <-  100.0*x[1]/(1+x[2]*10.*exp(-0.1*x[3]*tt)) - y
}

shobbs.jac  <-  function(x) { # scaled Hobbs weeds problem -- Jacobian
  jj  <-  matrix(0.0, 12, 3)
  tt  <-  1:12
  yy  <-  exp(-0.1*x[3]*tt)
  zz  <-  100.0/(1+10.*x[2]*yy)
  jj[tt,1]   <-   zz
  jj[tt,2]   <-   -0.1*x[1]*zz*zz*yy
  jj[tt,3]   <-   0.01*x[1]*zz*zz*yy*x[2]*tt
  attr(jj, "gradient") <- jj
  jj
}
RG <- resgr(st, shobbs.res, shobbs.jac)
RG
SS <- resss(st, shobbs.res)
SS

## ----c002---------------------------------------------------------------------
newDeriv() # a call with no arguments returns a list of available functions
           # for which derivatives are currently defined
newDeriv(sin(x)) # a call with a function that is in the list of available derivatives
                 # returns the derivative expression for that function
nlsDeriv(~ sin(x+y), "x") # partial derivative of this function w.r.t. "x"

## CAUTION !! ##
newDeriv(joe(x)) # but an undefined function returns NULL
newDeriv(joe(x), deriv=log(x^2)) # We can define derivatives, though joe() is meaningless.
nlsDeriv(~ joe(x+z), "x")

# Some examples of usage
f <- function(x) x^2
newDeriv(f(x), 2*x*D(x))
nlsDeriv(~ f(abs(x)), "x")
nlsDeriv(~ x^2, "x")
nlsDeriv(~ (abs(x)^2), "x")

# derivatives of distribution functions
nlsDeriv(~ pnorm(x, sd=2, log = TRUE), "x")
# get evaluation code from a formula
codeDeriv(~ pnorm(x, sd = sd, log = TRUE), "x")
# wrap it in a function call
fnDeriv(~ pnorm(x, sd = sd, log = TRUE), "x")
f <- fnDeriv(~ pnorm(x, sd = sd, log = TRUE), "x", args = alist(x =, sd = 2))
f
f(1)
100*(f(1.01) - f(1))  # Should be close to the gradient
                      # The attached gradient attribute (from f(1.01)) is
                      # meaningless after the subtraction.

## ----c008---------------------------------------------------------------------
## nlsSimplify
nlsSimplify(quote(a + 0))
nlsSimplify(quote(exp(1)), verbose = TRUE)

nlsSimplify(quote(sqrt(a + b)))  # standard rule
## sysSimplifications
# creates a new environment whose parent is emptyenv()  Why?
str(sysSimplifications)
myrules <- new.env(parent = sysSimplifications)
## newSimplification
newSimplification(sqrt(a), TRUE, a^0.5, simpEnv = myrules)
nlsSimplify(quote(sqrt(a + b)), simpEnv = myrules)

## isFALSE
print(isFALSE(1==2))
print(isFALSE(2==2))

## isZERO
print(isZERO(0))
x <- -0
print(isZERO(x))
x <- 0
print(isZERO(x))
print(isZERO(~(-1)))
print(isZERO("-1"))
print(isZERO(expression(-1)))


## isONE
print(isONE(1))
x <- 1
print(isONE(x))
print(isONE(~(1)))
print(isONE("1"))
print(isONE(expression(1)))


## isMINUSONE
print(isMINUSONE(-1))
x <- -1
print(isMINUSONE(x))
print(isMINUSONE(~(-1)))
print(isMINUSONE("-1"))
print(isMINUSONE(expression(-1)))

## isCALL  ?? don't have good understanding of this
x <- -1
print(isCALL(x,"isMINUSONE"))
print(isCALL(x, quote(isMINUSONE)))
      
## findSubexprs
findSubexprs(expression(x^2, x-y, y^2-x^2))

## sysDerivs
# creates a new environment whose parent is emptyenv()  Why? 
# Where are derivative definitions are stored?
str(sysDerivs)

## ----nlfbnumex, echo=TRUE-----------------------------------------------------
## Test with functional spec. of problem
## Default call WITH jacobian function
ans1 <- nlfb(st, resfn=shobbs.res, jacfn=shobbs.jac)
ans1
## No jacobian function -- and no japprox control setting
ans1n <- try(nlfb(st, shobbs.res)) # NO jacfn
ans1n
## Force jafwd approximation
ans1nf <- nlfb(st, shobbs.res, control=list(japprox="jafwd")) # NO jacfn, but specify fwd approx
ans1nf
## Alternative specification
ans1nfa <- nlfb(st, shobbs.res, jacfn="jafwd") # NO jacfn, but specify fwd approx in jacfn
ans1nfa
## Coeff differences from analytic: 
ans1nf$coefficients-ans1$coefficients
## Force jacentral approximation
ans1nc <- nlfb(st, shobbs.res, control=list(japprox="jacentral")) # NO jacfn
ans1nc
## Coeff differences from analytic:
ans1nc$coefficients-ans1$coefficients
## Force jaback approximation
ans1nb <- nlfb(st, shobbs.res, control=list(japprox="jaback")) # NO jacfn
ans1nb
## Coeff differences from analytic:
ans1nb$coefficients-ans1$coefficients
## Force jand approximation
ans1nn <- nlfb(st, shobbs.res, control=list(japprox="jand"), trace=FALSE) # NO jacfn
ans1nn
## Coeff differences from analytic:
ans1nc$coefficients-ans1$coefficients

## ----nlxbnumex1, echo=TRUE----------------------------------------------------
## rm(list=ls()) # clear workspace for each major section
## Use the Hobbs Weed problem
weed <- c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443,
          38.558, 50.156, 62.948, 75.995, 91.972)
tt <- 1:12
weeddf <- data.frame(y=weed, tt=tt)
st1 <- c(b1=1, b2=1, b3=1) # a default starting vector (named!)
wmods <- y ~ 100*b1/(1+10*b2*exp(-0.1*b3*tt)) ## Scaled model
print(wmods)
anlxb1  <-  try(nlxb(wmods, start=st1, data=weeddf))
print(anlxb1)

anlxb1fwd  <-  (nlxb(wmods, start=st1, data=weeddf,
    control=list(japprox="jafwd")))
print(anlxb1fwd)

## fwd - analytic
print(anlxb1fwd$coefficients - anlxb1$coefficients)

anlxb1bak  <- (nlxb(wmods, start=st1, data=weeddf,
    control=list(japprox="jaback")))
print(anlxb1bak)
## back - analytic
print(anlxb1bak$coefficients - anlxb1$coefficients)

anlxb1cen  <- (nlxb(wmods, start=st1, data=weeddf,
    control=list(japprox="jacentral")))
print(anlxb1cen)
## central - analytic
print(anlxb1cen$coefficients - anlxb1$coefficients)

anlxb1nd  <-  (nlxb(wmods, start=st1, data=weeddf,
    control=list(japprox="jand")))
print(anlxb1nd)
## numDeriv - analytic
print(anlxb1nd$coefficients - anlxb1$coefficients)

## ----wts1, echo=TRUE----------------------------------------------------------
## weighted nonlinear regression using inverse squared variance of the response
require(minpack.lm)
Treated <- Puromycin[Puromycin$state == "treated", ]
# We want the variance of each "group" of the rate variable 
# for which the conc variable is the same. We first find
# this variance by group using the tapply() function.
var.Treated <- tapply(Treated$rate, Treated$conc, var)
var.Treated <- rep(var.Treated, each = 2)
Pur.wt1nls <- nls(rate ~ (Vm * conc)/(K + conc), data = Treated, 
               start = list(Vm = 200, K = 0.1), weights = 1/var.Treated^2)
Pur.wt1nlm <- nlsLM(rate ~ (Vm * conc)/(K + conc), data = Treated, 
               start = list(Vm = 200, K = 0.1), weights = 1/var.Treated^2)
Pur.wt1nlx <- nlxb(rate ~ (Vm * conc)/(K + conc), data = Treated, 
               start = list(Vm = 200, K = 0.1), weights = 1/var.Treated^2)
rnls <- summary(Pur.wt1nls)$residuals
ssrnls<-as.numeric(crossprod(rnls))
rnlm <- summary(Pur.wt1nlm)$residuals
ssrnlm<-as.numeric(crossprod(rnlm))
rnlx <- Pur.wt1nlx$resid
ssrnlx<-as.numeric(crossprod(rnlx))
cat("Compare nls and nlsLM: ", all.equal(coef(Pur.wt1nls), coef(Pur.wt1nlm)),"\n")
cat("Compare nls and nlsLM: ", all.equal(coef(Pur.wt1nls), coef(Pur.wt1nlx)),"\n")
cat("Sumsquares nls - nlsLM: ", ssrnls-ssrnlm,"\n")
cat("Sumsquares nls - nlxb: ", ssrnls-ssrnlx,"\n")

## ----wtsfunction, eval=FALSE--------------------------------------------------
# ## minpack.lm::wfct
# ### Examples from 'nls' doc where wfct used ###
# ## Weighting by inverse of response 1/y_i:
# wtt1nlm<-nlsLM(rate ~ Vm * conc/(K + conc), data = Treated,
#               start = c(Vm = 200, K = 0.05), weights = wfct(1/rate))
# print(wtt1nlm)
# 
# wtt1nls<-nls(rate ~ Vm * conc/(K + conc), data = Treated,
#           start = c(Vm = 200, K = 0.05), weights = wfct(1/rate))
# print(wtt1nls)
# 
# wtt1nlx<-nlxb(rate ~ Vm * conc/(K + conc), data = Treated,
#           start = c(Vm = 200, K = 0.05), weights = wfct(1/rate))
# print(wtt1nlx)
# 
# ## Weighting by square root of predictor \sqrt{x_i}:
# # ?? why does try() not work
# wtt2nlm<-nlsLM(rate ~ Vm * conc/(K + conc), data = Treated,
#        start = c(Vm = 200, K = 0.05), weights = wfct(sqrt(conc)))
# print(wtt2nlm)
# 
# wtt2nls<-nls(rate ~ Vm * conc/(K + conc), data = Treated,
#        start = c(Vm = 200, K = 0.05), weights = wfct(sqrt(conc)))
# print(wtt2nls)
# 
# wtt2nlx<-nlxb(rate ~ Vm * conc/(K + conc), data = Treated,
#        start = c(Vm = 200, K = 0.05), weights = wfct(sqrt(conc)))
# print(wtt2nlx)
# 
# ## Weighting by inverse square of fitted values 1/\hat{y_i}^2:
# wtt3nlm<-nlsLM(rate ~ Vm * conc/(K + conc), data = Treated,
#        start = c(Vm = 200, K = 0.05), weights = wfct(1/fitted^2))
# print(wtt3nlm)
# 
# wtt3nls<-nls(rate ~ Vm * conc/(K + conc), data = Treated,
#        start = c(Vm = 200, K = 0.05), weights = wfct(1/fitted^2))
# print(wtt3nls)
# 
# # These don't work -- there is no fitted() function available in nlsr
# # wtt3nlx<-try(nlxb(rate ~ Vm * conc/(K + conc), data = Treated,
# #        start = c(Vm = 200, K = 0.05), weights = wfct(1/fitted^2)))
# ##  (nlxb(rate ~ Vm * conc/(K + conc), data = Treated,
# ##        start = c(Vm = 200, K = 0.05), weights = wfct(1/(resid+rate)^2)))
# ## (wrapnlsr(rate ~ Vm * conc/(K + conc), data = Treated,
# ##        start = c(Vm = 200, K = 0.05), weights = wfct(1/(resid+rate)^2)))
# 
# ## Weighting by inverse variance 1/\sigma{y_i}^2:
# wtt4nlm<-nlsLM(rate ~ Vm * conc/(K + conc), data = Treated,
#        start = c(Vm = 200, K = 0.05), weights = wfct(1/error^2))
# print(wtt4nlm)
# 
# wtt4nls<-nls(rate ~ Vm * conc/(K + conc), data = Treated,
#        start = c(Vm = 200, K = 0.05), weights = wfct(1/error^2))
# print(wtt4nls)
# 
# wtt4nlx<-nlxb(rate ~ Vm * conc/(K + conc), data = Treated,
#        start = c(Vm = 200, K = 0.05), weights = wfct(1/error^2))
# print(wtt4nlx)
# 
# wtt5nls<-nls(rate ~ Vm * conc/(K + conc), data = Treated,
#               start = c(Vm = 200, K = 0.05), weights = wfct(1/resid^2))
# print(wtt5nls)
# wtt5nlm<-nlsLM(rate ~ Vm * conc/(K + conc), data = Treated,
#               start = c(Vm = 200, K = 0.05), weights = wfct(1/resid^2))
# print(wtt5nlm)
# ## Won't work! No resid function available from nlxb.
# ## wtt5nlx<-nlxb(rate ~ Vm * conc/(K + conc), data = Treated,
# ##            start = c(Vm = 200, K = 0.05), weights = wfct(1/resid^2))
# ## print(wtt5nlx)

## ----nolhs, echo=TRUE---------------------------------------------------------
st<-c(b1=1, b2=1, b3=1)
frm <- ~ b1/(1+b2*exp(-b3*tt))-y
nolhs <- nlxb(formula=frm, data=weeddf, start=st)
nolhs

## ----wts1side-----------------------------------------------------------------
## weighted nonlinear regression using 1-sided model functions
Treated <- Puromycin[Puromycin$state == "treated", ]
weighted.MM <- function(resp, conc, Vm, K)
{
  ## Purpose: exactly as white book p. 451 -- RHS for nls()
  ##  Weighted version of Michaelis-Menten model
  ## ----------------------------------------------------------
  ## Arguments: 'y', 'x' and the two parameters (see book)
  ## ----------------------------------------------------------
  ## Author: Martin Maechler, Date: 23 Mar 2001
  
  pred <- (Vm * conc)/(K + conc)
  (resp - pred) / sqrt(pred)
}

Pur.wtMMnls <- nls( ~ weighted.MM(rate, conc, Vm, K), data = Treated,
               start = list(Vm = 200, K = 0.1))
print(Pur.wtMMnls)

Pur.wtMMnlm <- nlsLM( ~ weighted.MM(rate, conc, Vm, K), data = Treated,
                    start = list(Vm = 200, K = 0.1))
print(Pur.wtMMnlm)

Pur.wtMMnlx <- nlxb( ~ weighted.MM(rate, conc, Vm, K), data = Treated,
                    start = list(Vm = 200, K = 0.1), control=list(japprox="jafwd"))
print(Pur.wtMMnlx)

# Another approach
## Passing arguments using a list that can not be coerced to a data.frame
lisTreat <- with(Treated, list(conc1 = conc[1], conc.1 = conc[-1], rate = rate))
weighted.MM1 <- function(resp, conc1, conc.1, Vm, K){
  conc <- c(conc1, conc.1)
  pred <- (Vm * conc)/(K + conc)
  (resp - pred) / sqrt(pred)
}
Pur.wtMM1nls <- nls( ~ weighted.MM1(rate, conc1, conc.1, Vm, K),
                data = lisTreat, start = list(Vm = 200, K = 0.1))
print(Pur.wtMM1nls)
# Should we put in comparison of coeffs / sumsquares??
# stopifnot(all.equal(coef(Pur.wt), coef(Pur.wt1)))
Pur.wtMM1nlm <- nlsLM( ~ weighted.MM1(rate, conc1, conc.1, Vm, K),
                     data = lisTreat, start = list(Vm = 200, K = 0.1))
print(Pur.wtMM1nlm)
Pur.wtMM1nlx <- nlxb( ~ weighted.MM1(rate, conc1, conc.1, Vm, K),
                     data = lisTreat, start = list(Vm = 200, K = 0.1), control=list(japprox="jafwd"))
print(Pur.wtMM1nlx)
# yet another approach
## Chambers and Hastie (1992) Statistical Models in S  (p. 537):
## If the value of the right side [of formula] has an attribute called
## 'gradient' this should be a matrix with the number of rows equal
## to the length of the response and one column for each parameter.
weighted.MM.grad <- function(resp, conc1, conc.1, Vm, K) {
  conc <- c(conc1, conc.1)
  K.conc <- K+conc
  dy.dV <- conc/K.conc
  dy.dK <- -Vm*dy.dV/K.conc
  pred <- Vm*dy.dV
  pred.5 <- sqrt(pred)
  dev <- (resp - pred) / pred.5
  Ddev <- -0.5*(resp+pred)/(pred.5*pred)
  attr(dev, "gradient") <- Ddev * cbind(Vm = dy.dV, K = dy.dK)
  dev
}
Pur.wtMM.gradnls <- nls( ~ weighted.MM.grad(rate, conc1, conc.1, Vm, K),
                    data = lisTreat, start = list(Vm = 200, K = 0.1))
print(Pur.wtMM.gradnls)
Pur.wtMM.gradnlm <- nlsLM( ~ weighted.MM.grad(rate, conc1, conc.1, Vm, K),
                         data = lisTreat, start = list(Vm = 200, K = 0.1))
print(Pur.wtMM.gradnlm)
Pur.wtMM.gradnlx <- nlxb( ~ weighted.MM.grad(rate, conc1, conc.1, Vm, K),
                         data = lisTreat, start = list(Vm = 200, K = 0.1),
                         control=list(japprox="jafwd"))
print(Pur.wtMM.gradnlx)
#3 To display the coefficients for comparison:
## rbind(coef(Pur.wtMMnls), coef(Pur.wtMM1nls), coef(Pur.wtMM.gradnls))
## In this example, there seems no advantage to providing the gradient.
## In other cases, there might be.

## ----lls01--------------------------------------------------------------------
# QRsolveEx.R -- example of solving least squares with QR
# J C Nash 2020-12-06
# Enter matrix from Compact Numerical Methods for Computers, Ex04
text <- c(563, 262, 461, 221, 1, 305,
658, 291, 473, 222, 1, 342,
676, 294, 513, 221, 1, 331,
749, 302, 516, 218, 1, 339,
834, 320, 540, 217, 1, 354,
973, 350, 596, 218, 1, 369,
1079, 386, 650, 218, 1, 378,
1151, 401, 676, 225, 1, 368,
1324, 446, 769, 228, 1, 405,
1499, 492, 870, 230, 1, 438,
1690, 510, 907, 237, 1, 438,
1735, 534, 932, 235, 1, 451,
1778, 559, 956, 236, 1, 485)
mm<-matrix(text, ncol=6, byrow=TRUE) # byrow is critical!
A <- mm[,1:5] # select the matrix
A # display
b <- mm[,6] # rhs
b
oss<-as.numeric(crossprod(b))
cat("Overall sumsquares of b =",oss,"\n")
cat("(n-1)*variance=",(length(b)-1)*var(b),"= tss = ",as.numeric(crossprod(b-mean(b))),"\n")
Aqr<-qr(A) # compute the QR decomposition of A
QAqr<-qr.Q(Aqr) # extract the Q matrix of this decomposition
RAqr<-qr.R(Aqr) # This is how to get the R matrix
XAqr<-qr.X(Aqr) # And this is the reconstruction of A from the decomposition
cat("max reconstruction error = ", max(abs(XAqr-A)),"\n")
cat("check the orthogonality Q' * Q: ", max(abs(t(QAqr)%*%QAqr-diag(rep(1,dim(A)[2])))),"\n")
cat("note failure of orthogonality of Q * Q': ",max(abs(QAqr%*%t(QAqr)-diag(rep(1,dim(A)[1])))),"\n")
Absol<-qr.solve(A,b, tol=1000*.Machine$double.eps) # solve the linear ls problem
Absol # and display it 
Abres<-qr.resid(Aqr,b)# and get the residual
rss<-as.numeric(crossprod(Abres))
cat("residual sumsquares=",rss,"\n")
Qtb<-as.vector(t(QAqr)%*%b) # Q' * b -- turns out to be reduction in sumsquares
ssQtb<-as.numeric(crossprod(Qtb))
cat("oss-ssQtb = overall sumsquares minus sumsquares Q' * b = ",oss-ssQtb,"\n")
cat("Thus reduction in sumsquares is ", ssQtb,"\n")
Qtr<-as.vector(t(QAqr)%*%Abres) # projection of residuals onto range space of A
Qtr

## ----maskRvmmin---------------------------------------------------------------
require(optimx)
sq<-function(x){
   nn<-length(x)
   yy<-1:nn
   f<-sum((yy-x)^2)
#   cat("Fv=",f," at ")
#   print(x)
   f
}
sq.g <- function(x){
   nn<-length(x)
   yy<-1:nn
   gg<- 2*(x - yy)
}
xx <- c(.3, 4)
uncans <- Rvmmin(xx, sq, sq.g)
uncans
mybm <- c(0,1) # fix parameter 1
cans <- Rvmmin(xx, sq, sq.g, bdmsk=mybm)
cans

## ----masknlfb, echo=TRUE------------------------------------------------------
## rm(list=ls()) # clear workspace??
library(nlsr)
sqres<-function(x){
   nn<-length(x)
   yy<-1:nn
   res <- (yy-x)
}
sqjac <- function(x){
   nn<-length(x)
   yy<-1:nn
   JJ <- matrix(-1, nrow=nn, ncol=nn)
   attr(JJ, "gradient") <- JJ ## !! Critical statement
   JJ
}
xx <- c(.3, 4) # repeat for completeness
anlfbu<-nlfb(xx, sqres, sqjac)
anlfbu
anlfbc<-nlfb(xx, sqres, sqjac, lower=c(xx[1], -Inf), upper=c(xx[1], Inf), trace=FALSE)
anlfbc
# Following will warn All parameters are masked
anlfbe<-try(nlfb(xx, sqres, sqjac, lower=c(xx[1], xx[2]), upper=c(xx[1],xx[2]))) 

## ----masknlxb, echo=TRUE------------------------------------------------------
nn<-length(xx) # Also length(yy)
if (nn != 2) stop("This example has nn=2 only!")
yy<-1:2
v1 <- c(1,0); v2 <- c(0,1)
anlxbu <- nlxb(yy~v1*p1+v2*p2, start=c(p1=0.3, p2=4) )
anlxbu  
anlxbc <- nlxb(yy~v1*p1+v2*p2, start=c(p1=0.3, p2=4), masked=c("p1") )
anlxbc  

## ----masknlxb2, echo=TRUE-----------------------------------------------------
weed <- c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443,
   38.558, 50.156, 62.948, 75.995, 91.972)
ii <- 1:12
wdf <- data.frame(weed, ii)
weedux <- nlxb(weed~b1/(1+b2*exp(-b3*ii)), start=c(b1=200, b2=50, b3=0.3)) 
weedux
#- Old mechanism of 'nlsr' NO longer works
#- weedcx <- nlxb(weed~b1/(1+b2*exp(-b3*ii)), start=c(b1=200, b2=50, b3=0.3), masked=c("b1")) 
weedcx <- nlxb(weed~b1/(1+b2*exp(-b3*ii)), start=c(b1=200, b2=50, b3=0.3), 
               lower=c(b1=200, b2=0, b3=0), upper=c(b1=200, b2=100, b3=40)) 
weedcx
rfn <- function(bvec, weed=weed, ii=ii){
  res <- rep(NA, length(ii))
  for (i in ii){
    res[i]<- bvec[1]/(1+bvec[2]*exp(-bvec[3]*i))-weed[i]
  }
  res
}
weeduf <- nlfb(start=c(200, 50, 0.3),resfn=rfn,weed=weed, ii=ii, control=list(japprox="jacentral"))
weeduf
#- maskidx method to specify masks no longer works
#- weedcf <- nlfb(start=c(200, 50, 0.3),resfn=rfn,weed=weed, ii=ii, 
#-                maskidx=c(1), control=list(japprox="jacentral"))
weedcf <- nlfb(start=c(200, 50, 0.3),resfn=rfn,weed=weed, ii=ii, lower=c(200, 0, 0), 
               upper=c(200, 100, 40), control=list(japprox="jacentral"))
weedcf

## ----ssmicmen, echo=TRUE------------------------------------------------------
cat("=== SSmicmen ===\n")
dat <- PurTrt <- Puromycin[ Puromycin$state == "treated", ]
frm <- rate ~ SSmicmen(conc, Vm, K)
strt<-getInitial(frm, data = dat)
print(strt)
fm1x <- nlxb(frm, data = dat, start=strt,
             trace=FALSE, control=list(japprox="SSJac"))
fm1x

## ----partlin1-----------------------------------------------------------------
## Comment in example is "## using conditional linearity"
DNase1 <- subset(DNase, Run == 1) # data
## using nls with plinear
# Using a formula that explicitly includes the asymptote. (Default algorithm.)
fm2DNase1orig <- nls(density ~ Asym/(1 + exp((xmid - log(conc))/scal)),
                 data = DNase1,
                 start = list(Asym = 10, xmid = 0, scal = 1))
summary(fm2DNase1orig)
# Using conditional linearity. Note the linear paraemter does not appear explicitly.
#  And we must change the starting vector.
fm2DNase1 <- nls(density ~ 1/(1 + exp((xmid - log(conc))/scal)),
                 data = DNase1,
                 start = list(xmid = 0, scal = 1),
                 algorithm = "plinear")
summary(fm2DNase1)
# How to run with "port" algorithm -- NOT RUN
# fm2DNase1port <- nls(density ~ Asym/(1 + exp((xmid - log(conc))/scal)),
#                  data = DNase1,
#                  start = list(Asym = 10, xmid = 0, scal = 1),
#                  algorithm = "port")
# summary(fm2DNase1port)

# require(minpack.lm) # Does NOT offer VARPRO. First example is WRONG. NOT RUN.
# fm2DNase1mA <- nlsLM(density ~ 1/(1 + exp((xmid - log(conc))/scal)),
#                  data = DNase1,
#                  start = list(xmid = 0, scal = 1))
# summary(fm2DNase1mA)
# fm2DNase1mB <- nlsLM(density ~ Asym/(1 + exp((xmid - log(conc))/scal)),
#                  data = DNase1,
#                  start = list(Asym=10, xmid = 0, scal = 1))
# summary(fm2DNase1mB)
# 
# # This gets wrong answer. NOT RUN
# fm2DNase1xA <- nlxb(density ~ 1/(1 + exp((xmid - log(conc))/scal)),
#                  data = DNase1,
#                  start = list(xmid = 0, scal = 1))
# print(fm2DNase1xA)
# Original formula with nlsr::nlxb().
fm2DNase1xB <- nlxb(density ~ Asym/(1 + exp((xmid - log(conc))/scal)),
                 data = DNase1,
                 start = list(Asym=10, xmid = 0, scal = 1))
print(fm2DNase1xB)

## ----nls-indexed-ex-----------------------------------------------------------
## The muscle dataset in MASS is from an experiment on muscle
## contraction on 21 animals.  The observed variables are Strip
## (identifier of muscle), Conc (Cacl concentration) and Length
## (resulting length of muscle section).
utils::data(muscle, package = "MASS")

## The non linear model considered is
##       Length = alpha + beta*exp(-Conc/theta) + error
## where theta is constant but alpha and beta may vary with Strip.

with(muscle, table(Strip)) # 2, 3 or 4 obs per strip

## We first use the plinear algorithm to fit an overall model,
## ignoring that alpha and beta might vary with Strip.

musc.1 <- nls(Length ~ cbind(1, exp(-Conc/th)), muscle,
              start = list(th = 1), algorithm = "plinear")
summary(musc.1)

## Then we use nls' indexing feature for parameters in non-linear
## models to use the conventional algorithm to fit a model in which
## alpha and beta vary with Strip.  The starting values are provided
## by the previously fitted model.
## Note that with indexed parameters, the starting values must be
## given in a list (with names):
b <- coef(musc.1)
musc.2 <- nls(Length ~ a[Strip] + b[Strip]*exp(-Conc/th), muscle,
              start = list(a = rep(b[2], 21), b = rep(b[3], 21), th = b[1]))
## IGNORE_RDIFF_BEGIN
summary(musc.2)
## IGNORE_RDIFF_END

## ----nls-indexed-start--------------------------------------------------------
istart = list(a = rep(b[2], 21), b = rep(b[3], 21), th = b[1])
str(istart)

## ----appx1, cache=FALSE, echo=TRUE--------------------------------------------
# NOTE: This is OLD material and not consistent in usage with rest of vignette
library(knitr)
# try different ways of supplying data to R nls stuff
ydata <- c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443, 38.558,
           50.156, 62.948, 75.995, 91.972)
ttdata <- seq_along(ydata) # for testing
mydata <- data.frame(y = ydata, tt = ttdata)
hobsc <- y ~ 100*b1/(1 + 10*b2 * exp(-0.1 * b3 * tt))
ste <- c(b1 = 2, b2 = 5, b3 = 3)

# let's try finding the variables
findmainenv <- function(formula, prm) {
  vn <- all.vars(formula)
  pnames <- names(prm)
  ppos <- match(pnames, vn)
  datvar <- vn[-ppos]
  cat("Data variables:")
  print(datvar)
  cat("Are the variables present in the current working environment?\n")
  for (i in seq_along(datvar)){
    cat(datvar[[i]]," : present=",exists(datvar[[i]]),"\n")
  }
}

findmainenv(hobsc, ste)
y <- ydata
tt <- ttdata
findmainenv(hobsc, ste)
rm(y)
rm(tt)

# ===============================


# let's try finding the variables in dotargs
finddotargs <- function(formula, prm, ...) {
  dots <- list(...)
  cat("dots:")
  print(dots)
  cat("names in dots:")
  dtn <- names(dots)
  print(dtn)
  vn <- all.vars(formula)
  pnames <- names(prm)
  ppos <- match(pnames, vn)
  datvar <- vn[-ppos]
  cat("Data variables:")
  print(datvar)
  cat("Are the variables present in the dot args?\n")
  for (i in seq_along(datvar)){
    dname <- datvar[[i]]
    cat(dname," : present=",(dname %in% dtn),"\n")
  }
}

finddotargs(hobsc, ste, y=ydata, tt=ttdata)
# ===============================

y <- ydata
tt <- ttdata
tryq <- try(nlsquiet <- nls(formula=hobsc, start=ste))
if (class(tryq) != "try-error") {print(nlsquiet)} else {cat("try-error\n")}
#- OK
rm(y); rm(tt)

## this will fail
tdots1<-try(nlsdots <- nls(formula=hobsc, start=ste, y=ydata, tt=ttdata))
# But ...
y <- ydata # put data in globalenv
tt <- ttdata
tdots2<-try(nlsdots <- nls(formula=hobsc, start=ste))
tdots2
rm(y); rm(tt)
## but ...
mydata<-data.frame(y=ydata, tt=ttdata)
tframe <- try(nlsframe <- nls(formula=hobsc, start=ste, data=mydata))
if (class(tframe) != "try-error") {print(nlsframe)} else {cat("try-error\n")}
#- OK

y <- ydata
tt <- ttdata
tquiet <- try(nlsrquiet <- nlxb(formula=hobsc, start=ste))
if ( class(tquiet) != "try-error") {print(nlsrquiet)}
#- OK
rm(y); rm(tt)
test <- try(nlsrdots <- nlxb(formula=hobsc, start=ste, y=ydata, tt=ttdata))
#- Note -- does NOT work -- do we need to specify the present env. in nlfb for y, tt??
if (class(test) != "try-error") { print(nlsrdots) } else {cat("Try error\n") }
test2 <- try(nlsframe <- nls(formula=hobsc, start=ste, data=mydata))
if (class(test) != "try-error") {print(nlsframe) } else {cat("Try error\n") }
#- OK


library(minpack.lm)
y <- ydata
tt <- ttdata
nlsLMquiet <- nlsLM(formula=hobsc, start=ste)
print(nlsLMquiet)
#- OK
rm(y)
rm(tt)
## Dotargs
## Following fails
## tdots <- try(nlsLMdots <- nlsLM(formula=hobsc, start=ste, y=ydata, tt=ttdata))
## but ... 
tdots <- try(nlsLMdots <- nlsLM(formula=hobsc, start=ste, data=mydata))
if (class(tdots) != "try-error") { print(nlsLMdots) } else {cat("try-error\n") }
#-  Note -- does NOT work
## dataframe
tframe <- try(nlsLMframe <- nlsLM(formula=hobsc, start=ste, data=mydata) )
if (class(tdots) != "try-error") {print(nlsLMframe)} else {cat("try-error\n") }
#- does not work

## detach("package:nlsr", unload=TRUE)
## Uses nlmrt here for comparison
## library(nlmrt)
## txq <- try( nlxbquiet <- nlxb(formula=hobsc, start=ste))
## if (class(txq) != "try-error") {print(nlxbquiet)} else { cat("try-error\n")}
#- Note -- does NOT work
## txdots <- try( nlxbdots <- nlxb(formula=hobsc, start=ste, y=y, tt=tt) )
## if (class(txdots) != "try-error") {print(nlxbdots)} else {cat("try-error\n")}
#- Note -- does NOT work
## dataframe
## nlxbframe <- nlxb(formula=hobsc, start=ste, data=mydata)
## print(nlxbframe)
#- OK

## ----rndr, code=readLines("../R/numericDerivR.R"), echo=TRUE------------------
numericDerivR <- function(expr, theta, rho = parent.frame(), dir = 1,
                 eps = .Machine$double.eps ^ (1/if(central) 3 else 2), central = FALSE)
## Note: Must this expr must be set up as a call to work properly?
## We set eps conditional on central. But central set AFTER eps. Is this OK?
{   ## cat("numericDeriv-Alt\n")
    dir <- rep_len(dir, length(theta))
    stopifnot(is.finite(eps), eps > 0)
##    rho1 <- new.env(FALSE, rho, 0)
    if (!is.character(theta) ) {stop("'theta' should be of type character")}
    if (is.null(rho)) {
            stop("use of NULL environment is defunct")
            #        rho <- R_BaseEnv;
    } else {
          if(! is.environment(rho)) {stop("'rho' should be an environment")}
          #    int nprot = 3;
    }
    if( ! ((length(dir) == length(theta) ) & (is.numeric(dir) ) ) )
              {stop("'dir' is not a numeric vector of the correct length") }
    if(is.na(central)) { stop("'central' is NA, but must be TRUE or FALSE") }
    res0 <- eval(expr, rho) # the base residuals. C has a check for REAL ANS=res0
    nt <- length(theta) # number of parameters
    mr <- length(res0) # number of residuals
    JJ <- matrix(NA, nrow=mr, ncol=nt) # Initialize the Jacobian
    for (j in 1:nt){
       origPar<-get(theta[j],rho) # This is parameter by NAME, not index
       xx <- abs(origPar)
       delta <- if (xx == 0.0) {eps} else { xx*eps }
       ## JN: I prefer eps*(xx + eps)  which is simpler?
       prmx<-origPar+delta*dir[j]
       assign(theta[j],prmx,rho)
       res1 <- eval(expr, rho) # new residuals (forward step)
#       cat("res1:"); print(res1)
       if (central) { # compute backward step resids for central diff
          prmb <- origPar - dir[j]*delta
          assign(theta[j], prmb, envir=rho) # may be able to make more efficient later?
          resb <- eval(expr, rho)
          JJ[, j] <- dir[j]*(res1-resb)/(2*delta) # vectorized
       } else { ## forward diff
          JJ[, j] <- dir[j]*(res1-res0)/delta
       }  # end forward diff
       assign(theta[j],origPar, rho) # reset value 
    } # end loop over the parameters
    attr(res0, "gradient") <- JJ
    return(res0)
}

## ----ndx, echo=TRUE-----------------------------------------------------------
library(nlsr) # so we have numericDerivR code
# Data for Hobbs problem
weed  <-  c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443, 
            38.558, 50.156, 62.948, 75.995, 91.972) # for testing
tt  <-  seq_along(weed) # for testing
# A simple starting vector -- must have named parameters for nlxb, nls, wrapnlsr.
st  <-  c(b1=1, b2=1, b3=1)
wmodu  <-   weed ~ b1/(1+b2*exp(-b3*tt))
weeddf  <-  data.frame(weed=weed, tt=tt)
weedenv <- list2env(weeddf)
weedenv$b1 <- st[[1]]
weedenv$b2 <- st[[2]]
weedenv$b3 <- st[[3]]
rexpr<-call("-",wmodu[[3]], wmodu[[2]])
r0<-eval(rexpr, weedenv)
cat("Sumsquares at 1,1,1 is ",sum(r0^2),"\n")## Another way
expr <- wmodu
rho <- weedenv
rexpr<-call("-",wmodu[[3]], wmodu[[2]])
res0<-eval(rexpr, rho) # the base residuals
res0
## Try the numericDeriv option
theta<-names(st)
nDnls<-numericDeriv(rexpr, theta, weedenv)
nDnls
nDnlsR<-numericDerivR(rexpr, theta, weedenv)
nDnlsR

## ----derivexp, eval=TRUE------------------------------------------------------
partialderiv <- D(expression(a * (x ^ b)),"b")
print(partialderiv)

## ----timingshobbsmodel--------------------------------------------------------
require(microbenchmark)
## nls on Hobbs scaled model
wmods  <-   weed ~ 100*b1/(1+10*b2*exp(-0.1*b3*tt))
stx<-c(b1=2, b2=5, b3=3)
tnls<-microbenchmark((anls<-nls(wmods, start=stx, data=weeddf)), unit="us")
tnls
## nlsr::nlfb() on Hobbs scaled model
tnlxb<-microbenchmark((anlxb<-nlsr::nlxb(wmods, start=stx, data=weeddf)), unit="us")
tnlxb
## minpack.lm::nlsLM() on Hobbs scaled model
tnlsLM<-microbenchmark((anlsLM<-minpack.lm::nlsLM(start=stx, formula=wmods, data=weeddf)), 
                       unit="us")
tnlsLM

## ----timingnlsrminpackfn------------------------------------------------------
require(microbenchmark)
## nlsr::nlfb() on Hobbs scaled
tnlfb<-microbenchmark((anlfb<-nlsr::nlfb(start=st1, resfn=shobbs.res, jacfn=shobbs.jac)), unit="us")
tnlfb
## minpack.lm::nls.lm() on Hobbs scaled
tnls.lm<-microbenchmark((anls.lm<-minpack.lm::nls.lm(par=st1, fn=shobbs.res, jac=shobbs.jac)))
tnls.lm

## ----hiddenexprob1, echo=FALSE------------------------------------------------
# Here we set up an example problem with data
# Define independent variable
t0 <- 0:19
t0a<-t0
t0a[1]<-1e-20 # very small value
# Drop first value in vectors
t0t<-t0[-1]
y1 <- 4 * (t0^0.25)
y1t<-y1[-1]
n <- length(t0)
fuzz <- rnorm(n)
range <- max(y1)-min(y1)
## add some "error" to the dependent variable
y1q <- y1 + 0.2*range*fuzz
edta <- data.frame(t0=t0, t0a=t0a, y1=y1, y1q=y1q)
edtat <- data.frame(t0t=t0t, y1t=y1t)

start1 <- c(a=1, b=1)
try(nlsy0t0ax <- nls(formula=y1~a*(t0a^b), start=start1, data=edta, control=nls.control(maxiter=10000)))
nlsry1t0a <- nlxb(formula=y1~a*(t0a^b), start=start1, data=edta)
library(minpack.lm)
nlsLMy1t0a <- nlsLM(formula=y1~a*(t0a^b), start=start1, data=edta)

## ----nlsstr-------------------------------------------------------------------
str(nlsy0t0ax)

## ----nlsrstr------------------------------------------------------------------
str(nlsry1t0a)

## ----pred1--------------------------------------------------------------------
nudta <- colMeans(edta)
predict(nlsy0t0ax, newdata=nudta)
predict(nlsLMy1t0a, newdata=nudta)
predict(nlsry1t0a, newdata=nudta)

## ----exprob1------------------------------------------------------------------
# Here we set up an example problem with data
# Define independent variable
t0 <- 0:19
t0a<-t0
t0a[1]<-1e-20 # very small value
# Drop first value in vectors
t0t<-t0[-1]
y1 <- 4 * (t0^0.25)
y1t<-y1[-1]
n <- length(t0)
fuzz <- rnorm(n)
range <- max(y1)-min(y1)
## add some "error" to the dependent variable
y1q <- y1 + 0.2*range*fuzz
edta <- data.frame(t0=t0, t0a=t0a, y1=y1, y1q=y1q)
edtat <- data.frame(t0t=t0t, y1t=y1t)

## ----extrynls-----------------------------------------------------------------
cprint <- function(obj){
   # print object if it exists
  sobj<-deparse(substitute(obj))
  if (exists(sobj)) {
      print(obj)
  } else {
      cat(sobj," does not exist\n")
  }
#  return(NULL)
}
start1 <- c(a=1, b=1)
try(nlsy0t0 <- nls(formula=y1~a*(t0^b), start=start1, data=edta))
cprint(nlsy0t0)
# Since this fails to converge, let us increase the maximum iterations
try(nlsy0t0x <- nls(formula=y1~a*(t0^b), start=start1, data=edta,
                    control=nls.control(maxiter=10000)))
cprint(nlsy0t0x)
try(nlsy0t0ax <- nls(formula=y1~a*(t0a^b), start=start1, data=edta, 
                     control=nls.control(maxiter=10000)))
cprint(nlsy0t0ax)
try(nlsy0t0t <- nls(formula=y1t~a*(t0t^b), start=start1, data=edtat))
cprint(nlsy0t0t)

## ----extry1nlsr---------------------------------------------------------------
nlsry1t0 <- try(nlxb(formula=y1~a*(t0^b), start=start1, data=edta))
cprint(nlsry1t0)
nlsry1t0a <- nlxb(formula=y1~a*(t0a^b), start=start1, data=edta)
cprint(nlsry1t0a)
nlsry1t0t <- nlxb(formula=y1t~a*(t0t^b), start=start1, data=edtat)
cprint(nlsry1t0t)

## ----extry1minlm--------------------------------------------------------------
library(minpack.lm)
nlsLMy1t0 <- nlsLM(formula=y1~a*(t0^b), start=start1, data=edta)
nlsLMy1t0
nlsLMy1t0a <- nlsLM(formula=y1~a*(t0a^b), start=start1, data=edta)
nlsLMy1t0a
nlsLMy1t0t <- nlsLM(formula=y1t~a*(t0t^b), start=start1, data=edtat)
nlsLMy1t0t

## ----brownden-----------------------------------------------------------------
m <- 20
t <- seq(1, m) / 5
y <- rep(0,m)
library(nlsr)
library(minpack.lm)

bddata <- data.frame(t=t, y=y)
bdform <- y ~ ((x1 + t * x2 - exp(t))^2 + (x3 + x4 * sin(t) - cos(t))^2)
prm0 <- c(x1=25, x2=5, x3=-5, x4=-1)
fbd <-model2ssgrfun(bdform, prm0, bddata)
cat("initial sumsquares=",as.numeric(crossprod(fbd(prm0))),"\n")
nlsrbd <- nlxb(bdform, start=prm0, data=bddata, trace=FALSE)
nlsrbd

nlsbd10k <- nls(bdform, start=prm0, data=bddata, trace=FALSE, 
                control=nls.control(maxiter=10000))
nlsbd10k

nlsLMbd10k <- nlsLM(bdform, start=prm0, data=bddata, trace=FALSE, 
                    control=nls.lm.control(maxiter=10000, maxfev=10000))
nlsLMbd10k

## ----browndenpred-------------------------------------------------------------
ndata <- data.frame(t=c(5,6), y=c(0,0))
predict(nlsLMbd10k, newdata=ndata)

# now nls
predict(nlsbd10k, newdata=ndata)

# now nlsr
predict(nlsrbd, newdata=ndata)

## ----brownden2----------------------------------------------------------------
bdf2 <-  (x1 + t * x2 - exp(t))^2 ~ - (x3 + x4 * sin(t) - cos(t))^2

## ----bdcheck, eval=FALSE------------------------------------------------------
# #' Brown and Dennis Function
# #'
# #' Test function 16 from the More', Garbow and Hillstrom paper.
# #'
# #' The objective function is the sum of \code{m} functions, each of \code{n}
# #' parameters.
# #'
# #' \itemize{
# #'   \item Dimensions: Number of parameters \code{n = 4}, number of summand
# #'   functions \code{m >= n}.
# #'   \item Minima: \code{f = 85822.2} if \code{m = 20}.
# #' }
# #'
# #' @param m Number of summand functions in the objective function. Should be
# #'   equal to or greater than 4.
# #' @return A list containing:
# #' \itemize{
# #'   \item \code{fn} Objective function which calculates the value given input
# #'   parameter vector.
# #'   \item \code{gr} Gradient function which calculates the gradient vector
# #'   given input parameter vector.
# #'   \item \code{fg} A function which, given the parameter vector, calculates
# #'   both the objective value and gradient, returning a list with members
# #'   \code{fn} and \code{gr}, respectively.
# #'   \item \code{x0} Standard starting point.
# #' }
# #' @references
# #' More', J. J., Garbow, B. S., & Hillstrom, K. E. (1981).
# #' Testing unconstrained optimization software.
# #' \emph{ACM Transactions on Mathematical Software (TOMS)}, \emph{7}(1), 17-41.
# #' \url{https://doi.org/10.1145/355934.355936}
# #'
# #' Brown, K. M., & Dennis, J. E. (1971).
# #' \emph{New computational algorithms for minimizing a sum of squares of
# #' nonlinear functions} (Report No. 71-6).
# #' New Haven, CT: Department of Computer Science, Yale University.
# #'
# #' @examples
# #' # Use 10 summand functions
# #' fun <- brown_den(m = 10)
# #' # Optimize using the standard starting point
# #' x0 <- fun$x0
# #' res_x0 <- stats::optim(par = x0, fn = fun$fn, gr = fun$gr, method =
# #' "L-BFGS-B")
# #' # Use your own starting point
# #' res <- stats::optim(c(0.1, 0.2, 0.3, 0.4), fun$fn, fun$gr, method =
# #' "L-BFGS-B")
# #'
# #' # Use 20 summand functions
# #' fun20 <- brown_den(m = 20)
# #' res <- stats::optim(fun20$x0, fun20$fn, fun20$gr, method = "L-BFGS-B")
# #' @export
# #`
# brown_den <- function(m = 20) {
#   list(
#     fn = function(par) {
#       x1 <- par[1]
#       x2 <- par[2]
#       x3 <- par[3]
#       x4 <- par[4]
# 
#       ti <- (1:m) * 0.2
#       l <- x1 + ti * x2 - exp(ti)
#       r <- x3 + x4 * sin(ti) - cos(ti)
#       f <- l * l + r * r
#       sum(f * f)
#     },
#     gr = function(par) {
#       x1 <- par[1]
#       x2 <- par[2]
#       x3 <- par[3]
#       x4 <- par[4]
# 
#       ti <- (1:m) * 0.2
#       sinti <- sin(ti)
#       l <- x1 + ti * x2 - exp(ti)
#       r <- x3 + x4 * sinti - cos(ti)
#       f <- l * l + r * r
#       lf4 <- 4 * l * f
#       rf4 <- 4 * r * f
#       c(
#         sum(lf4),
#         sum(lf4 * ti),
#         sum(rf4),
#         sum(rf4 * sinti)
#       )
#     },
#     fg = function(par) {
#       x1 <- par[1]
#       x2 <- par[2]
#       x3 <- par[3]
#       x4 <- par[4]
# 
#       ti <- (1:m) * 0.2
#       sinti <- sin(ti)
#       l <- x1 + ti * x2 - exp(ti)
#       r <- x3 + x4 * sinti - cos(ti)
#       f <- l * l + r * r
#       lf4 <- 4 * l * f
#       rf4 <- 4 * r * f
# 
#       fsum <- sum(f * f)
#       grad <- c(
#         sum(lf4),
#         sum(lf4 * ti),
#         sum(rf4),
#         sum(rf4 * sinti)
#       )
# 
#       list(
#         fn = fsum,
#         gr = grad
#       )
#     },
#     x0 = c(25, 5, -5, 1)
#   )
# }
# mbd <- brown_den(m=20)
# mbd
# mbd$fg(mbd$x0)
# bdsolnm <- optim(mbd$x0, mbd$fn, control=list(trace=0))
# bdsolnm
# bdsolbfgs <- optim(mbd$x0, mbd$fn, method="BFGS", control=list(trace=0))
# bdsolbfgs
# 
# library(optimx)
# methlist <- c("Nelder-Mead","BFGS","Rvmmin","L-BFGS-B","Rcgmin","ucminf")
# 
# solo <- opm(mbd$x0, mbd$fn, mbd$gr, method=methlist, control=list(trace=0))
# summary(solo, order=value)
# 
# ## A failure above is generally because a package in the 'methlist' is not installed.

