## ----ex01, echo=TRUE----------------------------------------------------------
weed <- c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443,
          38.558, 50.156, 62.948, 75.995, 91.972)
tt <- 1:12
weeddf <- data.frame(tt, weed)
plot(weeddf, main="Hobbs weed infestation data")

## ----ex02set, echo=TRUE-------------------------------------------------------
# model formulas
frmu <- weed ~ b1/(1+b2*exp(-b3*tt))
frms <- weed ~ 100*c1/(1+10*c2*exp(-0.1*c3*tt))
frmt <- weed ~ Asym /(1 + exp((xmid-tt)/scal))
#
# Starting parameter sets
stu1<-c(b1=1, b2=1, b3=1)
sts1<-c(c1=1, c2=1, c3=1)
stt1<-c(Asym=1, xmid=1, scal=1)

## ----ex02nls------------------------------------------------------------------
unls1<-try(nls(formula=frmu, start=stu1, data=weeddf))
summary(unls1)
snls1<-try(nls(formula=frms, start=sts1, data=weeddf))
summary(snls1)
tnls1<-try(nls(formula=frmt, start=stt1, data=weeddf))
summary(tnls1)
cat("\n")

## ----ex02nlsr-----------------------------------------------------------------
library(nlsr)
unlx1<-try(nlxb(formula=frmu, start=stu1, data=weeddf))
print(unlx1) 
pshort(unlx1) # A short form output
snlx1<-try(nlxb(formula=frms, start=sts1, data=weeddf))
# pshort(snlx1)
print(snlx1)
tnlx1<-try(nlxb(formula=frmt, start=stt1, data=weeddf))
# pshort(tnlx1)
print(tnlx1)
ct<-as.list(tnlx1$coefficients) # Need a list to use ct$... in next line
cat("exp(xmid/scal)=",exp(ct$xmid/ct$scal),"\n")
cat("\n")
# explicit residuals (weighted if there are weights)
rtnlx1<-residuals(tnlx1)
print(rtnlx1)
cat("explicit sumsquares =", sum(rtnlx1^2),"\n")

## ----ex02minpack--------------------------------------------------------------
library(minpack.lm)
unlm1<-try(nlsLM(formula=frmu, start=stu1, data=weeddf))
summary(unlm1) # Use summary() to get display
unlm1 # Note the difference. Use this form to get sum of squares
# pnls(unlm1)  # Short form of output
snlm1<-try(nlsLM(formula=frms, start=sts1, data=weeddf))
# pnls(snlm1)
summary(snlm1)
tnlm1<-try(nlsLM(formula=frmt, start=stt1, data=weeddf))
if (inherits(tnlm1, "try-error")) {
   cat("Failure to compute solution -- likely singular Jacobian\n")
} else {  
   pnls(tnlm1) # short form to give sum of squares
   summary(tnlm1)
}   

## ----ex02wrapnlsr-------------------------------------------------------------
## Try the wrapper. Calling wrapnlsr() instead of nlsr() is equivalent
##
unlw1<-try(nlsr(formula=frmu, start=stu1, data=weeddf))
print(unlw1) # using 'unlx1' gives name of object as 'x' only
snlw1<-try(nlsr(formula=frms, start=sts1, data=weeddf))
# pshort(snlx1)
print(snlw1)
tnlw1<-try(nlsr(formula=frmt, start=stt1, data=weeddf))
print(tnlw1)

## ----ex03, echo=TRUE----------------------------------------------------------
stspecial<- c(Asym = 35.532,  xmid = 43376,  scal = -2935.4)
# force to exit at once by setting femax to 1 (maximum number of sum of squares evaluations)
getsvs<-nlxb(formula=frmt, start=stspecial, data=weeddf, control=list(femax=1))
print(getsvs)

## ----ex04, echo=TRUE, eval=FALSE----------------------------------------------
# if (inherits(tnlm1, "try-error")) {
#    cat("Cannot compute solution -- likely singular Jacobian\n")
# } else {
#    Jtm <- tnlm1$m$gradient()
#    svd(Jtm)$d # Singular values
# }

## ----ex05, echo=TRUE----------------------------------------------------------
stspecial<- c(Asym = 35.532,  xmid = 43376,  scal = -2935.4)
# force to exit at once by setting femax to 1 (maximum number of sum of squares evaluations)
badstart<-nlxb(formula=frmt, start=stspecial, data=weeddf)
print(badstart)

## ----ex06, echo=TRUE----------------------------------------------------------
frmtss <- weed ~ SSlogis(tt, Asym, xmid, scal)
ssts1<-nls(formula=frmtss, data=weeddf)
summary(ssts1)
require(minpack.lm)
sstm1<-nlsLM(formula=frmtss, data=weeddf)
summary(sstm1)

## ----ex07, echo=TRUE----------------------------------------------------------
Asymt <- 2*max(weed)
Asymt
pw <- Asymt/weed
pw1<-pw-1
lpw1<-log(pw1)
clin <- coef(lm(lpw1~tt))
clin<-as.list(clin)
scalt <- -1/clin$tt
scalt
xmidt <- scalt*as.numeric(clin[1])
xmidt
library(nlsr)
try1<-nlxb(frmt, data=weeddf, start=c(Asym=1, xmid=1, scal=1))
pshort(try1)
try2<-nlxb(frmt, data=weeddf, start=c(Asym=Asymt, xmid=xmidt, scal=scalt))
pshort(try2)

## ----ex08, eval=TRUE----------------------------------------------------------
frmtssJ <- weed ~ SSlogisJN(tt, Asym, xmid, scal) # The model formula
lstrt <- getInitial(frmtssJ, weeddf) # Here we get the starting parameters
cat("starting parameters:\n")
print(lstrt)
# Next line uses the start. We indicate that the "approximate" Jacobian code is in
#  the selfStart code, though in fact it is not an approximation in this case.
sstx1<-nlxb(formula=frmtssJ, start=lstrt, data=weeddf, control=list(japprox="SSJac"))
print(sstx1)
# If no Jacobian code is included in selfStart module, we can actually use an
#  approximate Jacobian. See "gradient" elements for differences in result.
sstx1a<-nlxb(formula=frmtssJ, start=lstrt, data=weeddf, control=list(japprox="jacentral"))
print(sstx1a)
## We can automate the above in function nlsrSS()
sstSS<-nlsrSS(formula=frmtssJ, data=weeddf)
print(sstSS)

## ----ex09, echo=TRUE----------------------------------------------------------
b1t <- 2*max(weed)
b1t
lpw1<-log(b1t/weed-1)
clin <- as.list(coef(lm(lpw1~tt)))
b3t <- -clin$tt
b3t
b2t <- exp(as.numeric(clin[1]))
b2t
library(nlsr)
try1<-nlxb(frmu, data=weeddf, start=c(b1=1, b2=1, b3=1))
pshort(try1)
try2<-nlxb(frmu, data=weeddf, start=c(b1=b1t, b2=b2t, b3=b3t))
pshort(try2)

## ----ex10, echo=TRUE----------------------------------------------------------
# Start MUST be feasible i.e. on or within bounds
anlshob1b <- nls(frms, start=sts1, data=weeddf, lower=c(0,0,0),
             upper=c(2,6,3), algorithm='port')
pnls(anlshob1b) #  check the answer (short form)
# nlsLM seems NOT to work with bounds in this example
anlsLM1b <- nlsLM(frms, start=sts1, data=weeddf, lower=c(0,0,0), upper=c(2,6,3))
pnls(anlsLM1b)
# also no warning if starting out of bounds, but gets good answer!!
st4<-c(c1=4, c2=4, c3=4)
anlsLMob <- nlsLM(frms, start=st4, data=weeddf, lower=c(0,0,0), upper=c(2,6,3))
pnls(anlsLMob)
# Try nlsr::nlxb()
anlx1b <- nlxb(frms, start=sts1, data=weeddf, lower=c(0,0,0), upper=c(2,6,3))
pshort(anlx1b)

## ----exrosen, eval=TRUE-------------------------------------------------------
frres<-function(x) { ## Rosenbrock as residuals
    x1 <- x[1]
    x2 <- x[2]
    res<-c( 10 * (x2 - x1 * x1), (1 - x1) )
}

frjac<-function(x) { ## Rosenbrock residual derivatives matrix
    x1 <- x[1]
    x2 <- x[2]
   J<-matrix(NA, nrow=2, ncol=2)
   J[1,1]<- -20*x1
   J[1,2]<- 10
   J[2,1]<- -1
   J[2,2]<- 0
   attr(J,"gradient")<-J # NEEDED for nlfb()
   J
}
require(minpack.lm)
require(nlsr)
strt<-c(-1.2,1)
rosnlf<-nlfb(strt, resfn=frres, jacfn=frjac)
print(rosnlf)
rosnlm<-nls.lm(par=strt, fn=frres, jac=frjac)
summary(rosnlm)  

## ----staticwts----------------------------------------------------------------
wts <- 1/weeddf$tt # wts are reciprocal of time value
tnlx1w<-try(nlxb(formula=frmt, start=stt1, data=weeddf, weights=wts))
# pshort(tnlx1)
print(tnlx1w)
ct<-as.list(tnlx1w$coefficients) # Need a list to use ct$... in next line
cat("exp(xmid/scal)=",exp(ct$xmid/ct$scal),"\n")
cat("\n")

rtnlx1w<-residuals(tnlx1w)
print(rtnlx1w)
cat("explicit sumsquares =", sum(rtnlx1w^2),"\n")

## ----puromycinwfct, echo=TRUE, eval=FALSE-------------------------------------
# library(minpack.lm)
# library(nlsr) # for pnls
# Treated <- Puromycin[Puromycin$state == "treated", ]
# # First get raw estimates
# wtt3nlm0<-nlsLM(rate ~ Vm * conc/(K + conc), data = Treated, trace=TRUE,
#                start = c(Vm = 200, K = 0.05))
# fit0<-fitted(wtt3nlm0)
# # Static wts using fit0
# wtt3nlm0s<-nlsLM(rate ~ Vm * conc/(K + conc), data = Treated, trace=TRUE,
#                start = c(Vm = 200, K = 0.05), weights = wfct(1/fit0^2))
# pnls(wtt3nlm0s)
# # and run directly, noting the 2 phase operation
# wtt3nlm<-nlsLM(rate ~ Vm * conc/(K + conc), data = Treated, trace=TRUE,
#                start = c(Vm = 200, K = 0.05), weights = wfct(1/fitted^2))
# pnls(wtt3nlm)
# cat("weights from wtt3nlm\n")
# as.numeric(wtt3nlm$weights)

## ----formulawts---------------------------------------------------------------
library(nlsr)
Treated <- Puromycin[Puromycin$state == "treated", ]
# "regression" formula
w1frm <- rate ~ Vm * conc/(K + conc)
# Explicit dynamic weighted nlxb
w3frm <- rate/(Vm * conc/(K + conc)) ~ 1
dynweighted <- nlxb(w3frm, data = Treated, start = c(Vm = 201.003, K = 0.04696),
                    trace=TRUE, control=list(prtlvl=0))
pshort(dynweighted)
# Compute the model from the ORIGINAL model (w1frm) using parameters from dynweighted
dyn0 <- nlxb(w1frm, start=coefficients(dynweighted), data=Treated, control=list(jemax=0))
wfdyn0 <- 1/fitted(dyn0)^2 # weights
print(as.numeric(wfdyn0))
pshort(dyn0) # Shows sumsquares without weights, but computes weights we used
formulaweighted <- nlxb(w1frm, data = Treated, start = c(Vm = 201.003, K = 0.04696),
                        weights = ~ 1/fitted^2, trace=TRUE, control=list(prtlvl=0))
pshort(formulaweighted)
wtsfromfits<-1/fitted(formulaweighted)^2
print(as.numeric(wtsfromfits))
print(as.numeric(formulaweighted$weights))

