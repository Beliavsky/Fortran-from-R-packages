## ----setup, include = FALSE---------------------------------------------------
library(calibrar)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----echo=TRUE, eval=FALSE, results='markup'----------------------------------
#  optim2(
#    par,
#    fn,
#    gr = NULL,
#    ...,
#    method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent", "nlm", "nlminb",
#      "Rcgmin", "Rvmmin", "hjn", "spg", "LBFGSB3", "AHR-ES"),
#    lower = -Inf,
#    upper = +Inf,
#    active = NULL,
#    control = list(),
#    hessian = FALSE,
#    parallel = FALSE
#  )

## ----echo=TRUE, results='markup'----------------------------------------------
library(calibrar)
optim(par=rep(1, 5), fn=function(x) sum(x^2))

## ----echo=TRUE, results='markup'----------------------------------------------
optim2(par=rep(1, 5), fn=function(x) sum(x^2))

## ----echo=TRUE, results='markup'----------------------------------------------
optim2(par=rep(1, 5), fn=function(x) sum(x^2), method="nlm")

## ----echo=TRUE, results='markup'----------------------------------------------
set.seed(880820) # for reproducibility
optim2(par=rep(1, 5), fn=function(x) sum(x^2), method="AHR-ES")

## ----echo=TRUE, results='markup'----------------------------------------------
optim2(par=rep(1, 5), fn=function(x) sum(x^2), 
       active=c(TRUE, TRUE, FALSE, FALSE, TRUE))

## ----echo=TRUE, results='markup'----------------------------------------------
optim2(par=rep(1, 5), fn=function(x) sum(x^2), parallel=TRUE)

## ----echo=TRUE, results='hide', eval=FALSE------------------------------------
#  optim2(par=rep(0.5, 5), fn=function(x) sum(2*x^(3.1*x)), control=list(gr.method="richardson"))
#  optim2(par=rep(0.5, 5), fn=function(x) sum(2*x^(3.1*x)), control=list(gr.method="central"))
#  optim2(par=rep(0.5, 5), fn=function(x) sum(2*x^(3.1*x)), control=list(gr.method="forward"))
#  

## ----echo=TRUE, eval=FALSE, results='markup'----------------------------------
#  optimh(
#    par,
#    fn,
#    gr = NULL,
#    ...,
#    method = c("AHR-ES", "Nelder-Mead", "SANN", "hjn", "CMA-ES", "genSA", "DE", "soma",
#      "genoud", "PSO", "hybridPSO", "mads", "hjk", "hjkb", "nmk", "nmkb"),
#    lower = -Inf,
#    upper = +Inf,
#    active = NULL,
#    control = list(),
#    hessian = FALSE,
#    parallel = FALSE
#  )

## ----echo=TRUE, results='markup'----------------------------------------------
# Covariance Matrix Adaptation Evolutionary Strategy
set.seed(880820) # for reproducibility
optimh(par=rep(1, 5), fn=function(x) sum(x^2), method="CMA-ES",
       control=list(maxit=200))

## ----echo=TRUE, results='markup'----------------------------------------------
# Generalized Simulated Anneling
set.seed(880820) # for reproducibility
optimh(par=rep(1, 5), fn=function(x) sum(x^2), method="genSA", 
       lower=rep(-100, 5), upper=rep(100, 5),
       control=list(maxit=200, temperature=6000))

## ----echo=TRUE, results='markup'----------------------------------------------
# Self-Organising Migrating Algorithm
set.seed(880820) # for reproducibility
optimh(par=rep(1, 5), fn=function(x) sum(x^2), method="soma",
       lower=rep(-100, 5), upper=rep(100, 5),
       control=list(maxit=200))

## ----echo=FALSE, results='hide', eval=TRUE, message=FALSE---------------------
library(parallel)

## ----echo=TRUE, results='markup', eval=FALSE----------------------------------
#  library(parallel)
#  ncores = detectCores() - 1 # number of cores to be used
#  cl = makeCluster(ncores)
#  # this is slower than sequential for very fast models (like this one)
#  optim2(par=rep(0.5, 5), fn=function(x) sum(x^2),
#                 control=list(ncores=ncores), parallel=TRUE)
#  stopCluster(cl) # close the parallel connections

## ----echo=TRUE, eval=FALSE, results='markup'----------------------------------
#  calibrate(
#    par,
#    fn,
#    gr = NULL,
#    ...,
#    method = NULL,
#    lower = NULL,
#    upper = NULL,
#    phases = NULL,
#    control = list(),
#    hessian = FALSE,
#    replicates = 1,
#    parallel = FALSE
#  )

## ----echo=TRUE, results='markup'----------------------------------------------
calibrate(par=c(1,2,3,NA,NA), fn=function(x) sum(x^2))

## ----echo=TRUE, results='markup'----------------------------------------------
calibrate(par=c(1,2,3,NA,5), fn=function(x) sum(x^2),
          lower=rep(-100, 5), upper=rep(100, 5))

## ----echo=TRUE, results='markup'----------------------------------------------
calibrate(par=c(1,2,3,NA,5), fn=function(x) sum(x^2),
          lower=rep(-100, 5), upper=rep(100, 5),
          phases=c(1,2,3,2,1))

## ----echo=TRUE, results='markup'----------------------------------------------
calibrate(par=c(1,2,3,NA,5), fn=function(x) sum(x^2),
          lower=rep(-100, 5), upper=rep(100, 5),
          phases=c(1,2,-1,2,1))

## ----echo=TRUE, results='markup'----------------------------------------------
calibrate(par=c(1,2,3,NA,5), fn=sphereN,
          lower=rep(-100, 5), upper=rep(100, 5),
          phases=c(1,2,3,2,1), replicates=3, control=list(maxit=1000))

## ----echo=TRUE, results='markup'----------------------------------------------
calibrate(par=c(1,2,3,NA,5), fn=sphereN,
          lower=rep(-100, 5), upper=rep(100, 5),
          phases=c(1,2,3,2,1), replicates=c(1,1,5), control=list(maxit=1000))

## ----echo=TRUE, results='markup'----------------------------------------------
calibrate(par=list(par1=c(1,2,3), par2=NA, par3=5), fn=sphereN,
          lower=rep(-100, 5), upper=rep(100, 5),
          phases=c(1,2,-3,2,1), replicates=c(1,5), control=list(maxit=1000))

## ----echo=TRUE, results='hide', eval=FALSE------------------------------------
#  library(parallel)
#  ncores = detectCores() - 1 # number of cores to be used
#  cl = makeCluster(ncores)
#  # this is slower than sequential for very fast models (like this one)
#  calib = calibrate(par=rep(0.5, 5), fn=sphereN,
#                    replicates=3,
#                    lower=rep(-5, 5),
#                    upper=rep(+5, 5),
#                    phases=c(1,1,1,2,3),
#                    control=list(parallel=TRUE, ncores=ncores))
#  stopCluster(cl) # close the parallel connections

