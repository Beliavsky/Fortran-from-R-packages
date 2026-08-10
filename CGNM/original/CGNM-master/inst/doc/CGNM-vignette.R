## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(CGNM)
library(knitr)

## ----define the model function------------------------------------------------
model_function=function(x){

  observation_time=c(0.1,0.2,0.4,0.6,1,2,3,6,12)
  Dose=1000
  F=1

  ka=x[1]
  V1=x[2]
  CL_2=x[3]
  t=observation_time

  Cp=ka*F*Dose/(V1*(ka-CL_2/V1))*(exp(-CL_2/V1*t)-exp(-ka*t))

  log10(Cp)
}


## -----------------------------------------------------------------------------
observation=log10(c(4.91, 8.65, 12.4, 18.7, 24.3, 24.5, 18.4, 4.66, 0.238))

## ----warning = FALSE----------------------------------------------------------
CGNM_result=Cluster_Gauss_Newton_method(nonlinearFunction=model_function,targetVector = observation,
initial_lowerRange =rep(0.01,3),initial_upperRange =  rep(100,3),lowerBound = rep(0,3), saveLog=TRUE, num_minimizersToFind = 500, ParameterNames = c("Ka","V1","CL"))

## -----------------------------------------------------------------------------
kable(head(acceptedApproximateMinimizers(CGNM_result)))

## -----------------------------------------------------------------------------
kable(table_parameterSummary(CGNM_result))

## -----------------------------------------------------------------------------
CGNM_bootstrap=Cluster_Gauss_Newton_Bootstrap_method(CGNM_result, nonlinearFunction=model_function)

## -----------------------------------------------------------------------------
kable(table_parameterSummary(CGNM_bootstrap))

## -----------------------------------------------------------------------------
library(ggplot2)

## ----fig.width=6, fig.height=3.5----------------------------------------------
plot_Rank_SSR(CGNM_result)

## ----fig.width=6, fig.height=3.5----------------------------------------------
plot_paraDistribution_byHistogram(CGNM_bootstrap, bins = 50)+scale_x_continuous(trans="log10")

## ----fig.width = 7------------------------------------------------------------
plot_goodnessOfFit(CGNM_result, plotType = 1, independentVariableVector = c(0.1,0.2,0.4,0.6,1,2,3,6,12), plotRank = seq(1,50))

## ----fig.width = 7------------------------------------------------------------
plot_goodnessOfFit(CGNM_bootstrap, plotType = 1, independentVariableVector = c(0.1,0.2,0.4,0.6,1,2,3,6,12))

## ----fig.width = 7------------------------------------------------------------
plot_profileLikelihood(c("CGNM_log","CGNM_log_bootstrap"))+scale_x_continuous(trans="log10")

## -----------------------------------------------------------------------------
kable(table_profileLikelihoodConfidenceInterval(c("CGNM_log","CGNM_log_bootstrap"), alpha = 0.25))

## ----fig.width = 7, fig.height = 6--------------------------------------------
plot_2DprofileLikelihood(CGNM_result, showInitialRange=FALSE, alpha = 0.05)+scale_x_continuous(trans="log10")+scale_y_continuous(trans="log10")

## -----------------------------------------------------------------------------
model_matrix_function=function(x){
  Y_list=lapply(split(x, rep(seq(1:nrow(x)),ncol(x))), model_function)
  Y=t(matrix(unlist(Y_list),ncol=length(Y_list)))
}

testX=t(matrix(c(rep(0.01,3),rep(10,3),rep(100,3)), nrow = 3))
print("testX")
print(testX)
print("model_matrix_function(testX)")
print(model_matrix_function(testX))

print("model_matrix_function(testX)-rbind(model_function(testX[1,]),model_function(testX[2,]),model_function(testX[3,]))")
print(model_matrix_function(testX)-rbind(model_function(testX[1,]),model_function(testX[2,]),model_function(testX[3,])))


## -----------------------------------------------------------------------------

# library(parallel)
#
# obsLength=length(observation)
# 
## Given CGNM searches through wide range of parameter combination, it can encounter
## parameter combinations that is not feasible to evaluate. This try catch function
## is implemented within CGNM for regular functions but for the matrix functions
## user needs to implement outside of CGNM
#
# modelFunction_tryCatch=function(x_in){
#  out=tryCatch({model_function(x_in)},
#               error=function(cond) {rep(NA, obsLength)}
#  )
#  return(out)
# }
# 
#  model_matrix_function=function(X){
#   Y_list=mclapply(split(X, rep(seq(1:nrow(X)),ncol(X))), modelFunction_tryCatch,mc.cores = (parallel::detectCores()-1), mc.preschedule = FALSE)
#
#   Y=t(matrix(unlist(Y_list),ncol=length(Y_list)))
# 
#   return(Y)
#  }


## -----------------------------------------------------------------------------
#library(foreach)
#library(doParallel)

#numCore=8
#registerDoParallel(numCore-1)
#cluster=makeCluster(numCore-1, type = "PSOCK")
#registerDoParallel(cl=cluster)

# obsLength=length(observation)

## Given CGNM searches through wide range of parameter combination, it can encounter
## parameter combinations that is not feasible to evaluate. This try catch function
## is implemented within CGNM for regular functions but for the matrix functions
## user needs to implement outside of CGNM

# modelFunction_tryCatch=function(x_in){
#  out=tryCatch({model_function(x_in)},
#               error=function(cond) {rep(NA, obsLength)}
#  )
#  return(out)
# }

#model_matrix_function=function(X){
#  Y_list=foreach(i=1:dim(X)[1], .export = c("model_function", "modelFunction_tryCatch"))%dopar%{ #make sure to include all related functions in .export and all used packages in .packages for more information read documentation of dopar
#      modelFunction_tryCatch((X[i,]))
#    }
  
#  Y=t(matrix(unlist(Y_list),ncol=length(Y_list)))
#}



## -----------------------------------------------------------------------------
CGNM_result=Cluster_Gauss_Newton_method(nonlinearFunction=model_matrix_function,
targetVector = observation,
initial_lowerRange =rep(0.01,3),initial_upperRange =  rep(100,3),lowerBound = rep(0,3), saveLog=TRUE, num_minimizersToFind = 500, ParameterNames = c("Ka","V1","CL"))


#stopCluster(cluster) #make sure to close the created cluster if needed

## -----------------------------------------------------------------------------
unlink("CGNM_log", recursive=TRUE)
unlink("CGNM_log_bootstrap", recursive=TRUE)

