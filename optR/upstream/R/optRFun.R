#' Function based optimization module
#' 
#' @description Function based optimization module
#' @param formula : Function to optimize
#' @param x0 : Initial Solution
#' @param iteration : Number of Iterations
#' @param method  : Method for solving the optimization
#' @param tol : Tolerance
#' @return optRFun : Optimal Solution class
#' @export
optRFun<-function(formula, x0, iteration=30, method=c("newtonrapson"), tol=1e-9){
  if(length(method)>1){
    method<-"newtonrapson"
  }
  
  if(method=="newtonrapson"){
    optRFun<-optRFun.newtonRapson(formula, x0, iteration=iteration, tol=tol)  
  } else
  {
    stop("Method not defined!!")
  }
  
  optRFun$method<-method
  optRFun$iteration<-iteration
  optRFun$tol<-tol
  return(optRFun)
}


#' Function based optimization module
#' 
#' @description Newton Rapshon based optimization
#' @param formula : Function to optimize
#' @param x0 : Initial Solution
#' @param iteration : Number of Iterations
#' @param tol : Tolerance
#' @return optRFun : Optimal Solution
#' @export
optRFun.newtonRapson<-function(formula, x0, iteration=30, tol=1e-9){
  optRFun<-list()
  optRFun$optsol<-newtonRapson(formula, x0, iteration=iteration, tol=tol)
  optRFun$x0<-x0
  return(optRFun)
}