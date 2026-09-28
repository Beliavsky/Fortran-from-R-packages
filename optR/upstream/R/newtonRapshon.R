#' Function for Newton Rapson roots for given equations
#' 
#' newtonRapson function perform optimization
#' @param f     : function to optimize
#' @param x     : Initial Solution
#' @param iteration : Iterations
#' @param tol   : Tolerance
#' @return x    : optimal roots
newtonRapson<-function(f, x, iteration=30, tol=1e-9){
  for(i in 1:iteration){
    jacobResult<-jacobian(f, x)
    if(sqrt((jacobResult$f0%*%jacobResult$f0)/length(x)) < tol){
      return(x)
    }
    dx=optR(x = jacobResult$jacobianMatrix, y = -jacobResult$f0, method = c("gauss"))
    dx<-dx$beta
    x=x+dx
    if(sqrt(sum(dx*dx)) < tol*max(c(abs(x), 1))){
      return(x)
    }
  }
  print("Solution not converged...")
  return(x)
}


#' Function to evaluate jacobian matrix from functions
#' 
#' jacobian is function to determine the jacobian matrix for function f for input x
#' @param f     : function to optimize
#' @param x     : Initial Solution
#' @return jacobiabMatrix: Jacobian matrix
#' @return f0: Intial solution
#' @export
jacobian<-function(f, x){
  h<-1e-4
  n<-length(x)
  jac<-matrix(0,n,n)
  f0<-f(x)
  for(i in 1:n){
    temp<-x[i]
    x[i]<-temp+h
    f1<-f(x)
    x[i]<-temp
    jac[,i]<-(f1-f0)/h
  }
  return(list("jacobianMatrix"=jac, "f0"=f0))
}