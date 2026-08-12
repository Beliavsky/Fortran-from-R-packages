
get.knapsack <- function(objective,a,n,problem="uknap", bounds=NULL){
## function to solve knapsack problem: unbounded, 0-1 or bounded problem
## 'a' is an l-vector of positive integers (coefficients of the left-hand-side of NLDE) with l>= 2
## 'n' is positive integer of the right-hand-side of NLDE
## vector of s are unknowns to be found  
## 'objective' is a vector of coefficients of the objective function to be maximized when solving knapsack problem 
## 'problem' is one of the following names of the problems to be solved: 
##          "uknap" for an unbounded knapsack problem,
##          "knap01" for a 0-1 knapsack problem, 
##          "bknap" a bounded knapsack problem,
## 'bounds' is an l-vector of positive integers, bounds of s_i, i.e. 0 <= s_i <= b_i

 if (length(a) < 2) {stop("length of vector 'a' has to be more than 1")}
 if (length(objective) < 2) {stop("length of vector 'objective' has to be more than 1")}
 if (!isTRUE(all(a == floor(a))) || !isTRUE(all(a > 0))) stop("'a' must only contain positive integer values")
 if (length(n) >1) {stop("'n' has to be a positive integer")}
 if (!isTRUE(n == floor(n)) || !isTRUE(n > 0)) {stop("'n' has to be a positive integer")}
 if (!(problem %in% c("uknap", "knap01", "bknap")))
                stop("unknown problem is used") 
  if (is.null(objective)) stop("no coefficients of the objective function, 'objective', supplied to solve the knapsack problem") 
  if (length(a)!= length(objective)) stop("length of vector 'a' must be equal to the length of 'objective'")
  if (problem=="bknap" & is.null(bounds)) stop("no upper limits for the set of indices, 'bounds', supplied to solve the bounded problem") 
  if (problem=="bknap" & length(bounds)!=length(a)) stop("lengths of vectors 'bounds' and 'a' must be the same")
 
  ra <- rank(a, ties.method= "first")
  a <- sort(a)
 # bounds <- bounds[ra]
  l <- length(a)
  out <-numeric(0)
  a1 <- c(a[ra],1) ## adding a slack variable 
  ra <- rank(a1, ties.method= "first")
  M <- floor(n/min(a1))
  a1 <- sort(a1)
  l1 <- length(a1)
  b <- c(floor(n/a1[l1]),rep(NA,l1-2))     
  out <- recursive.fn1(numeric(0), b,a1,n,M)
  if (length(out)==0) {
           out <- NULL; n.sol <-0          
      } else {
          dim(out) <- c(l1,length(out)/l1)
          out <- as.matrix(out[ra,],l1,length(out)/l1) ## going back to original unsorted coefficients
          out <- out[1:l,] ## remove the last row of slacks
         # out1 <- out
         # obj.values <- colSums(out1*objective)
         # ind <- which(obj.values==max(obj.values))
         # out <- as.matrix(out1[,ind])                  
          if (problem=="knap01"){
              check01 <- function(vec) all(vec== 0 | vec== 1)
              ind <- apply(out,2,check01)
              out <- as.matrix(out[,ind])
              obj.values <- colSums(out*objective)
              out <- as.matrix(out[,which(obj.values==max(obj.values))]) 
          } else if (problem=="bknap"){
               checkb <- function(vec,bounds) all(vec<= bounds)
               ind <- apply(out,2,checkb,bounds=bounds)
               out <- as.matrix(out[,ind])
               obj.values <- colSums(out*objective)
               out <- as.matrix(out[,which(obj.values==max(obj.values))]) 
             }
          if (length(out)==0) {
              out <- NULL; n.sol <-0           
          } else {
               rownames(out) <- paste("s", c(1:l), sep="")
               colnames(out) <- paste(c("sol."), seq(1:dim(out)[2]), sep="")
               n.sol <- ncol(out)
            }  
       }
  object <- list()
  object$p.n <- n.sol
  object$solutions <- out
  class(object)<-"knapsack"
  object
} ## end get.knapsack


## printing knapsack solutions...
print.knapsack <- function (x,...) 
## default print function for "knapsack" objects
{
    cat("\n")
    if (is.null(x$solutions)) cat("no solutions", "\n")
      else {
         cat("The number of solutions: ", x$p.n, "\n", sep = "")
         cat("\nSolutions:\n")
         printCoefmat(x$solutions,  ...)
         cat("\n")
      }
    invisible(x)
}



