
 fn1 <- function(ss,a,n,M){ ## getting s_1 and checking for solution for 'at most' case...
     ## ss has to be in order from 2 till 'l'
     res <- numeric(0)
     s1 <- (n-sum(ss*a[-1]))/a[1]
     if (s1>=0 && s1==floor(s1) && sum(c(s1,ss))<= M) 
        res <- c(s1,ss)
     res
 } 

 ## function to work through the loops for 'at most' case...... 
 recursive.fn1 <- function(w,b,a,n,M){
     S <- rep(0,0)
     d <- 0
     if(length(b)) {
        for( i in seq(0,b[1]) ) {
            if (length(b) >1){
                d<-sum(w*a[length(a):(length(a)-length(w)+1)])+i*a[length(b)+1] 
                b[2] <- floor((n-d)/a[length(b)])  
            }
            S <- c(S, Recall( c(w,i), b[-1],a,n,M))
        }
     } else {
          return(fn1(rev(w),a,n,M))
     }
     S 
 }
 
 fn2 <- function(ss,a,n,M){ ## getting s_{l-1} and s_l and checking for solution for 'exactly' case...
     ## ss has to be in order from 3 till 'l'
     res <- numeric(0)
     l <- length(a)
     s2 <- (n-M*a[1] - sum(ss[length(ss):1]*(a[c(l:3)]-a[1])))/(a[2]-a[1])
     if (s2>=0 && s2==floor(s2) && sum(c(s2,ss))<= M ) {
        s1 <- M- sum(c(s2,ss))
        res <- c(s1,s2,ss)
     }
     res
 } 

 ## function to work through the loops for 'exactly' case...... 
 recursive.fn2 <- function(w,b,a,n,M){
     S <- rep(0,0)
     d <- 0
     if(length(b)) {
       if (sum(is.finite(b))==length(b)) {
         for( i in seq(0,b[1]) ) {
            if (length(b) >1){
                d<-sum(w*(a[length(a):(length(a)-length(w)+1)]-a[1]))+i*(a[length(b)+1] -a[1])
                b[2] <- max(0, floor((n-M*a[1]-d)/(a[length(b)]-a[1])) ) 
            }
            S <- c(S, Recall( c(w,i), b[-1],a,n,M))
         }
       } else S <- c(S, Recall( c(w,0), b[-1],a,n,M)) # return(fn2(rev(w),a,n,M)) ?????? most probably not correct here!!!!
     } else {
          return(fn2(rev(w),a,n,M))
     }
     S 
 }


nlde <- function(a,n,M=NULL,at.most=TRUE,option=0){
## function to solve nonnegative linear diophantine equation (NLDE): 
## a_1*s_1 +a_2*s_2 +...+ a_l*s_l =n
## 'a' is an l-vector of positive integers (coefficients of the left-hand-side of NLDE) with l>= 2
## 'n' is positive integer of the right-hand-side of NLDE
## vector of s are unknowns to be found
## 'at.most' is logical with TRUE standing for constructing all partitions of n into at most M parts and FALSE for exactly M parts  
## 'M' is a positive integer, M <= n
## 'option' when set to '1' solves nlde 0-1 problem (nlde has only 0 -1 solutions),
##          when 'option=2' (or > 2) solves 0-1 NLD inequality

 if (length(a) < 2) {stop("length of vector 'a' has to be more than 1")}
 if (!isTRUE(all(a == floor(a))) || !isTRUE(all(a > 0))) stop("'a' must only contain positive integer values")
 if (length(n) >1) {stop("'n' has to be a positive integer")}
 if (!isTRUE(n == floor(n)) || !isTRUE(n > 0)) {stop("'n' has to be a positive integer")}
 if (!is.logical(at.most)) stop("'at.most' must be TRUE or FALSE")
 if (is.null(M)) M <- floor(n/min(a))
 else { if (!isTRUE(M == floor(M)) || !isTRUE(M > 0)) {stop("'M' has to be a positive integer")}
       if (M > n) stop("'M' has to be less or equal to 'n'")
   }
 if (option < 0) warning("'option' must be 0 or positive, ignored") 


  ra <- rank(a, ties.method= "first")
  a <- sort(a)
  l <- length(a)
  if (option > 1) { ## solving 01 inequality...
     a1 <- c(a[ra],1) ## adding a slack variable 
     ra <- rank(a1, ties.method= "first")
     M <- floor(n/min(a1))
     a <- sort(a1)
     l <- length(a1)   
  }

  out <-numeric(0)
  ## solving nlde...
  if (at.most){ ## getting partitions of n into AT MOST M parts...
             b <- c(floor(n/a[l]),rep(NA,l-2)) 
             out <- recursive.fn1(numeric(0), b,a,n,M)
    } else { ## getting partitions of n into EXACTLY M parts...
                b <- c(floor((n-M*a[1])/(a[l]-a[1])), rep(NA, l-3))
                out <- recursive.fn2(numeric(0), b,a,n,M)
           }

  if (length(out)==0) 
       {out<- NULL; n.sol <-0 } 
  else {
        dim(out) <- c(l,length(out)/l)
        out <- as.matrix(out[ra,],l,length(out)/l) ## going back to original unsorted coefficients
        if (option > 1)
              out <- out[1:(l-1),] ## remove the last row of slacks
        if (option > 0){
              check01 <- function(vec) all(vec== 0 | vec== 1)
              ind <- apply(out,2,check01)
              out <- as.matrix(out[,ind])
        } 
        if (length(out)==0) {
              out <- NULL; n.sol <-0           
          } else {
             rownames(out) <- paste("s", c(1:dim(out)[1]), sep="")
             colnames(out) <- paste(c("sol."), seq(1:dim(out)[2]), sep="")    
             n.sol <- ncol(out)
            }
      }
  object <- list()
  object$p.n <- n.sol
  object$solutions <- out
  class(object)<-"nlde"
  object
} ## end nlde



## printing nlde solutions...
print.nlde <- function (x,...) 
## default print function for "nlde" objects
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









