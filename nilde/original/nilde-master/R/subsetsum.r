
 fn3 <- function(ss,a,n,M){ ## getting s_1 and checking for solution for a 0-1 subset sum problem, where s_i is either 0 or 1...
     ## ss has to be in order from 2 till 'l'
     res <- numeric(0)
     s1 <- (n-sum(ss*a[-1]))/a[1]
     if ((s1==0 | s1 == 1) && sum(c(s1,ss))<= M) 
        res <- c(s1,ss)
     res
 } 


 ## function to work through the loops for a '0-1 subset sum' problem...... 
 recursive.fn3 <- function(w,b,a,n,M){
     S <- rep(0,0)
     d <- 0
     if(length(b)) {
        for( i in seq(0,b[1]) ) {
            if (length(b) >1){
                d<-sum(w*a[length(a):(length(a)-length(w)+1)])+i*a[length(b)+1] 
                b[2] <- min(1,floor((n-d)/a[length(b)]))  
            }
            S <- c(S, Recall( c(w,i), b[-1],a,n,M))
        }
     } else {
          return(fn3(rev(w),a,n,M))
     }
     S 
 }

 fn4 <- function(ss,a,n,M,bounds){ ## getting s_1 and checking for solution for a bounded subset sum problem, where s_i is either 0 or 1...
     ## ss has to be in order from 2 till 'l'
     res <- numeric(0)
     s1 <- (n-sum(ss*a[-1]))/a[1]
     if (s1>=0 && s1 <= bounds[1] && s1==floor(s1) && sum(c(s1,ss))<= M) 
        res <- c(s1,ss)
     res
 } 

 ## function to work through the loops for a bounded subset sum problem...... 
 recursive.fn4 <- function(w,b,a,n,M,bounds){
     S <- rep(0,0)
     d <- 0
     if(length(b)) {
        for( i in seq(0,b[1]) ) {
            if (length(b) >1){
                d<-sum(w*a[length(a):(length(a)-length(w)+1)])+i*a[length(b)+1] 
                b[2] <- min(bounds[length(b)],floor((n-d)/a[length(b)]))  
            }
            S <- c(S, Recall( c(w,i), b[-1],a,n,M,bounds))
        }
     } else {
          return(fn4(rev(w),a,n,M,bounds))
     }
     S 
 }



get.subsetsum <- function(a,n,M=NULL,problem="subsetsum01", bounds=NULL){
## function to solve  0-1 or bounded subset sum problems:
## given as set of pos integers (a) and a pos integer n, does any non-empty subset sum to n?
## 'a' is an l-vector of positive integers with l>= 2
## 'n' is positive integer  
## 'M' is a positive integer, the maximum number of summands, M <= l
## 'problem' is one of the two problems to be solved: "subsetsum01" (default) for a 0-1 subset sum problem, 
##         or "bsubsetsum" a bounded subset sum problem, 
## 'bounds' is an l-vector of positive integers, bounds of s_i, i.e. 0 <= s_i <= b_i

 if (length(a) < 2) {stop("length of vector 'a' has to be more than 1")}
 if (!isTRUE(all(a == floor(a))) || !isTRUE(all(a > 0))) stop("'a' must only contain positive integer values")
 if (length(n) >1) {stop("'n' has to be a positive integer")}
 if (!isTRUE(n == floor(n)) || !isTRUE(n > 0)) {stop("'n' has to be a positive integer")}
 l <- length(a)
 if (is.null(M)) M <- floor(n/min(a))
 else { if (!isTRUE(M == floor(M)) || !isTRUE(M > 0)) {stop("'M' has to be a positive integer")}
       if (M > l) stop("'M' has to be less or equal to the length of 'a'")
   }

 if (!(problem %in% c("subsetsum01", "bsubsetsum")))  stop("unknown problem is used") 
 if (problem=="bsubsetsum" & is.null(bounds)) stop("no upper limits for the set of indices, 'bounds', supplied to solve the bounded problem") 
 if (problem=="bsubsetsum" & length(bounds)!=length(a)) stop("lengths of vectors 'bounds' and 'a' must be the same")
 
  ra <- rank(a, ties.method= "first")
  a <- sort(a)
  bounds <- bounds[ra]
  out <-numeric(0)

  if (problem=="subsetsum01"){
               b <- c(min(1,floor(n/a[l])),rep(NA,l-2)) 
               out <- recursive.fn3(numeric(0), b,a,n,M)
               if (length(out)==0) {out <- NULL }
               else {
                  dim(out) <- c(l,length(out)/l)
                  out <- as.matrix(out[ra,],l,length(out)/l) ## going back to original unsorted coefficients
                  rownames(out) <- paste("s", c(1:l), sep="")
                  colnames(out) <- paste(c("sol."), seq(1:dim(out)[2]), sep="")    
               } 
      } else if (problem=="bsubsetsum"){
               b <- c(min(bounds[l],floor(n/a[l])),rep(NA,l-2)) 
               out <- recursive.fn4(numeric(0), b,a,n,M,bounds)
               if (length(out)==0) {out <- NULL }
               else {
                 dim(out) <- c(l,length(out)/l)
                 out <- as.matrix(out[ra,],l,length(out)/l) ## going back to original unsorted coefficients
                 rownames(out) <- paste("s", c(1:l), sep="")
                 colnames(out) <- paste(c("sol."), seq(1:dim(out)[2]), sep="")    
               } 
       } 
  if (is.null(out))  n.sol <-0           
      else n.sol <- ncol(out)
           
  object <- list()
  object$p.n <- n.sol
  object$solutions <- out
  class(object)<-"subsetsum"
  object   
} ## end get.subsetsum


## printing subsetsum solutions...
print.subsetsum <- function (x,...) 
## default print function for "subsetsum" objects
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


























