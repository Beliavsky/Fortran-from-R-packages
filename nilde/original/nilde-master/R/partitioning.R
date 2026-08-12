######################################
## partitioning...
######################################
fn1p <- function(l,n,M){ ## getting partitions for 'at most' case...
     ## 'l' is a (n-1)-vector
     res <- numeric(0)
     s <- rep(NA,n+1) ## initialize repetitions ('powers' of the partition)
     s[1] <- M-n+l[1]
     s[2] <- n-2*l[1] +l[2]
     ind <- c(3:(n-1))
     s[ind] <- l[ind-2]-2*l[ind-1]+l[ind]
     s[n] <- l[n-2]-2*l[n-1]
     s[n+1] <- l[n-1]    
     for (i in 1:(n+1)) 
          if (s[i]!=0)  res <- c(res, rep(i-1,s[i]))
     res
  } 
  
  ## function to work through the loops (for 'at most' case)...
  recursive.fn1p <- function(w,ll,lu,n,M){
     S <- rep(0,0)
     if(length(lu)) {
        for( i in seq(ll[1],lu[1]) ) {
            if (length(lu) ==(n-1)){
                ll[2] <- max(0,2*i-n)
                lu[2] <- floor((n-2)*i/(n-1))  
            } else if (length(lu) >2) {
                  nw <- length(w)
                  ll[2] <- max(0,2*i-w[nw])
                  lu[2] <- floor((n-nw-2)*i/(n-nw-1)) 
              } else if (length(lu) == 2) {
                    nw <- length(w)
                    ll[2] <- max(0,2*i-w[nw])
                    lu[2] <- floor(i/2)  
                }
            S <- c(S, Recall(c(w,i), ll[-1],lu[-1],n,M))
        }
     } else 
          return(fn1p(w,n,M))
     S 
  }

  ##------------
  fn2p <- function(l,n,M){ ## getting partitions for 'exactly' case...
     ## 'l' is a (n-M-1)-vector
     res <- numeric(0)
     s <- rep(NA,n-M+1) ## initialize repetitions ('powers' of the partition)
     s[1] <- 2*M-n+l[1]
     s[2] <- n-M-2*l[1] +l[2]
     ind <- c(3:(n-M-1))
     s[ind] <- l[ind-2]-2*l[ind-1]+l[ind]
     s[n-M] <- l[n-M-2]-2*l[n-M-1]
     s[n-M+1] <- l[n-M-1]
     for (i in 1:(n-M+1)) 
          if (s[i]!=0)  res <- c(res, rep(i,s[i]))
     res
  } 
  
  ## function to work through the loops (for 'exactly' case)...
  recursive.fn2p <- function(w,ll,lu,n,M){
     S <- rep(0,0)
     if(length(lu)) {
        for( i in seq(ll[1],lu[1]) ) {
            if (length(lu) ==(n-M-1)){
                ll[2] <- max(0,2*i-n+M)
                lu[2] <- max(0,i-1)  
            } else if (length(lu) >1) {
                  ll[2] <- max(0,2*i-w[length(w)])
                  lu[2] <- max(0,i-1) 
              } 
            S <- c(S, Recall(c(w,i), ll[-1],lu[-1],n,M))
        }
     } else 
          return(fn2p(w,n,M))
     S 
  }


get.partitions <- function(n, M, at.most=TRUE){
## function to enumerate additive partitions of integer 'n' on at most or exactly 'M' parts, M <= n
## 'n' is a positive integer
## 'M' is a positive integer, M <= n
## 'at.most' is logical with TRUE standing for partitioning of n into at most M parts and FALSE for exactly M parts  
  if (length(n) >1 || length(M) >1) {stop("'n' and 'M' have to be positive integers")}
  if (!isTRUE(n == floor(n)) || !isTRUE(n > 0)) {stop("'n' has to be a positive integer")}
  if (!isTRUE(M == floor(M)) || !isTRUE(M > 0)) {stop("'M' has to be a positive integer")}
  if (M > n) {stop("'M' has to be less or equal to 'n'")}
  if (!is.logical(at.most)) {stop("'at.most' must be TRUE or FALSE")}

  out <-numeric(0)
  if (at.most){ ## getting partitions of n into AT MOST M parts...
     L <- n
     lu <- c(floor((L-1)*n/L),rep(NA,L-2)) 
     ll <- c(n-M,rep(NA,L-2))
     out <- recursive.fn1p(numeric(0), ll,lu,n,M)
  } else { ## getting partitions of n into EXACTLY M parts...
     ll <- c(max(0,n-2*M),rep(NA,n-M-2)) 
     lu <- c(n-M-1,rep(NA,n-M-2))
     if (length(ll)==1){
         l1 <- c(ll:lu)
         out <- numeric(0)
         for (i in l1)
               out <- c(out, c(rep(1,2*M-n+i),rep(2,n-M-2*i),rep(3,i)))         
     } else
         out <- recursive.fn2p(numeric(0), ll,lu,n,M)
  }
  if (length(out)==0) 
       {out<- NULL; n.sol <-0 } ## print("no solutions")}
  else {
     dim(out) <- c(M,length(out)/M)
   #  out <- out[c(nrow(out):1),c(ncol(out):1)]
     out <- out[c(1:nrow(out)),c(1:ncol(out))]
     rownames(out) <- paste("a", c(1:M), sep="")
     colnames(out) <- paste(c("p."), seq(1:dim(out)[2]), sep="")
     n.sol <- ncol(out)
  }
 object <- list()
 object$p.n <- n.sol
 object$partitions <- out
 class(object)<-"partitions"
 object
} ## end get.partitions


###################################################################
## printing the results of "partitions"...    ##
###################################################################

print.partitions <- function (x,...) 
## default print function for "partitions" objects
{
    cat("\n")
    cat("The number of partitions: ", x$p.n, "\n", sep = "")
    cat("\nPartitions:\n")
    printCoefmat(x$partitions,  ...)
    cat("\n")
    invisible(x)
}





#########################
## loading functions...
##########################

print.nilde.version <- function()
{ library(help=nilde)$info[[1]] -> version
  version <- version[pmatch("Version",version)]
  um <- strsplit(version," ")[[1]]
  version <- um[nchar(um)>0][2]
  hello <- paste("This is nilde ",version,".",sep="")
  packageStartupMessage(hello)
}


.onAttach <- function(...) { 
  print.nilde.version()
}




