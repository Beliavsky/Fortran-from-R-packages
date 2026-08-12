
########################################
## _v1.1-6: solving TSP using get.subsetsum()...
########################################

tsp_solver <- function(data, labels=NULL,cluster=0,upper_bound=NULL, lower_bound=NULL,  ## symmetric=FALSE,
              method="cheapest_insertion", no_go=NULL){
## function to solve a TSP using the following steps:
## 1.(Initialization) Solve an assignment problem to obtain an initial lower bound on the value of the optimal TSP solution. Apply heuristics to obtain an initial upper bound.
## 2.(Subproblem solution) Given the initial lower bound construct all 0-1 solutions to a linear Diophantine equation.
## 3.(Degree constraints check) Remove solutions  that do not satisfy the degree constraints. 
## 4.(Subtour elimination) Remove solutions that contain subtours by applying a new simple 
## subtour elimination routine. If there is a solution(s) that contains no subtours, it forms the optiomal solution(s): stop. Otherwise, increase the initial lower bound by one and go to step 2. Repeat until the upper bound is reached.
##---------------------------------------------------------
## 'data' is an n x n matrix of costs/distances of the TSP (with 0's or NAs on the main diagonal). Costs/distances of the unconnected edges must be supplied as NA.
## 'labels' is an n vector of optional city labels. If not given, labels are taken from 'data'.
## 'symmetrics' is logical indicating if TSP is symmetric or not (default).
## 'cluster': Degree constraints can be checked in parallel using parLapply from the parallel package.
## 'cluster' is either '0'(default) for no parallel computing to be used; or '1' for one less than
## the number of cores; or user-supplied cluster on which to do checking.(a cluster here can be some cores of a single machine).
## 'method' is heuristic method used in TSP::solve_TSP() (default: cheapest_insertion)
## 'upper_bound' is a positive integer, an upper bound of the tour length (cost function), if not supplied (default: NULL) heuristic solution is obtained using TSP::solve_TSP(data,method).
## 'lower_bound' is a positive integer, a lower possible value of the tour lenght (cost function); 
## if not supplied (default: NULL), obtained by solving a corresponding assignment problem
## using lpSolve::lp.assign(data) 
## 'no_go' is a suitably large value used in the distance/cost matrix to make related edges infeasible, 
## if NULL (default) no_go=max(data)*10^5. It can be set to Inf for TSP(). However, lpSolve() is very sensitive to too large values and can result in high values of the lower_bound.
##===========================================================================================
 
 if (!is.matrix(data))
        stop("'data' must be a matrix ")
 ## check if matrix is square...
 if(dim(data)[1] != dim(data)[2]) stop("tsp_solver() requires a square matrix")

 if(!is.null(labels)) dimnames(data) <- list(labels, labels)
 ## make sure we have labels
 if(is.null(dimnames(data)))
       dimnames(data) <- list(1:dim(data)[1], 1: dim(data)[1])
 if(is.null(colnames(data)))  colnames(data) <- rownames(data)
 if(is.null(rownames(data)))  rownames(data) <- colnames(data)

 ## make sure data is numeric
 mode(data) <- "numeric"

 if (!is.null(upper_bound)&& length(upper_bound) >1) 
    stop("'upper_bound' has to be a positive integer")
 if (!is.null(lower_bound)&& length(lower_bound) >1) 
    stop("'lower_bound' has to be a positive integer") 
 if (!is.null(upper_bound)&& upper_bound <= 0) 
    stop("'upper_bound' has to be a positive integer")
 if (!is.null(lower_bound)&& lower_bound <= 0) 
    stop("'lower_bound' has to be a positive integer") 
 if (!is.null(upper_bound)&& !is.null(lower_bound)&& (upper_bound< lower_bound))
    stop("'lower_bound' has to be equal to or less than 'upper_bound'") 

 
 if (cluster==1){
    if (parallel::detectCores()>1)  ## no point otherwise
        cl <- parallel::makeCluster(parallel::detectCores()-1)
      else cl <- NULL
 } else if (cluster==0){
     cl <- NULL
 } else if (!is.null(cl)&&inherits(cluster,"cluster")){
     cl <- cluster
 } else {
       warning("Supplied cluster is unknown - ignored.")
       cl <- NULL
   } 
  
 d <- data ## copy data
 if (is.null(no_go)) no_go <- max(d,na.rm=TRUE)*10^5

 ## 1.(Initialization) i) apply TSP::solve_TSP(data,method) to obtain heuristic solution/upper_bound...
 ##                    ii) apply lpSolve::lp.assign(data) to get a lower bound of the tour length by solving the corresponding assignment problem (AP)...
 if (is.null(upper_bound)){ ## upper_bound is not supplied
    diag(data) <- 0 ## solve_TSP() sets 0's for diagonal elements
    data[is.na(data)] <- no_go ## TSP doesn't allow NAs
    atsp <- TSP::ATSP(data)
    tour<- TSP::solve_TSP(atsp,method=method)
    upper_bound <- as.integer(TSP::tour_length(tour))
 }

 ## if (is.null(lower_bound)){ ## lower_bound is not supplied
 ##   ## get 'badly' rough lower possible value of the tour length...
 ##   if (isSymmetric.matrix(d)){ ##  symmetric
 ##         sym_d <- d
 ##         sym_d[lower.tri(sym_d)] <- NA
 ##         lower_bound <- sum(sort(sym_d)[1:n_cities])
 ##   } else ## for asymmetric tsp 
 ##         lower_bound <-  sum(sort(dv)[1:n_cities])
 ##  if (!is.null(L0)&& L0 < lower_bound)
 ##       stop("supplied 'L0' is too small") 
 ## }

 if (is.null(lower_bound)){ ## lower_bound is not supplied
     ## get a lower bound of the tour length by solving the corresponding assignment problem...
     diag(data) <- no_go ## set suitably large values as lp.assign() treats 0's as 0 cost
     data[is.na(data)] <- no_go
     lower_bound <- as.integer(lpSolve::lp.assign(data)$objval) ## since set as.double(objval) in lp.assign
     if (!is.null(upper_bound)&& as.integer(upper_bound)< as.integer(lower_bound))
          stop("'upper_bound' must be larger that lower_bound") 
 }

 object <- list()
 object$lower_bound <- lower_bound
 object$upper_bound <- upper_bound

 posle_sol <- vector(mode = "list", length = upper_bound-lower_bound+1) ## initializing an empty list of solutions
 posle_length <- rep(NA, upper_bound-lower_bound+1) ## initializing a vector of lengths

 diag(d) <- NA ## setting NA's for diagonal elements to be removed in a vector below
 dv <- c(d)
 ind_na <- is.na(dv)
 dv <- dv[!is.na(dv)]
 n_cities <- ncol(d) ## number of cities/nodes in TSP

 iter <- 0
 while (lower_bound <= upper_bound){ ## loop within possible values of the tour length L0; starts at lower_bound and moves upward until a full tour is found. In this case a tour is an optimal tour.
   iter <- iter +1
   ## 2.(Subproblem solution) given lower_bound, the AP solution, check it for optimality
   ## by solving corresponding subset sum problem...
   ## first convert matrix c to a (n^2-n) vector of summands of the linear Diophantine equations,
   ## moving column-wise... 
   g <- get.subsetsum(a=dv,M=n_cities,n=lower_bound,problem="subsetsum01")
   if (!is.null(g$solutions)){ ## there is at least one feasible solution... 
     ## removing solutions of subset sum problem with the number of 1's not equal to n_cities...
     g_sol <- as.matrix(g$solutions[, colSums(g$solutions)==n_cities])
    if (ncol(g_sol)==0) { ## no solutions...
         list_tours <- NULL
         tour_length <- NA        
    } else { ## there are feasible solutions satisfying sum to n_cities constraint 
     ## 3.(Degree constraints check) Remove solutions  that do not satisfy the degree constraints...
     ## first, creating list of solutions as vectors... 
     list_sol <- apply(g_sol,2,list)  
     rm(g_sol) ## save space 
     
     if (!is.null(cl)&&inherits(cl,"cluster")) { ## use parallel computation
       # cl <- makeCluster(detectCores()-1)
       ## then creating list of solutions as matrices...
       ## could also use makeForkCluster, but read warnings first!
       list_sol_mats <- parallel::parLapply(cl,list_sol,function(x,ind_na,n_cities) {
                                   x0 <- rep(0,n_cities*n_cities)
                                   x0[ind_na] <- NA
                                   x0[!ind_na] <- x[[1]]
                                   matrix(x0,n_cities,n_cities)
                                 }, ind_na, n_cities)
       ## checking column-wise degree constraints on each solution in the list...
       list_sol_dconst <- parallel::parLapply(cl,list_sol_mats, function(x,n_cities){
                               if (sum(colSums(x,na.rm=T)==1)!=n_cities || sum(rowSums(x,na.rm=T)==1)!=n_cities)
                                    x <- NULL ## setting NULL is the solution doesn't satisfy degree const  
                                else x                                    
                                  }, n_cities)
       ind <- which(parallel::parSapply(cl,list_sol_dconst, is.null))
       if (length(ind)!=0)
           list_sol_dconst <- list_sol_dconst[-which(parallel::parSapply(cl,list_sol_dconst, is.null))] ## removing NULLs from the list
     } else{ ## cl= NULL
        ## then creating list of solutions as matrices...
        list_sol_mats <- lapply(list_sol,function(x,ind_na,n_cities) {
                                   x0 <- rep(0,n_cities*n_cities)
                                   x0[ind_na] <- NA
                                   x0[!ind_na] <- x[[1]]
                                   matrix(x0,n_cities,n_cities)
                                 }, ind_na, n_cities)
        ## checking column-wise degree constraints on each solution in the list...
        list_sol_dconst <- lapply(list_sol_mats, function(x,n_cities){
                               if (sum(colSums(x,na.rm=T)==1)!=n_cities || sum(rowSums(x,na.rm=T)==1)!=n_cities)
                                    x <- NULL ## setting NULL is the solution doesn't satisfy degree const  
                               else x                                    
                                  }, n_cities)
        ind <- which(sapply(list_sol_dconst, is.null))
        if (length(ind)!=0)
          list_sol_dconst <- list_sol_dconst[-which(sapply(list_sol_dconst, is.null))] ## removing NULLs from the list
       }
       rm(list_sol_mats) # save space
   ## 4.(Subtour elimination) Remove solutions that contain subtours...
   
       tours_org <- function(x,n_cities){ 
        ## function to re-arrange nodes to obtain tours or subtours,
        ## and return NULL in case of a subtour, or the tour...
          i <- 1
          while (i < n_cities){
             v1 <- c((i+1):n_cities)
             ind <- match(x[i,2],x[v1,1],nomatch = 0)          
             if (ind!=0) { ## there is a connecting node
           ##  ind <- which(x[i,2]== x[v1,1])
           ##  if (length(ind)!=0) { ## there is a connecting node
                oo <- c(ind+i,v1[-ind])
                x[v1,] <- x[oo,]
                if (x[1,1]==x[i+1,2]){
                     i<-i+1; break
                 }
               # else i <- i+1
             } else {
                 # i<- i-1; 
                 break
               }
             i <- i+1
           }
          nu_tour <- i 
          if (nu_tour!=n_cities) x <- NULL ## setting NULL is the solution has a subtour
          else x  
      }  ## end tours_org  

   if (length(list_sol_dconst)!=0) { ## there are some feasible solutions to check from...
       ## first, extract indexes of the passing nodes (nodes with 1's) of each feasible tour/solution...  
       list_tours_init <- lapply(list_sol_dconst, function(x) which(x==1, arr.ind=TRUE))
       if (!is.null(cl)&&inherits(cl,"cluster")) { ## use parallel computation to re-arrange nodes to obtain tours or subtours, and return NULL in case of a subtour, or the tour... 
             list_tours <- parallel::parLapply(cl,list_tours_init,tours_org, n_cities)
        } else{ ## cl= NULL, no parallelization...
             list_tours <- lapply(list_tours_init,tours_org, n_cities)
          }
       nulls_check <- sapply(list_tours, is.null)
       ## list_tours <- list_tours[-which(sapply(list_tours, is.null))] ## removing NULLs from the list
       rm(list_tours_init)
       if (sum(nulls_check)!=0)
            list_tours <- list_tours[-which(nulls_check)] ## removing NULLs from the list
       if (length(list_tours)!=0){ ## there are optimal solutions
          list_tours <- lapply(list_tours,function(x) x[,-2]) ## removing 2nd column, leaving simply tours
          names(list_tours) <- paste("tour", 1:length(list_tours), sep=".")
          tour_length <- sum(d[cbind(list_tours[[1]], c(list_tours[[1]][2:n_cities],list_tours[[1]][1]))])
       } else {
           list_tours <- NULL
           tour_length <- NA
          }
    } else{ ## no feasible solutions found, i.e. there is no optimal solutions for the given L0...
         list_tours <- NULL
         tour_length <- NA
       }
     rm(list_sol_dconst)  # save space     

     } ## end ncol(g_sol)!=0

   } else{ ## is.null(g$solutions): no solutions for subset sum problem...
       list_tours <- NULL
       tour_length <- NA
     }
   posle_sol[[iter]] <- list_tours
   posle_length[iter] <- tour_length
   if (is.null(list_tours)) ## no solution for current tour length, increase the tour length, L0
       lower_bound <- lower_bound +1
   else 
       break ## solution found with current lower_bound, L0, being optimal
 } ## end while loop within [lower_bound, upper_bound] for the tour length

  if (!is.null(cl)) parallel::stopCluster(cl) 

 # tours <- posle_sol[[which(posle_length==min(posle_length,na.rm=TRUE))]]
 # tour_length <- min(posle_length,na.rm=TRUE))
 # object <- list()
  object$coming_solutions <- posle_sol
  object$coming_tour_lengths <- posle_length
  object$iter <- iter ## number of L0 gone through

  if (sum(is.na(posle_length))!=length(posle_length)) { ## there are optimal solutions
     object$tour <- posle_sol[[which(posle_length==min(posle_length,na.rm=TRUE))]]
     object$tour_length <- min(posle_length,na.rm=TRUE)
  } else { ## no solutions
      object$tour <- NULL
      object$tour_length <- NULL
    }
  
  class(object) <- "tsp_solver"
  object
} ## end tsp_solver



## issues/questions....
## i) is it reasonable enough to set max(data)*10^5 in the cost/distance matrix for unconnected edges?



## printing subsetsum solutions...
print.tsp_solver <- function (x,...) 
## default print function for "tsp_solver" objects
{
    cat("\n")
    if (is.null(x$tour)) cat("no solutions", "\n")
      else {
         cat("The optimal tour length: ", x$tour_length, "\n", sep = "")
         cat("\n Optimal tours:\n")
         print(x$tour,  ...)
         cat("\n")
      }
    invisible(x)
}


























