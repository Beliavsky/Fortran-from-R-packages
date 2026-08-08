#' Global optimization of a blackbox function given a finite budget of noisy evaluations,
#' via the Stochastic-Simultaneous Optimistic Optimisation algorithm.
#' The deterministic-SOO method is available for noiseless observations.
#' The knowledge of the function's smoothness is not required.
#' @title StoSOO and SOO algorithms
#' @param par vector with length defining the dimensionality of the optimization problem.
#' Providing actual values of \code{par} is not necessary (\code{NA}s are just fine).
#' Included primarily for compatibility with \code{\link[stats]{optim}}.
#' @param fn scalar function to be minimized, with first argument to be optimized over.
#' @param lower,upper vectors of bounds on the variables.
#' @param nb_iter number of function evaluations allocated to optimization.
#' @param control list of control parameters:
#' \itemize{
#' \item verbose: verbosity level, either \code{0} (default), \code{1} or greater than \code{1},
#' \item type: either '\code{det}' for optimizing a deterministic function or '\code{sto}' for a stochastic one,
#' \item k_max: maximum number of evaluations per leaf (default: from analysis),
#' \item h_max: maximum depth of the tree (default: from analysis),
#' \item delta: confidence (default: \code{1/sqrt(nb_iter)} - from analysis),
#' \item light: set to \code{FALSE} to return the search tree,
#' \item max: if \code{TRUE}, performs maximization.
#' }
#' @param ... optional additional arguments to \code{fn}.
#' @return list with components:
#' \itemize{
#' \item \code{par} best set of parameters (for a stochastic function, it corresponds to the minimum reached over the deepest unexpanded node),
#' \item \code{value} value of \code{fn} at \code{par},
#' \item \code{tree} search tree built during the execution, not returned unless \code{control$light == TRUE}.
#' }
#' @details 
#' The optional \code{tree} element returned is a list, whose first element is the root node and the last element the deepest nodes.
#' A each level, \code{x} provides the center(s), one per row, whose corresponding bounds are given by \code{x_min} and \code{x_max}. Then:
#' \itemize{
#' \item \code{leaf} indicates if \code{x} is a leaf (1 if \code{TRUE});
#' \item \code{new} indicates if \code{x} has been sampled last;
#' \item \code{sums} gives the sum of values at \code{x};
#' \item \code{bs} is for the upper bounds at \code{x};
#' \item \code{ks} is the number of evaluations at \code{x};
#' \item \code{values} stores the values evaluated as they come (mostly useful in the deterministic case)
#' }
#' @export
#' @author M. Binois (translation in R code), M. Valko, A. Carpentier, R. Munos (Matlab code)
#' @references
#' R. Munos (2011), Optimistic optimization of deterministic functions without the knowledge of its smoothness,
#' \emph{NIPS}, 783-791. \cr \cr
#' M. Valko, A. Carpentier and R. Munos (2013), Stochastic Simultaneous Optimistic Optimization,
#' \emph{ICML}, 19-27 \url{https://inria.hal.science/hal-00789606}. Matlab code: \url{https://team.inria.fr/sequel/software/StoSOO/}. \cr \cr
#' P. Preux, R. Munos, M. Valko (2014), Bandits attack function optimization, \emph{IEEE Congress on Evolutionary Computation (CEC)}, 2245-2252.
#'
#'@examples
#' #------------------------------------------------------------
#' # Example 1 : Deterministic optimization with SOO
#' #------------------------------------------------------------
#' ## Define objective
#' fun1 <- function(x) return(-guirland(x))
#'
#' ## Optimization
#' Sol1 <- StoSOO(par = NA, fn = fun1, nb_iter = 1000, control = list(type = "det", verbose = 1))
#'
#' ## Display objective function and solution found
#' curve(fun1, n = 1001)
#' abline(v = Sol1$par, col = 'red')
#'
#' #------------------------------------------------------------
#' # Example 2 : Stochastic optimization with StoSOO
#' #------------------------------------------------------------
#' set.seed(42)
#'
#' ## Same objective function with uniform noise
#' fun2 <- function(x){return(fun1(x) + runif(1, min = -0.1, max = 0.1))}
#'
#' ## Optimization
#' Sol2 <- StoSOO(par = NA, fn = fun2, nb_iter = 1000, control = list(type = "sto", verbose = 1))
#'
#' ## Display solution
#' abline(v = Sol2$par, col = 'blue')
#'
#' #------------------------------------------------------------
#' # Example 3 : Stochastic optimization with StoSOO, 2-dimensional function
#' #------------------------------------------------------------
#'
#' set.seed(42)
#'
#' ## 2-dimensional noisy objective function, defined on [0, pi/4]^2
#' fun3 <- function(x){return(-sin1(x[1]) * sin1(1 - x[2]) + runif(1, min = -0.05, max = 0.05))}
#'
#' ## Optimizing
#' Sol3 <- StoSOO(par = rep(NA, 2), fn = fun3, upper = rep(pi/4, 2), nb_iter = 1000)
#'
#' ## Display solution
#' xgrid <- seq(0, pi/4, length.out = 101)
#' Xgrid <- expand.grid(xgrid, xgrid)
#' ref <- apply(Xgrid, 1, function(x){(-sin1(x[1]) * sin1(1 - x[2]))})
#' filled.contour(xgrid, xgrid, matrix(ref, 101), color.palette  = terrain.colors,
#' plot.axes = {axis(1); axis(2); points(Xgrid[which.min(ref),, drop = FALSE], pch = 21);
#'              points(Sol3$par[1],Sol3$par[2], pch = 13)})
#'
#'
#' #------------------------------------------------------------
#' # Example 4 : Deterministic optimization with StoSOO, 2-dimensional function with plots
#' #------------------------------------------------------------
#' set.seed(42)
#'
#' ## 2-dimensional noiseless objective function, defined on [0, pi/4]^2
#' fun4 <- function(x){return(-sin1(x[1]) * sin1(1 - x[2]))}
#'
#' ## Optimizing
#' Sol4 <- StoSOO(par = rep(NA, 2), fn = fun4, upper = rep(pi/4, 2), nb_iter = 1000,
#'   control = list(type = 'det', light = FALSE))
#'
#' ## Display solution
#' xgrid <- seq(0, pi/4, length.out = 101)
#' Xgrid <- expand.grid(xgrid, xgrid)
#' ref <- apply(Xgrid, 1, function(x){(-sin1(x[1]) * sin1(1 - x[2]))})
#' filled.contour(xgrid, xgrid, matrix(ref, 101), color.palette  = terrain.colors,
#' plot.axes = {axis(1); axis(2); plotStoSOO(Sol4, add = TRUE, upper = rep(pi/4, 2));
#'              points(Xgrid[which.min(ref),, drop = FALSE], pch = 21);
#'              points(Sol4$par[1], Sol4$par[2], pch = 13, col = 2)})
StoSOO <- function(par, fn, ..., lower = rep(0, length(par)), upper = rep(1, length(par)), nb_iter, control = list(verbose = 0, type = "sto", max = FALSE, light = TRUE)){
  ## Default values of the control.
  control$nb_iter <- nb_iter
  
  if(is.null(control$verbose)){
    control$verbose <- 0
  }
  verbose <- control$verbose
  
  if(is.null(control$k_max)){
    control$k_max <- ceiling(nb_iter/log(nb_iter)^3)
  }
  
  control$sample_when_created <- 1
  
  if(is.null(control$type)){
    control$type <- "sto"
  }
  
  if(is.null(control$delta)){
    control$delta <- 1/sqrt(nb_iter)
  }
  
  if(control$type == "det"){
    control$k_max <- 1
    if(is.null(control$h_max)){
      control$h_max <- ceiling(sqrt(nb_iter))
    }
  }
  
  if(is.null(control$h_max)){
    control$h_max <- ceiling(sqrt(nb_iter/control$k_max))
  }
  
  d <- length(par)
  UCBK <- log((control$nb_iter)^2/control$delta)/2
  
  if(is.null(control$max)) control$max <- FALSE
  if(is.null(control$light)) control$light <- TRUE
  
  if(control$max) fnscale <- 1
  else fnscale <- -1
  
  f <- function(par){
    fn(par * (upper - lower) + lower, ...)/fnscale
  }
  
  ## Initialisation of the tree
  struct_list <- list()
  struct_list$x_min <- numeric(0)
  struct_list$x_max <- numeric(0)
  struct_list$x <- numeric(0)
  struct_list$leaf <- numeric(0)
  struct_list$new <- numeric(0)
  struct_list$sums <- numeric(0)
  struct_list$bs <- numeric(0)
  struct_list$ks <- numeric(0)
  struct_list$values <- list()
  
  t <- rep(list(struct_list), control$h_max)
  
  t[[1]]$x_min <- matrix(0, 1, d)
  t[[1]]$x_max <- matrix(1, 1, d)
  t[[1]]$x <- matrix(0.5, 1, d)
  t[[1]]$leaf <- 1
  t[[1]]$new <- 0
  t[[1]]$sums <- f(t[[1]]$x)
  t[[1]]$ks <- 1
  t[[1]]$bs <- t[[1]]$sums + sqrt(UCBK)
  t[[1]]$values <- list(t[[1]]$sums)
  
  Xs <- t[[1]]$x
  ys <- t[[1]]$sums
  
  ## execution
  finaly <- -Inf # for deterministic case
  at_least_one <- 1 # at least one leaf was selected
  n <- 1
  while(n < control$nb_iter){
    if(at_least_one != 1)
      break
    if(verbose > 1){
      cat("----- new pass ", n, "of ", control$nb_iter, " evaluations used ..", "\n")
    }
    
    v_max <- -Inf
    at_least_one <- 0
    
    for(h in 1:control$h_max){ # traverse the whole tree, depth by depth
      if(n > control$nb_iter)
        break
      i_max <- -1
      b_hi_max <- -Inf
      
      if(!is.null(nrow(t[[h]]$x))){
        for(i in 1:nrow(t[[h]]$x)){ # find max UCB at depth h
          if(t[[h]]$leaf[i] == 1 & t[[h]]$new[i] == 0){
            if(control$type == "sto"){
              b_hi <- t[[h]]$bs[i]
            }else{
              b_hi <- t[[h]]$sums[i]/t[[h]]$ks[i]
            }
            if(b_hi > b_hi_max){
              b_hi_max <- b_hi
              i_max <- i
            }
          }
        }
      }
      
      
      # we found a maximum open the leaf (h,i_max)
      if(i_max > -1){
        if(verbose > 2)
          cat("max b-value for:", b_hi_max, "(", i_max, " of ", nrow(t[[h]]$x), ")..\n")
        
        # Animations (TODO)
        
        # check maximum depth constraint
        if((h + 1) > control$h_max){
          if(verbose > 0) cat("Attempt to go beyond maximum depth refused. \n")
        }else{
          if(b_hi_max >= v_max){
            at_least_one <- 1
            # Sample the state and collect the reward
            xx <- t[[h]]$x[i_max,]
            if(t[[h]]$ks[i_max] < control$k_max){
              sampled_value <- f(xx)
              Xs <- rbind(Xs, xx)
              ys <- c(ys, sampled_value)
              if(sampled_value > finaly){
                finalx <- xx
                finaly <- sampled_value
              }
              t[[h]]$values[[i_max]] <- c(t[[h]]$values[[i_max]], sampled_value) # just for tracing
              t[[h]]$sums[i_max] <- t[[h]]$sums[i_max] + sampled_value # sample the function at xx
              t[[h]]$ks[i_max] <- t[[h]]$ks[i_max] + 1 # increment the count
              t[[h]]$bs[i_max] <- t[[h]]$sums[i_max]/t[[h]]$ks[i_max] + sqrt(UCBK/t[[h]]$ks[i_max]) # update b
              
              n <- n + 1
              
              if(verbose > 0){
                cat(n, ": sampling (", h, ", ", i_max, "), for the ", t[[h]]$ks[i_max],
                    " time  (max = ", control$k_max, "), x =", xx, " f(x) = ", fnscale * sampled_value, "\n")
              }
            }else{
              
              # the leaf becomes an inner node
              t[[h]]$leaf[i_max] <- 0
              
              # we find the dimension to split, it will be the one with largest range
              splitd <- which.max(t[[h]]$x_max[i_max,] - t[[h]]$x_min[i_max,])
              x_g <- xx
              x_g[splitd] <- (5 * t[[h]]$x_min[i_max, splitd] + t[[h]]$x_max[i_max, splitd])/6
              x_d <- xx
              x_d[splitd] <- (t[[h]]$x_min[i_max, splitd] + 5 * t[[h]]$x_max[i_max, splitd])/6
              
              # splits the leaf of the tree, if dim >1, splits along the largest dimension
              
              # left node
              t[[h+1]]$x <- rbind(t[[h+1]]$x, as.vector(x_g))
              if(control$sample_when_created){
                sampled_value <- f(x_g)
                Xs <- rbind(Xs, x_g)
                ys <- c(ys, sampled_value)
                if(sampled_value > finaly){
                  finalx <- x_g
                  finaly <- sampled_value
                }
                t[[h + 1]]$ks <- c(t[[h + 1]]$ks, 1)
                t[[h + 1]]$sums <- c(t[[h + 1]]$sums, sampled_value)
                t[[h + 1]]$bs <- c(t[[h + 1]]$bs, sampled_value + sqrt(UCBK))
                t[[h + 1]]$values <- c(t[[h + 1]]$values, sampled_value)
                
                n <- n + 1
                
                if(verbose > 0){
                  cat(n, ": sampling (", h, ", ", i_max, "), for the ", t[[h]]$ks[i_max],
                      " time  (max = ", control$k_max, "), x =", xx, " f(x) = ", fnscale * sampled_value, "\n")
                }
              }else{
                # not sampled yet
                t[[h + 1]]$ks <- c(t[[h + 1]]$ks, 0)
                t[[h + 1]]$sums <- c(t[[h + 1]]$sums, 0)
                t[[h + 1]]$bs <- c(t[[h + 1]]$bs, Inf)
                # t[[h + 1]]$values <- c(t[[h + 1]]$values, sampled_value)
              }
              
              t[[h + 1]]$x_min <- rbind(t[[h + 1]]$x_min, t[[h]]$x_min[i_max,])
              newmax <- t[[h]]$x_max[i_max,]
              newmax[splitd] <- (2 * t[[h]]$x_min[i_max, splitd] + t[[h]]$x_max[i_max, splitd])/3
              t[[h + 1]]$x_max <- rbind(t[[h + 1]]$x_max, as.vector(newmax))
              t[[h + 1]]$leaf <- c(t[[h + 1]]$leaf, 1)
              t[[h + 1 ]]$new <- c(t[[h + 1]]$new, 1)
              
              # right node
              t[[h+1]]$x <- rbind(t[[h+1]]$x, as.vector(x_d))
              if(control$sample_when_created){
                sampled_value <- f(x_d)
                Xs <- rbind(Xs, x_d)
                ys <- c(ys, sampled_value)
                if(sampled_value > finaly){
                  finalx <- x_d
                  finaly <- sampled_value
                }
                t[[h + 1]]$ks <- c(t[[h + 1]]$ks, 1)
                t[[h + 1]]$sums <- c(t[[h + 1]]$sums, sampled_value)
                t[[h + 1]]$bs <- c(t[[h + 1]]$bs, sampled_value + sqrt(UCBK))
                t[[h + 1]]$values <- c(t[[h + 1]]$values, sampled_value)
                
                n <- n + 1
                
                if(verbose > 0){
                  cat(n, ": sampling (", h, ", ", i_max, "), for the ", t[[h]]$ks[i_max],
                      " time  (max = ", control$k_max, "), x =", xx, " f(x) = ", fnscale * sampled_value, "\n")
                }
              }else{
                # not sampled yet
                t[[h + 1]]$ks <- c(t[[h + 1]]$ks, 0)
                t[[h + 1]]$sums <- c(t[[h + 1]]$sums, 0)
                t[[h + 1]]$bs <- c(t[[h + 1]]$bs, Inf)
                # t[[h + 1]]$values <- c(t[[h + 1]]$values, sampled_value)
              }
              
              newmin <- t[[h]]$x_min[i_max,]
              newmin[splitd] <- (t[[h]]$x_min[i_max, splitd] + 2*t[[h]]$x_max[i_max, splitd])/3
              t[[h + 1]]$x_min <- rbind(t[[h + 1]]$x_min, as.vector(newmin))
              t[[h + 1]]$x_max <- rbind(t[[h + 1]]$x_max, t[[h]]$x_max[i_max,])
              t[[h + 1]]$leaf <- c(t[[h + 1]]$leaf, 1)
              t[[h + 1]]$new <- c(t[[h + 1]]$new, 1)
              
              # central node
              t[[h + 1]]$x <- rbind(t[[h + 1]]$x, as.vector(xx))
              t[[h + 1]]$ks <- c(t[[h + 1]]$ks, t[[h]]$ks[i_max])
              t[[h + 1]]$sums <- c(t[[h + 1]]$sums, t[[h]]$sums[i_max])
              t[[h + 1]]$bs <- c(t[[h + 1]]$bs, t[[h]]$bs[i_max])
              newmin <- t[[h]]$x_min[i_max,]
              newmax <- t[[h]]$x_max[i_max,]
              newmin[splitd] <- (2*t[[h]]$x_min[i_max, splitd] + t[[h]]$x_max[i_max, splitd])/3
              newmax[splitd] <- (t[[h]]$x_min[i_max, splitd] + 2*t[[h]]$x_max[i_max, splitd])/3
              t[[h + 1]]$x_min <- rbind(t[[h + 1]]$x_min, as.vector(newmin))
              t[[h + 1]]$x_max <- rbind(t[[h + 1]]$x_max, as.vector(newmax))
              t[[h + 1]]$new <- c(t[[h + 1]]$new, 1)
              t[[h + 1]]$leaf <- c(t[[h + 1]]$leaf, 1)
              t[[h + 1]]$values <- c(t[[h + 1]]$values, list(t[[h]]$values[[i_max]]))
              
              # set the max Bvalue and increment the number of iteration
              v_max <- b_hi_max
            }
          }
        }
      }
    }
    
    # mark old just created leafs as not new anymore
    for(h in 1:control$h_max){
      if(!is.null(nrow(t[[h]]$x)))
        t[[h]]$new <- rep(0, nrow(t[[h]]$x))
    }
  }
  
  # get the deepest unexpanded node (and among all of those, pick a maximum)
  if(control$type == "sto"){
    for(h in control$h_max:1){
      if(length(t[[h]]$leaf) == 0){
        t[[h]] <- NULL
      }else{
        final_idx <- which(t[[h]]$leaf == 0)
        if(length(final_idx) != 0){
          final_idx <- final_idx[which.max(t[[h]]$sums[final_idx])]
          finalx <- t[[h]]$x[final_idx,]
          finaly <- t[[h]]$sums[final_idx]/t[[h]]$ks[final_idx]
          break
        }
      }
    }
  }
  if(control$light)
    return(list(par = finalx * (upper - lower) + lower, value = fnscale * finaly))
  
  return(list(par = finalx * (upper - lower) + lower, value = fnscale * finaly, tree = t,
              Xs = Xs %*% diag(upper - lower) + matrix(lower, nrow(Xs), d, byrow = TRUE), ys = fnscale * ys))
}

#' Plot bivariate tree structure obtained when running \code{\link[OOR]{StoSOO}}
#' @title Plot 2-D \code{\link[OOR]{StoSOO}} result
#' @param sol outcome of running \code{\link[OOR]{StoSOO}}, with \code{control$light} set to \code{FALSE}
#' @param levels which levels to print. Default to all levels
#' @param add if \code{TRUE}, use existing plot
#' @param cpch \code{\link[graphics]{points}} \code{pch} code for the centers
#' @param lcols color at each level, or a single color for all levels (default)
#' @param lower,upper vectors of bounds on the variables.
#' @param ylim vector of bounds, required in the 1d case
#' @export
#' @importFrom graphics lines points
#' @examples 
#' #------------------------------------------------------------
#' # Example 1 : Deterministic optimization with SOO, 1-dimensional function
#' #------------------------------------------------------------
#' ## Define objective
#' fun1 <- function(x) return(-guirland(x))
#' 
#' ## Optimization
#' Sol1 <- StoSOO(par = NA, fn = fun1, nb_iter = 1000, 
#'   control = list(type = "det", verbose = 1, light = FALSE))
#' 
#' ## Display objective function and solution found
#' curve(fun1, n = 1001, ylim = c(-1.3, 0))
#' abline(v = Sol1$par, col = 'red')
#' plotStoSOO(Sol1, ylim = c(-1.3, -1.1), add = TRUE)
#' 
#' #------------------------------------------------------------
#' # Example 2 : Deterministic optimization with SOO, 2-dimensional function
#' #------------------------------------------------------------
#' set.seed(42)
#' 
#' ## 2-dimensional noiseless objective function, defined on [0, pi/4]^2
#' fun <- function(x){return(-sin1(x[1]) * sin1(1 - x[2]))}
#' 
#' ## Optimizing
#' Sol <- StoSOO(par = rep(NA, 2), fn = fun, upper = rep(pi/4, 2), nb_iter = 1000,
#'   control = list(type = 'det', light = FALSE))
#'   
#' ## Display solution
#' plotStoSOO(Sol, upper = rep(pi/4, 2))
plotStoSOO <- function(sol, lower = rep(0, length(sol$par)), upper = rep(1, length(sol$par)), levels = NULL,
                       add = FALSE, cpch = '.', lcols = 1, ylim = NULL){
  if(length(sol$par) > 2) stop(paste("Cannot plot results in dimension", length(sol$par)))
  if(is.null(sol$tree)) warning("No tree in sol, perhaps use light = TRUE in control.")
  
  if(is.null(levels)) levels <- 1:length(sol$tree)
  if(length(lcols) == 1) lcols <- rep(lcols, length(levels))
  
  # 1D case
  if(length(sol$par) == 1){
    if(is.null(ylim)) stop("ylim is needed to define the y-scale.")
    if(!add){
      plot(NA, xlim = c(lower[1], upper[1]), ylim = ylim, 
           xlab =  expression(x[1]), ylab = expression(f))
    }
    
    for(level in levels){
      if(is.null(dim(sol$tree[[level]]$x))) break
      points(sol$tree[[level]]$x[,1] * (upper - lower) + lower, rep(ylim[1] + (level-0.5) * (ylim[2] - ylim[1]) / length(sol$tree), length(sol$tree[[level]]$x)), pch = cpch)
      
      for(i in 1:nrow(sol$tree[[level]]$x)){
        lines(x = c(sol$tree[[level]]$x_min[i], sol$tree[[level]]$x_min[i]) * (upper[1] - lower[1]) + lower[1],
              y = ylim[1] + (level - 1):level * (ylim[2] - ylim[1]) / length(sol$tree), 
              col = lcols[level])
        lines(x = c(sol$tree[[level]]$x_max[i], sol$tree[[level]]$x_max[i]) * (upper[1] - lower[1]) + lower[1],
              y = ylim[1] + (level - 1):level * (ylim[2] - ylim[1]) / length(sol$tree), 
              col = lcols[level])
        lines(x = c(sol$tree[[level]]$x_min[i], sol$tree[[level]]$x_max[i]) * (upper[1] - lower[1]) + lower[1],
              y = rep(ylim[1] + (level) * (ylim[2] - ylim[1]) / length(sol$tree), 2),
              col = lcols[level])
      }
    }
  }
  
  # 2D case
  if(length(sol$par) == 2){
    if(!add){
      plot(NA, xlim = c(lower[1], upper[1]), ylim = c(lower[2], upper[2]), 
           xlab =  expression(x[1]), ylab = expression(x[2]))
    }
    
    for(level in levels){
      if(is.null(dim(sol$tree[[level]]$x))) break
      points(sol$tree[[level]]$x %*% diag(upper - lower) + matrix(lower, nrow(sol$tree[[level]]$x), 2, byrow = TRUE), 
             pch = cpch)
      
      for(i in 1:nrow(sol$tree[[level]]$x)){
        lines(x = c(sol$tree[[level]]$x_min[i, 1], sol$tree[[level]]$x_min[i, 1], sol$tree[[level]]$x_max[i,1],
                    sol$tree[[level]]$x_max[i, 1], sol$tree[[level]]$x_min[i, 1]) * (upper[1] - lower[1]) + lower[1],
              y = c(sol$tree[[level]]$x_min[i, 2], sol$tree[[level]]$x_max[i, 2], sol$tree[[level]]$x_max[i, 2],
                    sol$tree[[level]]$x_min[i, 2], sol$tree[[level]]$x_min[i, 2]) * (upper[2] - lower[2]) + lower[2],
              col = lcols[level])
      }
    }
  }
}




