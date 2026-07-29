## Demo "numerical_experiments"  (06-11-2020)
## By Erik Hintz, Marius Hofert and Christiane Lemieux 

## Numerical experiments for 'pnvmix()', 'dnvmix()'and 'fitnvmix()' ############

## Table of contents ###########################################################
##
## 1.    Helper functions to perform the experiments
## 
## 2.    Numerical experiments for 'pnvmix'
## 
## 3.    Numerical experiments for 'dnvmix'
## 
## 4.    Numerical experiments for 'fitnmvix'
## 
## 5.    Data analysis of 5 Stock portfolio
## 
## 6.    Plots
## 
################################################################################


## Load packages
library(nvmix) 
library(mvtnorm) # for comparison with pmvt()
library(qrng) # to generate sobol points 
library(sensitivity) # for sobol indices
library(RColorBrewer) # for colors
library(microbenchmark) # for accurate timing 
library(QRM) # for 'fit.mst()' (EM algorithm for MVT dist'n) and 'returns()'
library(qrmdata) # for the dataset SP500
library(zoo) # for working with time series
library(rugarch) # fitting ARMA-GARCH in data analysis
library(qrmtools) # for fit_ARMA_GARCH

## Defaults for non-interactive demo
doPLOT   <- TRUE # generate plots?
doPDF    <- FALSE # generate .pdfs? (ignored if doPLOT = FALSE)
doRUN    <- FALSE # run all experiments?
doSTORE  <- FALSE # store result arrays via 'save(...)'?

## Global variables for plotting
lwds <- c(1, 1.3, 1.8, 1.6, 1.3, 1.5) # lwd for lty = 'l'
# 'solid', 'dashed', 'dotted', 'dotdash', 'longdash', 'twodash'

## Ask user if experiments shall be re-performed
answer <- 
   readline(cat("Press 'Y' if all numerical experiments shall be re-run (~90 hrs) before plotting or", 
                " press any other key if plots shall be generated from the files in ../data.", sep="\n"))
if(answer == "Y" || answer == "y") doRUN <- TRUE 

## Ask user if data shall be stored after being generated
if(doRUN){
   answer <- 
      readline(cat("Press 'Y' if all numerical results shall be stored via 'save(...)' in the current directory or",
                   " press any other key otherwise.", sep="\n"))
   if(answer == "Y" || answer == "y") doSTORE <- TRUE 
}


## Load data if necessary
if(!doRUN){
   data("numerical_experiments_data", package = "nvmix")
   ## Grab individual datasets
   if(!exists("numerical_experiments_data")) error("Could not find the list 'numerical_experiments_data'")
   fit.estimated       <- numerical_experiments_data$fit.estimated
   fit.analytical      <- numerical_experiments_data$fit.analytical
   fitnvmix.results    <- numerical_experiments_data$fitnvmix.results
   qqplots             <- numerical_experiments_data$qqplots
   pnvmix.t.variances  <- numerical_experiments_data$pnvmix.t.variances
   pnvmix.t.sobolind   <- numerical_experiments_data$pnvmix.t.sobolind
   pnvmix.t.timing     <- numerical_experiments_data$pnvmix.t.timing
   tailprobs           <- numerical_experiments_data$tailprobs
   dnvmix.results      <- numerical_experiments_data$dnvmix.results
   pnvmix.abserrors    <- numerical_experiments_data$pnvmix.abserrors
} 



## 1. Helper functions to perform the experiments ##############################

## Global variables
lwds <- c(1, 1.3, 1.8, 1.6, 1.3, 1.5) # lwd for lty = 'l'
# 'solid', 'dashed', 'dotted', 'dotdash', 'longdash', 'twodash'
## 1.1  Experiments for 'pnvmix()'  ############################################

#' Title: Data generation for numerical experiments for 'pnvmix()': 
#'        Estimate absolute error as a function of total number of fun.evals
#' 
#' @param qmix either a (vector of) strings ("inverse.gamma" or "pareto") or
#'        a (list of) function(s) which have to be of the form 
#'        function(u, nu) that are then interpreted as the quantile
#'        function of the mixing variable 'W'
#' @param nu numeric vector of length(qmix) containing parameter value 'nu' of 
#'        the underlying mixing variable 'W'       
#' @param d dimension of the normal variance mixture, can be vector 
#' @param max.fun.evals vector of maximal number of function-evaluations 
#'        to be used in each setting
#' @param n_init corresponds to control$fun.eval[1], see ?get_set_param
#'        Can be a vector. 
#' @return Array with dimensions c(length(d), length(qmix), length(max.fun.evals),
#'         2, 2, length(rep)) containing estimated absolute error in each setting. 
#'         with methods "sobol" and "PRNG" and with/without preconditioning;
#'         see also dimnames of the return. 
#' @author Erik Hintz
#' 
pnvmix_testing_abserr <- function(qmix, nu, d, n, max.fun.evals, n_init = 2^6,
                                  plot = FALSE)
{
   start <- Sys.time()
   set.seed(271) # for reproducibility 
   names.qmix  <- as.character(1:length(qmix))
   ## Make 'qmix' a list
   if(is.function(qmix)){
      ## 'qmix' is one function => now becomes a list with one element 
      qmix <- list(qmix)
   } else if(is.vector(qmix)){ # includes also is.character(qmix)
      ## Obtain names, if available
      for(i in seq_along(qmix)){
         names.qmix[i] <- if(is.character(qmix[i])) qmix[i] else as.character(i)
      }
      qmix <- as.list(qmix)
   } 
   ## Settings for control$method and control$precond
   precond <-  c(TRUE, FALSE)
   method  <-  c("sobol", "PRNG")
   ## Result object to store estimated errors 
   pnvmix.abserrors <- array(0, dim = c(length(d), length(qmix),
                                        length(max.fun.evals), length(method),
                                        length(precond), n),
                             dimnames = list(d = d, qmix = names.qmix, 
                                             n = max.fun.evals, method = method,
                                             precond = precond, rep = 1:n))
   ## Perform the simulation:
   for(i in seq_along(d)){
      ## One progress bar for each dimension
      pb. <- txtProgressBar(max = n, style = 3)
      ## Current dimension
      dim. <- d[i]
      ## Generate 'n' random wishart matrices:
      Wish.mat <- rWishart(n, dim., diag(dim.))
      for(rep in 1:n){
         ## Sample (random) upper limit and a (random) correlation matrix
         upper <- runif(dim.) * sqrt(dim.) * 3
         scale <- as.matrix(nearPD(cov2cor(Wish.mat[,,rep]))$mat)
         for(j in seq_along(qmix)){ 
            for(k in seq_along(max.fun.evals)){ 
               for(l in seq_along(method)){ 
                  for(m in seq_along(precond)){
                     ## pStudent() calls pnvmix(.., qmix = "inverse.gamma", ..) 
                     pnvmix.abserrors[i, j, k, l, m, rep] <- attr(
                        pnvmix(upper = upper, qmix = qmix[[j]], scale = scale, 
                               df = nu[j], alpha = nu[j], nu = nu[j], control = 
                                  list(pnvmix.abstol = 0, 
                                       fun.eval = c(n_init, max.fun.evals[k]), 
                                       method = method[l], precond = precond[m],
                                       max.iter.rqmc = 1e8), 
                               verbose = FALSE), "abs. error") # don't suppress warnings 
                  } # for(m in seq_along(preconds))
               } # for(l in seq_along(methods))
            } # for(k in seq_along(max.fun.evals)) 
         } # for(j in seq_along(qmixs))
         setTxtProgressBar(pb., rep) # update progress bar
      } # for(rep in reps)
      close(pb.)
   } # for(i in seq_along(dims))
   ## Total duration of the experiment
   duration <- Sys.time() - start
   attr(pnvmix.abserrors, "duration") <- duration 
   ## Return 
   pnvmix.abserrors
}

#' Title: Plot results obtained by 'pnvmix_testing_abserr()'
#'
#' @param pnvmix.abserrors a 6-dimensional result array exactly as produced 
#'    by 'pnvmix_testing_abserr()' 
#' @param index.qmix index of 'qmix' (second dimension 'pnvmix.abserrors') to be used 
#' @param index.d index of 'd' (first dimension of 'pnvmix.abserrors') to be used
#' @return plots estimated absolute errors as a function of 'n' (= max number of fun.evals)
#'         for methods "sobol" and "PRNG", with and without preconditioning,
#'         including regression coefficients in the legend.
#'         Invisibly returns input array 'pnvmix.abserrors'
#' @author Erik Hintz
#' 
pnvmix_testing_abserr_plot <- function(pnvmix.abserrors, index.qmix = 1,
                                       index.dim = 1){
   ## Basic checking
   stopifnot(length(dim(pnvmix.abserrors)) == 6)
   ## Grab various quantities from the dimnames of the input array
   dimnames.pnvmix.abserrors <- dimnames(pnvmix.abserrors)
   d             <- dimnames.pnvmix.abserrors$d
   names.qmix    <- dimnames.pnvmix.abserrors$qmix
   max.fun.evals <- as.numeric(dimnames.pnvmix.abserrors$n)
   precond       <- c(TRUE, FALSE) 
   method        <- c("sobol", "PRNG") 
   ## Get *mean* absolute errors in each setting over all length(reps) repetitions
   mean.abs.errors <- array(0, dim = c(length(max.fun.evals), length(method),
                                       length(precond)),
                            dimnames = list(n = max.fun.evals, method = method,
                                            precond = precond))
   ## Colors 
   pal <- colorRampPalette(c("#000000", brewer.pal(8, name = "Dark2")[c(7, 3, 5)]))
   cols <- pal(4) 
   ## Each plot has mean estimated errors as 
   ## fct of n in the 4 settings PRNG / Sobol + precond/!precond 
   for(k in seq_along(max.fun.evals)){
      for(l in seq_along(method)){
         for(m in seq_along(precond)){
            ## Mean absolute errors in dimension d[index.dim] and for qmix[index.qmix]
            mean.abs.errors[k, l, m] <- 
               mean(as.matrix(pnvmix.abserrors[index.dim, index.qmix, k, l, m, ]))
         }
      }
   }
   mean.abs.errors <- cbind(mean.abs.errors[, method = "sobol", precond = "TRUE"],
                            mean.abs.errors[, method = "sobol", precond = "FALSE"],
                            mean.abs.errors[, method = "PRNG",  precond = "TRUE"],
                            mean.abs.errors[, method = "PRNG",  precond = "FALSE"]) 
   ## Names for 'mean.abs.errors' and also names for the legend in the plot
   nms <- c("Sobol with reordering", "Sobol  w/o  reordering",
            "PRNG with reordering", "PRNG w/o  reordering") 
   colnames(mean.abs.errors) <- nms
   ## Compute regression coefficients
   coeff <- apply(mean.abs.errors, 2, function(y) lm(log(y) ~ log(max.fun.evals))$coeff[2])
   names(coeff) <- nms
   ## Plot
   plot(NA, log = "xy", xlim = range(max.fun.evals), ylim = range(mean.abs.errors),
        xlab = "Number of function evaluations", ylab = "Estimated error")
   lgnd <- character(4)
   lwds <- c(1, 1.1, 1.7, 1.6, 1.45, 1.3)[1:4] # lwds[k] for 'lty = k'
   for(k in 1:4) {
      lines(max.fun.evals, mean.abs.errors[,k], col = cols[k], lty = k, 
            lwd = lwds[k])
      lgnd[k] <- paste0(nms[k]," (",round(coeff[k], 2),")")
   }
   legend("bottomleft", bty = "n", lty = rev(1:4), col = rev(cols), 
          lwd = rev(lwds), legend = rev(lgnd))
   mtext(paste0("Dimension ", d[index.dim]),  side = 4) # Dimension on the 'right' axis
   ## Return
   invisible(pnvmix.abserrors)
}

#' Title: Estimate variance of integrand with and without reordering
#' 
#' @param qmix see ?pnvmix
#' @param N Number of runs: How often shall variance be estimated?
#' @param n number of points to estimate the variance of the integrand
#' @param mindim in each run, dimension is sampled unformly from betw 'mindim' and 'maxdim'
#' @param maxdim see 'mindim' 
#' @param nu.lower in each run, the mixing parameter 'nu' of 'qmix'
#'        is sampled unformly from betw 'nu.lower' and 'nu.lower'
#' @param nu.upper see 'nu.lower'
#' @return (N, 4) matrix, each row consists of the estimated variance with/
#'         without reordering along with the dimension and 'nu' used in that run 
#' @author Erik Hintz 
#' 
precond_testing_variance  <- function(qmix = "inverse.gamma", N = 1e3, n = 1e4, 
                                      mindim = 5, maxdim = 500,
                                      nu.lower = 0.1, nu.upper = 5){
   start <- Sys.time() # record duration
   ## Result matrix: Column 1/2: Variance of integrand with/without preconditioning,
   ## Column 3/4: dimension/df (not needed for plotting, but can be interesting)   
   pnvmix.variances <- matrix(NA, ncol = 4, nrow = N)
   colnames(pnvmix.variances) <- c("VarWith", "VarWithout", "d", "nu")
   pb. <- txtProgressBar(max = N, style = 3)
   for(i in 1:N){ # in each run...
      ## Sample dimension
      dim. <- sample(mindim:maxdim, 1)
      ## Sample df, upper limit, scale matrix (wlog correlation matrix)
      nu       <- runif(1, 0.1, 5)
      upper    <- runif(dim.)*sqrt(dim.)*(3)
      Wish.mat <- rWishart(1, dim., diag(dim.))
      scale    <- as.matrix(nearPD(cov2cor(Wish.mat[,,1]))$mat)
      nu       <- runif(1, nu.lower, nu.upper)
      ## Matrix of uniforms to estimate Var(g(U)) 
      U <- matrix(runif( dim. * n), nrow = n)
      ## Estimate the variances and store the results:
      ## 'pnvmix.g()' returns c( mean(g(U)), var(g(U)) ) if 'return.all = FALSE'
      var_precond   <- nvmix:::pnvmix_g(U, qmix = qmix, upper = upper, 
                                        scale = scale, df = nu, alpha = nu, nu = nu,
                                        precond = TRUE, return.all = FALSE)[2] 
      var_noprecond <- nvmix:::pnvmix_g(U, qmix = qmix, upper = upper, 
                                        scale = scale, df = nu, alpha = nu, nu = nu,
                                        precond = FALSE, return.all = FALSE)[2] 
      pnvmix.variances[i, ] <- c(var_precond, var_noprecond, dim., nu)
      setTxtProgressBar(pb., i) # update progress bar
   } # for(i in 1:N)
   ## Store results
   duration <- Sys.time() - start
   attr(pnvmix.variances, "duration") <- duration 
   ## Return
   pnvmix.variances
} 

#' Title: Plot results obtained by 'precond_testing_variance()'
#' 
#' @param pnvmix.variances (N, 4) matrix as created by 'precond_testing_variance()'
#' @param scatterplot logical; if TRUE (default) scatterplot of variances 
#'        with/without reordering is produduces, otherwise a histogram of
#'        variance ratios
#' @return see 'scatterplot'. Additionally invisibly returns the input 
#' @author Erik Hintz 
#' 
#' Title: Plot results obtained by 'precond_testing_variance()'
#' 
#' @param pnvmix.variances (N, 4) matrix as created by 'precond_testing_variance()'
#' @param scatterplot logical; if TRUE (default) scatterplot of variances 
#'        with/without reordering is produduces, otherwise a histogram/density
#'        plot of variance ratios
#' @param density logical; only used when 'scatterplot = FALSE'. If true,
#'        a density plot of the variance ratios is generated, otherwise a 
#'        histogram truncated to 30. 
#' @return see 'scatterplot'. Additionally invisibly returns the input 
#' @author Erik Hintz 
#' 
precond_testing_variance_plot  <- function(pnvmix.variances, scatterplot = TRUE,
                                           density = TRUE){
   ## Only take first two columns containing 'var_precond' and 'var_noprecond'
   vars  <- pnvmix.variances[, 1:2, drop = FALSE]
   N     <- dim(vars)[1] 
   ## Remove NAs
   notNA <- rowSums(is.na(vars)) == 0
   vars  <- vars[notNA,, drop = FALSE] 
   ## Sort 'vars' according to the ordering of first column (with preconditioning)
   ordering <- order(vars[, 1])
   vars     <- vars[ordering, ]
   ## Grab colors
   pal <- colorRampPalette(c("#000000", brewer.pal(8, name = "Dark2")[c(7, 3, 5)]))
   cols <- pal(2) # two colors: with/without reordering
   if(scatterplot){ # produce scatterplot of ordered pairs (var_with, var_without)
      plot(NA, xlim = c(1, N), ylim = range(vars), 
           # xlab = "Run (ordered according to Var(g(U)))", 
           xlab = expression(paste("Run (ordered according to Var(", tilde(g), "(U)))")),
           ylab = "Estimated variance")
      for(i in 2:1){
         points(vars[,i], col = cols[i], pch = i)
      }
      legend('topleft', c("Without reordering",  "With reordering"), 
             col = rev(cols), pch = 2:1, bty = "n")
   } else { # produce histogram /density plot of variance ratios
      ## Variance ratios for the histogram:
      vars.ratios <- vars[, 2] / vars[, 1] # non-precond / precond
      if(any(is.na(vars.ratios))) vars.ratios <- vars.ratios[-which(is.na(vars.ratios))]
      if(!density){
         end.hist <- 100 # any ratio > end.hist is set to end.hist as ow plot too wide
         vars.ratios.hist <- vars.ratios
         vars.ratios.hist[vars.ratios.hist > end.hist] <- end.hist
         hist(vars.ratios.hist, breaks = 50, freq = FALSE, main = NULL, 
              xlab = "Estimated variance ratio without versus with reordering")
         abline(v = 1, col = 'red', lty = 1, lwd = 1)
      } else {
         dens <- density(1/vars.ratios)
         plot(dens$x, dens$y, type = 'l', axes = T, 
              xlab = "Estimated variance ratio with versus without reordering",
              ylab = "Density")
         # abline(v = 1, col = 'red', lty = 2) 
         legend("topright", expression(paste("Var(", tilde(g), "(U))" / Var(g(U)))),
                lty = 1, col = "black", bty = 'n')
      }
   }
   invisible(pnvmix.variances)
}


#' Title: Estimate sobol indices with/without reordering
#'
#' @param qmix see ?pnvnix()
#' @param d dimension of the unerlying normal variance mixture
#' @param nu parameter of 'qmix()' (eg degree-of-freedom parameter)
#' @param n sample size of the 3 design matrices needed in 'sobolowen()'
#' @param seeds vector of seeds to use to generate 'upper', 'loc' and 'scale'
#' @param original.seed seed to use to generate design matrices 
#' @return Result array of dimension c(length(seeds), 2, 2, d + 1) with total
#'         effect index and sobol index in each component, once with and once
#'         without re-ordering. See also dimnames of the return. 
#' @author Erik Hintz
pnvmix_estimate_sobolind <- function(qmix = "inverse.gamma", d = 10, nu = 1,
                                     n = 5e5, seeds = c(10575, 96865, 30367),
                                     original.seed = 271)
{    
   start <- Sys.time()
   ## Result object
   which.index <- c("Total", "First") # total and first order index
   preconds    <- c(TRUE, FALSE) 
   pnvmix.t.sobolind <- array(0, dim = c(length(seeds), length(preconds), 
                                         length(which.index), d + 1),
                              dimnames = list(seed = seeds, precond = preconds, 
                                              index = which.index, 
                                              d = c(1:d, "variance")))
   set.seed(original.seed) 
   ## Generate design matrices:
   X1 <- data.frame(matrix(runif( d * n), nrow = n))
   X2 <- data.frame(matrix(runif( d * n), nrow = n))
   X3 <- data.frame(matrix(runif( d * n), nrow = n))
   ## Note: the first 'dim' elements are used to store the corresponding indices,
   ## element 'dim + 1' stores the corresponding variance of the integrand 
   pb. <- txtProgressBar(max = length(seeds), style = 3)
   for(i in seq_along(seeds)){
      set.seed(seeds[i])
      ## Generate 'scale' and 'upper'
      Wish.mat <- rWishart(1, d, diag(d))
      scale    <- as.matrix(nearPD(cov2cor(Wish.mat[,,1]))$mat)
      upper    <- runif(d)*sqrt(d)*(3)
      ## Sensitivity analysis via 'sobolowen()':
      sens_precond <- 
         sobolowen(model = nvmix:::pnvmix_g, X1 = X1, X2 = X2, X3 = X3,
                   nboot = 0, qmix = qmix, df = nu, nu = nu, alpha = nu,
                   upper = upper, scale = scale, precond = TRUE, 
                   return.all = TRUE)
      sens_noprecond <- 
         sobolowen(model = nvmix:::pnvmix_g, X1 = X1, X2 = X2, X3 = X3,
                   nboot = 0, qmix = qmix, df = nu, nu = nu, alpha = nu,
                   upper = upper, scale = scale, precond = FALSE,
                   return.all = TRUE)
      ## Estimate variance of the integrand (Var(g(U)))
      vars <- 
         c(nvmix:::pnvmix_g(X1, upper = upper, scale = scale, precond = TRUE,
                            qmix = "inverse.gamma", df = nu, return.all = FALSE)[2],
           nvmix:::pnvmix_g(X1, upper = upper, scale = scale, precond = FALSE,
                            qmix = "inverse.gamma", df = nu, return.all = FALSE)[2])
      ## Store results in the array
      pnvmix.t.sobolind[i, 1, 1, ] <- c(sens_precond$T$original,   vars[1])
      pnvmix.t.sobolind[i, 1, 2, ] <- c(sens_precond$S$original,   vars[1])
      pnvmix.t.sobolind[i, 2, 1, ] <- c(sens_noprecond$T$original, vars[2])
      pnvmix.t.sobolind[i, 2, 2, ] <- c(sens_noprecond$S$original, vars[2])
      setTxtProgressBar(pb., i) # update progress bar
   }
   ## Store array
   duration <- Sys.time() - start
   attr(pnvmix.t.sobolind, "duration") <- duration 
   pnvmix.t.sobolind
}

#' Title: Plot results obtained by 'pnvmix_estimate_sobolind()'
#'
#' @param pnvmix.t.sobolind array as output by 'pnvmix_estimate_sobolind()'
#' @param index.seed index of 'seed' to be plotted (first dimension of
#'         'pnvmix.t.sobolind')
#' @return Plot of total/first order indices with and without reordering.
#'         Additionally invisibly returns input array
#' @author Erik Hintz   
pnvmix_estimate_sobolind_plot <- function(pnvmix.t.sobolind, index.seed = 1){
   ## Grab dimension, first order indices and total indices 
   d <- length(dimnames(pnvmix.t.sobolind)$d) - 1
   first.indices <- 
      cbind(pnvmix.t.sobolind[index.seed, precond = "TRUE",  index = "First", 1:d], 
            pnvmix.t.sobolind[index.seed, precond = "FALSE", index = "First", 1:d]) 
   total.indices <- 
      cbind(pnvmix.t.sobolind[index.seed, precond = "TRUE",  index = "Total", 1:d], 
            pnvmix.t.sobolind[index.seed, precond = "FALSE", index = "Total", 1:d])
   ## Grab Var(g(U)) with/without reordering (index does not have an effect on Var(g(U)))
   vars <- 
      c(pnvmix.t.sobolind[index.seed, precond = "TRUE",  index = "First", d+1],
        pnvmix.t.sobolind[index.seed, precond = "FALSE", index = "First", d+1])
   nms <- c("With reordering", "Without reordering")
   colnames(first.indices) <- nms
   colnames(total.indices) <- nms
   names(vars)             <- nms
   ## Get the variance explained by all subsets of up to order = 1
   sums <- colSums(first.indices)
   ## Prepare plot
   pal  <- colorRampPalette(c("#000000", brewer.pal(8, name = "Dark2")[c(7, 3, 5)]))
   cols <- pal(2) # colors
   def.par <- par(no.readonly = TRUE) # save default, for resetting...
   layout(matrix(1:2, nrow = 2))
   ## Left: First order indices 
   plot(NA, xlab = "", ylab = "First order index", ylim = range(first.indices), 
        xlim = c(0, d - 1))
   for(k in 1:2){
      points(0:(d-1), first.indices[, k], col = cols[k], pch = k, type = "b", lty = k)
   }
   legend("topright", bty = "n", lty = rev(1:2), col = rev(cols), legend = rev(nms), 
          pch = rev(1:2))
   ## Right: Total effect indices 
   plot(NA, xlab = "", ylab = "Total effect index", ylim = range(total.indices), 
        xlim = c(0, d - 1))
   for(k in 1:2){
      points(0:(d-1), total.indices[,k], col = cols[k], pch = k, type = "b", lty = k)
   }
   legend("topright", bty = "n", lty = rev(1:2), col = rev(cols), legend = rev(nms),
          pch = rev(1:2))
   ## Text under the pot
   mtext(paste("With/without reordering: Var(g(U)) =", 
               toString(round(vars[1], 5)), "/", toString(round(vars[2], 5)), 
               "; sum of first order indices =", toString(round(sums[1], 2)),
               "/",toString(round(sums[2], 2))), side = 1, line = -2, outer = TRUE)
   par(def.par)
   ## Return
   invisible(pnvmix.t.sobolind)
}

#' Title: Estimate CPU times needed for 'pmvt()' (from mvtnorm) and 'pStudent()' 
#' 
#' @param d vector of dimensions to be used
#' @param n number of runs in each dimension
#' @param rep number of repetitions for each run (reduce noise in microbenchmar())
#' @param tol absolute error tolerance for pnvmix()
#' @param df degree-of-freedom parameter; has to be integer (o.w. pmvt() won't work)
#' @return array with dimensions c(2, length(d), n) with run-times for either
#'         'pmvt()' or 'pStudent()'; see also dimnames 
#' @author Erik Hintz 
#' 
pnvmix_timing_mvt <- function(d, n, rep, tol = 1e-3, df = 2){ 
   start <- Sys.time()
   algs <- c("pmvt", "pStudent") # algorithms under consideration 
   ## Set up result array
   pnvmix.t.timing <- array(0, dim = c( length(algs), length(d), n),
                            dimnames = list(alg = algs, d = d, run = 1:n))
   ## Set up progress bar
   pb. <- txtProgressBar(max = length(d), style = 3)
   for(i in seq_along(d)){ # in each dimension 
      dim. <- d[i]
      Wish.mat <- rWishart(n, dim., diag(dim.))
      for(j in 1:n){ # in each run
         ## Generate 'upper' and 'scale'
         upper <- runif(dim.) * sqrt(dim.)*(3)
         scale <- as.matrix(nearPD(cov2cor(Wish.mat[,,j]))$mat)
         ## pmvt needs df = 0 instead of df = Inf 
         df.pStudent <- df
         df.pmvt     <- if(is.infinite(df.pStudent)) 0 else df.pStudent
         
         ## Call pmvt()
         if(dim. <= 1000){  # pmvt() can only handle dimensions <= 1,000
            t <- microbenchmark(pmvt(upper = upper, sigma = scale, df = df.pmvt, 
                                     abseps = tol, maxpts= 1e9), times = rep)
            pnvmix.t.timing[1, i, j] <- mean(t$time)
         }
         ## Call pStudent()
         t <- microbenchmark(pStudent(upper = upper, scale = scale, df = df.pStudent, 
                                      control = list(pnvmix.abstol = tol,
                                                     fun.eval = c(2^7, 1e18), 
                                                     max.iter.rqmc = 200)), times = rep)
         pnvmix.t.timing[2, i, j] <- mean(t$time)
      } # for(j in seq_along(runs))
      setTxtProgressBar(pb., i) # update progress bar
   }
   ## Save results
   duration <- Sys.time() - start
   attr(pnvmix.t.timing, "duration") <- duration 
   attr(pnvmix.t.timing, "df")       <- df
   attr(pnvmix.t.timing, "abstol")   <- tol
   pnvmix.t.timing
} 

#' Title: Plot results obtained by 'pnvmix_timing_mvt()'
#' 
#' @param pnvmix.t.timing 3-dimensional array as output by 'pnvmix_timing_mvt()'
#' @param ratio logical if ratio CPU(pmvt)/CPU(pStudent) is to be plotted
#' @return plots estimated CPUs as fct of 'd' for 'pmvt()' and 'pStudent()'
#' @author Erik Hintz
#' 
pnvmix_timing_mvt_plot <- function(pnvmix.t.timing, ratio = FALSE){
   ## Input checking 
   stopifnot(length(dim(pnvmix.t.timing)) == 3)
   ## Grab dimensions used from dimnames:
   dims       <- as.numeric(dimnames(pnvmix.t.timing)$d)
   dims.pmvt  <- dims[dims<=1000] # pmvt() only works for dim <= 1,000
   ## Calculate mean/max/min CPU in each dimension for either method:
   CPU.pStudent <- rbind(
      apply(pnvmix.t.timing[method = "pStudent",,], 1, mean),
      apply(pnvmix.t.timing[method = "pStudent",,], 1, max),
      apply(pnvmix.t.timing[method = "pStudent",,], 1, min))
   rownames(CPU.pStudent) <- c("mean", "max", "min")
   CPU.pmvt     <- rbind(
      apply(pnvmix.t.timing[method = "pmvt",,], 1, mean),
      apply(pnvmix.t.timing[method = "pmvt",,], 1, max),
      apply(pnvmix.t.timing[method = "pmvt",,], 1, min))
   rownames(CPU.pmvt) <- c("mean", "max", "min")
   if(!ratio){
      ## Prepare plot
      pal  <- colorRampPalette(c("#000000", brewer.pal(8, name = "Dark2")[c(7, 3, 5)]))
      cols <- pal(2) # colors
      pchs <- c(18, 20) 
      nms  <- c("pStudent()", "pmvt()") # for the legend
      plot(NA, xlab = "Dimension", ylab = "CPU (Nanosec)", xlim = range(dims),
           ylim = c(0, max( max(CPU.pmvt), max(CPU.pStudent))))
      ## Plot pStudent() results:
      for(m in 1:3){
         points(dims, CPU.pStudent[m, ], col = cols[1], # 'max' and 'min' with same lty
                type = if(m == 1) 'p' else 'l', lty = min(m, 2), pch = pchs[1]) 
      }
      ## Plot pmvt() results:
      for(m in 1:3){
         points(dims.pmvt, CPU.pmvt[m, seq_along(dims.pmvt)], col = cols[2], # 'max' and 'min' with same lty
                type = if(m == 1) 'p' else 'l', lty = min(m, 2), pch = pchs[2])
      }
      ## Legend
      legend("topleft", rev(nms), col = rev(cols), pch = rev(pchs), bty = 'n')
   } else {
      ## Plot ratio pmvt() / pStudent() 
      plot(dims.pmvt, 
           CPU.pmvt[2, seq_along(dims.pmvt)]/CPU.pStudent[2, seq_along(dims.pmvt)], 
           xlab = "Dimension", ylab = "CPU (pmvt) / CPU (pStudent)", type = 'l', 
           lty = 1)
   }
   invisible(pnvmix.t.timing)
}


## 1.2  Experiments for 'dnvmix()'  ############################################

#' Title: Data generation for numerical experiments for 'dnvmix()'
#' 
#' @param d dimension
#' @param n sample size (= number of evaluation points of 'dnvmix()')
#' @param qmix either a (vector of) strings ("inverse.gamma" or "pareto") or
#'        a (list of) function(s) which have to be of the form 
#'        function(u, nu) that are then interpreted as the quantile
#'        function of the mixing variable 'W'
#' @param nu.sample numeric vector of length(qmix); n points are sampled via 
#'        'rnvmix()' with parameter 'nu' of 'qmix()' set to 'nu.sample'
#' @param nu.dens numeric vector of length(qmix); parameter value of 'nu' at 
#'        which the density shall be evaluated
#' @param control passed to 'dnvmix()', see ?dnvmix() and ?get_set_param()
#' @param dnvmix.doAdapt logical if adaptive procedure is to be used; can be a vector
#' @param seed either NA (then ignored) or an integer seed which is set by
#'        'set.seed(seed)' each time before calling 'rnvmix()' or 'dnvmix()'
#' @param plot logical if a plot shall be produced 
#' @param verbose logical if warnings should be thrown. Defaults to FALSE
#'        as non-adaptive procedure will almost surely not meet tolerance              
#' @return Array with dimensions c( length(qmix), 2, n, 6) and attributes 
#'         'd', 'MVT', 'nu.sample' and 'nu.dens' is returned. Entries of the array
#'         are 'maha', 'estimated log-density', 'true log-density', 'pgreater', 
#'        'estimated error' and 'CPU'.
#' @author Erik Hintz
dnvmix_testing <- function(d = 10, n = 1000, qmix = "inverse.gamma",
                           loc = rep(0, d), scale = diag(d), nu.sample = 1, 
                           nu.dens = 4, control = list(), 
                           dnvmix.doAdapt = c(TRUE, FALSE), seed = NA, 
                           plot = FALSE, plot.title = FALSE, verbose = FALSE)
{
   doAdapts        <- dnvmix.doAdapt
   control.doAdapt <- get_set_param(control)
   ## Deal with 'doAdapt' and maximum iterations accordingly 
   if(any(!doAdapts)){
      control.noAdapt <- control.doAdapt
      control.noAdapt$dnvmix.doAdapt <- FALSE
      ## In this case 'dnvmix()' only uses 'control$dnvmix.max.iter.rqmc.pilot' 
      ## many iterations which defaults to 4 (< control$max.iter.rqmc)
      control.noAdapt$dnvmix.max.iter.rqmc.pilot <- 800
   } else {
      control.noAdapt <- control.doAdapt
   }
   ## Prepare 'qmix' as list 
   if(is.function(qmix)){
      ## 'qmix' is one function => now becomes a list with one element 
      qmix <- list(qmix)
   } else if(is.vector(qmix)){ # includes also 'is.character(qmix) = TRUE'
      qmix <- as.list(qmix)
   } 
   ## To record if analytical density is available 
   special.mix <- rep("none", length(qmix))
   qmix.       <- vector("list", length(qmix)) 
   names.qmix  <- as.character(1:length(qmix))
   ## 'qmix.' is a list of quantile functions => force estimated densities
   for(i in 1:length(qmix)){
      qmix.[[i]] <- if(is.character(qmix[[i]])){
         ## qmix[[i]] character => analytical weights/densities available
         switch(qmix[[i]],
                "inverse.gamma" = {
                   special.mix[i] <- "inverse.gamma"
                   names.qmix[i]  <- "inverse.gamma"
                   function(u, nu) 1 / qgamma(1 - u, shape = nu/2, rate = nu/2)
                },
                "pareto" = {
                   special.mix[i] <- "pareto"
                   names.qmix[i]  <- "pareto"
                   function(u, nu) (1-u)^(-1/nu)
                },
                stop("only 'inverse.gamma' and 'pareto' allowed when 'qmix' is a string"))
      } else if(is.function(qmix[[i]])){
         ## Otherwise no further knowledge about 'qmix[[i]]' 
         if(!is.null(names(qmix))) names.qmix[i] <- names(qmix[[i]])
         if(any(fromalArgs(qmix[[i]]) == "nu")){
            function(u, nu) qmix[[i]](u, nu, ...)
         } else {
            function(u) qmix[[i]](u, ...)
         }
      } else stop("'qmix' has to be (a vector) of type 'character' or (a list of type) 'function'")
   }
   ## Check provided 'loc' and 'scale'
   stopifnot(length(loc) == d, dim(scale) == rep(d, 2))
   ## Check provided 'nu.sample' and 'nu.dens'
   if(length(nu.sample) != length(qmix)) nu.sample <- 
      rep(nu.sample, length.out = length(qmix))
   if(length(nu.dens) != length(qmix)) nu.dens <- 
      rep(nu.dens, length.out = length(qmix))
   ## Define result array to store data
   dnvmix.results <- array(NA, dim = c(length(qmix), 2, n, 6),
                           dimnames = list(qmix = names.qmix, 
                                           doAdapt = c("TRUE", "FALSE"), 1:n, 
                                           c('maha', 'estimated log-density', 
                                             'true log-density', 'pgreater', 
                                             'estimated error', 'CPU')))
   ## Already set attributes
   attr(dnvmix.results, 'd')           <- d
   attr(dnvmix.results, 'special.mix') <- special.mix
   attr(dnvmix.results, 'doAdapt')     <- dnvmix.doAdapt
   attr(dnvmix.results, 'nu.sample')   <- nu.sample
   attr(dnvmix.results, 'nu.dens')     <- nu.dens
   ## Loop over 'qmix'
   for(i in seq_along(qmix)){
      if(!is.na(seed)) set.seed(seed)
      ## Sample from nvmix-dist'n at which log-density is evaluated
      x <- rnvmix(n, loc = loc, scale = scale, qmix = qmix.[[i]], nu = nu.sample[i])
      ## Get (squared) mahalanobis distances of 'x' from 'loc' wrt 'scale'
      maha <- mahalanobis(x, loc, scale)
      ## Sort maha distances as well as x according to maha distances 
      ord  <- order(maha)
      maha <- sort(maha)
      x    <- x[ord, ]
      ## Call 'dnvmix' and measure CPU used:
      if(any(doAdapts)){
         if(!is.na(seed)) set.seed(seed)
         CPUused.doAdapt <- 
            system.time(ldens.est.doAdapt 
                        <- dnvmix(x, qmix = qmix.[[i]], loc = loc, scale = scale, 
                                  log = TRUE, control = control.doAdapt, 
                                  nu = nu.dens[i], verbose = verbose))[1] 
         error.doAdapt <- attr(ldens.est.doAdapt, "abs. error")
      } else {
         CPUused.doAdapt   <- NA
         ldens.est.doAdapt <- rep(NA, n)
         error.doAdapt     <- rep(NA, n)
      }
      if(any(!doAdapts)){
         if(!is.na(seed)) set.seed(seed)
         CPUused.noAdapt <- 
            system.time(ldens.est.noAdapt 
                        <- dnvmix(x, qmix = qmix.[[i]], loc = loc, scale = scale, 
                                  log = TRUE, control = control.noAdapt, 
                                  nu = nu.dens[i], verbose = verbose))[1]
         error.noAdapt <- attr(ldens.est.noAdapt, "abs. error")
      } else {
         CPUused.noAdapt   <- NA
         ldens.est.noAdapt <- rep(NA, n)
         error.noAdapt     <- rep(NA, n)
      }
      ## Check if analytical density avaible 
      ldens.true <- if(!(special.mix[i] == "none")){
         ## 'special.mix' is either "inverse.gamma" or "pareto"
         dnvmix(x, qmix = special.mix[i], loc = loc, scale = scale, 
                alpha = nu.dens[i], df = nu.dens[i], log = TRUE)
      } else rep(NA, n) # no analytical density function available 
      ## Estimate P(M>maha) where M = (X-mu)^T Sigma^{-1} (X-mu) 
      pgreater <- pmax(pgammamix(maha, qmix = qmix.[[i]], d = d, 
                                 lower.tail = FALSE, nu = nu.dens[i],
                                 control = control, verbose = verbose),
                       .Machine$double.xmin)
      ## Store results in the array
      dnvmix.results[i, "TRUE",,]  <- 
         cbind(maha, ldens.est.doAdapt, ldens.true, pgreater, error.doAdapt, 
               rep(CPUused.doAdapt[1], n))
      dnvmix.results[i, "FALSE",,] <- 
         cbind(maha, ldens.est.noAdapt, ldens.true, pgreater, error.noAdapt, 
               rep(CPUused.noAdapt[1], n))
      
   }
   ## Return
   if(plot) dnvmix_testing_plot(dnvmix.results) else dnvmix.results 
}

#' Title: Plot results of numerical experiments for 'dnvmix()'
#'
#' @param dnvmix.results a 4-dimensional result array exactly as produced 
#'        by 'dnvmix_testing()' 
#' @param index.qmix index of 'qmix' (first dimension 'dnvmix.results') to be used 
#' @return plots estimated log-densiites (with/without adaptive procedure and,
#'         if avail. true log-density) along with P( (X-mu)^T Sigma^{-1}(X-mu) > m) 
#'         where 'm' are the sampled mahalanobis distances and the run-time in the legend.
#'         Invisibly returns the input array 'dnvmix.results'
#' @author Erik Hintz
dnvmix_testing_plot <- function(dnvmix.results, index.qmix, plot.title = FALSE)
{
   i <- index.qmix # shorter 
   ## Basic input checking
   stopifnot(length(dim(dnvmix.results)) == 4, dim(dnvmix.results)[4] == 6)
   ## Grab sample size
   n <- dim(dnvmix.results)[3]
   ## Get names of 'qmix' used
   names.qmix <- unlist((dimnames(dnvmix.results)[1]), use.names = FALSE)
   ## Get some more quantities for plotting
   d       <- attr(dnvmix.results, "d")
   nu.dens <- attr(dnvmix.results, "nu.dens")
   ## Prepare colors for plotting
   pal  <- colorRampPalette(c("#000000", brewer.pal(8, name = "Dark2")[c(7, 3, 5)]))
   cols <- pal(4) # colors
   ## Grab (squared) mahalanobis distances
   maha <- dnvmix.results[i, "TRUE", , 1]
   ## Grab estimated/true log-densities with/without adaptive procedure
   ## as variables for better readability
   ldens.est.doAdapt    <- dnvmix.results[i, "TRUE", , 2]
   ldens.est.noAdapt    <- dnvmix.results[i, "FALSE",, 2]
   error.doAdapt        <- dnvmix.results[i, "TRUE",, 5]
   error.noAdapt        <- dnvmix.results[i, "FALSE",, 5]
   CPUused.doAdapt      <- dnvmix.results[i, "TRUE", 1, 6]
   CPUused.noAdapt      <- dnvmix.results[i, "FALSE", 1, 6]
   ldens.true           <- dnvmix.results[i, "TRUE", ,3]
   pgreater             <- dnvmix.results[i, "TRUE",,4]
   ## Prepare plot 
   rgX  <- range(sqrt(maha))
   rgY  <- range(ldens.true, ldens.est.doAdapt, na.rm = TRUE)
   def.par <- par(no.readonly = TRUE) # save default, for resetting...
   par(mar = c(4, 3, 3, 3) + 0.15)
   ## Initiallize legend (also lty, pch etc)
   lgnd <- expression(paste("P(", (X-mu)^T, Sigma^-1, (X-mu), ">", m^2,")"))
   ## Any errors NA? Will be set below
   anyerrorNA <- FALSE
   ## Is there a 'no.adapt' data-set provided?
   plot.noadapt <- !prod(is.na(ldens.est.noAdapt))
   ## Is there a 'do.adapt' data-set provided?
   plot.doadapt <- !prod(is.na(ldens.est.doAdapt))
   ## Plot P((X-mu)^T Sigma^{-1} (X-mu) > maha) 
   plot(sqrt(maha), pgreater, type = 'l', col = cols[1], xlab = "", 
        ylab = "", log = "y", axes = F, lty = 3, lwd = lwds[3],
        ylim = c(0.01*min(pgreater), 1))
   axis(2, ylim = range(pgreater), lwd = 1)
   # mtext(2, text = expression(paste("P(", D^2, ">", m^2,")")), line = 1.9)
   mtext(2, text = 
            expression(paste("P(", (X-mu)^T, Sigma^-1, (X-mu), ">", m^2,")")), 
         line = 1.9)
   lty.used <- c(3)
   lwd.used <- lwds[3]
   pch.used <- c(NA)
   col.used <- c(cols[1])
   par(new = T)
   ## Plot 'ldens.est.doadapt' as a function of 'sqrt(maha)'
   if(plot.doadapt){
      plot(sqrt(maha), ldens.est.doAdapt,  col = cols[2], lty = 1, xlab = "", 
           ylab = "", lwd = lwds[1], main = if(plot.title) 
              paste0("d = ", d, " ; n = ", n, " ; nu = ", nu.dens) else "", 
           xlim = rgX, ylim = rgY, axes = F, type = "l")
      ## Add to legend
      lgnd <- c(lgnd, paste0("est. log-density (adaptive, ", round(CPUused.doAdapt , 2), " sec)"))
      lty.used <- c(lty.used, 1)
      lwd.used <- c(lwd.used, lwds[1])
      pch.used <- c(pch.used, NA)
      col.used <- c(col.used, cols[2])
      ## Mark points for which the estimation is unreliable (ie 'error = NA')
      #if(any(is.na(error.doAdapt))){
      #   anyerrorNA <- TRUE
      #   whichNA <- which(is.na(error.doAdapt)) 
      #   points(sqrt(maha)[whichNA], ldens.est.doAdapt[whichNA], col = cols[2], pch = 9)
      #}
      ## Plot 'ldens.est.noadapt' as a function of 'sqrt(maha)'
      if(plot.noadapt){
         lines(sqrt(maha), ldens.est.noAdapt,  col = cols[3], lty = 2,
               lwd = lwds[2])
         ## Add to legend
         lgnd <- c(lgnd, paste0("est. log-density (non-adaptive, ", 
                                round(CPUused.noAdapt , 2), " sec)"))
         lty.used <- c(lty.used, 2)
         lwd.used <- c(lwd.used, lwds[2])
         pch.used <- c(pch.used, NA)
         col.used <- c(col.used, cols[3])
         ## Mark points for which the estimation is unreliable (ie 'error = NA')
         #if(any(is.na(error.noAdapt))){
         #   anyerrorNA <- TRUE
         #   whichNA <- which(is.na(error.noAdapt)) 
         #   points(sqrt(maha)[whichNA], ldens.est.noAdapt[whichNA], col = cols[3], pch = 9)
         #}
      }
   } else if(plot.noadapt){
      ## Plot *only* 'ldens.est.noadapt' as a function of 'sqrt(maha)'
      plot(sqrt(maha), ldens.est.noAdapt,  col = cols[3], lwd = lwds[1], 
           xlab = "", ylab = "", main = if(plot.title) 
              paste0("d = ", d, " ; n = ", n, " ; nu = ", nu.dens) else "", 
           xlim = rgX, ylim = rgY, axes = F, type = "l")
      lwd.used <- c(lwd.used, lwds[1])
      lty.used <- c(lty.used, 1)
      pch.used <- c(pch.used, NA)
      col.used <- c(col.used, cols[3])
      ## Add to legend
      lgnd <- c(lgnd, paste0("est. log-density (non-adaptive, ", round(CPUused.noAdapt , 2), " sec)"))
      ## Mark points for which the estimation is unreliable (ie 'error = NA')
      if(any(is.na(error.noAdapt))){
         anyerrorNA <- TRUE
         whichNA <- which(is.na(error.noAdapt)) 
         points(sqrt(maha)[whichNA], ldens.est.noAdapt[whichNA], col = cols[3], pch = 9)
      }
   } else {
      stop("Nothing to plot")
   }
   ## Also plot 'ldens.true' if provided
   if(!prod(is.na(ldens.true))){
      ## special.mix is either "inverse.gamma" or "pareto" => ldens.true available
      lines(sqrt(maha), ldens.true, type = "p", col = cols[4], lty = 2, pch = 4)
      ## Add to legend
      lgnd <- c(lgnd, "True log-density")
      lty.used <- c(lty.used, NA)
      lwd.used <- c(lwd.used, NA)
      pch.used <- c(pch.used, 4)
      col.used <- c(col.used, cols[4])
   } 
   ## Get both axis
   axis(4, ylim = rgY, lwd = 1, line = 0)
   mtext(4, text = "log-density", line = 2)
   axis(1, pretty(range(sqrt(maha)), 10))
   mtext("m (Mahalanobis distance)", side = 1, col = "black", line = 2)
   ## Plot legend
   ## Include 'errorNA' if necessary
   if(anyerrorNA){
      lgnd     <- c(lgnd, "Error NA")
      lty.used <- c(lty.used, NA)
      lwd.used <- c(lwd.used, NA)
      pch.used <- c(pch.used, 9)
      col.used <- c(col.used, "black")
   }
   legend('topright', lgnd, col = col.used, lty = lty.used, pch = pch.used,
          box.lty = 0, lwd = lwd.used)
   par(def.par)
   ## Return (invisbly)
   invisible(dnvmix.results)
}

## 1.3  Experiments for 'fitnvmix()'  ##########################################

#' Title: Data generation for numerical experiments for 'fitnvmix()'
#'
#' @param qmix either a (vector of) strings ("inverse.gamma" or "pareto") or
#'             a (list of) function(s) which have to be of the form 
#'             function(u, nu) that are then interpreted as the quantile
#'             function of the mixing variable 'W'
#' @param n    numeric vector. How large should the samples to be estimated from be?
#' @param d    numeric vector. In which dimension(s) shall experiment be performed?
#' @param nu   numeric vector of the same lenght as 'qmix' including parameter 
#'             values of 'nu' to be taken as true parameter
#' @param mix.param.bounds see ?fitnvmix()
#' @param nu.init see ?fitnvmix()
#' @param size.subsample vector of the same length as 'n', see ?fitnvmix()
#' @param control see ?fitnvmix()
#' @param analytical.weights shall results for analytical weights/densities also be 
#'         calculated if available? (only if qmix = "inverse.gamma" or "pareto")
#' @param use.EM.mvt shall QRM:fit.mst() (EM algorithm for MVT()) also be used?
#'         (only if 'qmix == "inverse.gamma"')
#' @param verbose see ?fitnvmix()
#' @param plot logical, if TRUE fitnvmix_testing_plot is called
#' @param ... additional parameters passed to 'qmix' if it is a function 
#' @author Erik Hintz
#' @return Result array c( length(dim), length(n), length(qmix), control$ECME.maxiter + 3),
#' including 'nu' estimates for each iteration and additional information,
#' see dimnames of the output for more details. 
fitnvmix_testing <- function(qmix = "inverse.gamma", n = 50, d = 10, nu = 2.5, 
                             mix.param.bounds = c(0.5, 9), nu.init = NA, 
                             size.subsample = n, control = list(), 
                             analytical.weights = TRUE, use.EM.mvt = TRUE, 
                             verbose = FALSE, plot = FALSE, ...)
{
   ## Record duration of experiment
   start <- Sys.time()
   stopifnot(n >= 2, d >= 1, is.list(control))
   control <- get_set_param(control)
   control$addReturns <- TRUE # force additional returns in 'fitnvmvix()' 
   if(length(size.subsample) != length(n)) size.subsample <- 
      rep(size.subsample, length.out = length(n))
   ## Prepare 'qmix' as list 
   if(is.function(qmix)){
      ## 'qmix' is one function => now becomes a list with one element 
      qmix <- list(qmix)
   } else if(is.vector(qmix)){ # includes also is.character(qmix)
      qmix <- as.list(qmix)
   } 
   ## To record if analytical weights are available 
   special.mix <- rep("none", length(qmix))
   qmix.       <- vector("list", length(qmix)) 
   names.qmix  <- as.character(1:length(qmix))
   ## 'qmix.' is a list of quantile functions => force estimated densities/weights
   for(i in 1:length(qmix)){
      qmix.[[i]] <- if(is.character(qmix[[i]])){
         ## qmix[[i]] character => analytical weights/densities available
         switch(qmix[[i]],
                "inverse.gamma" = {
                   special.mix[i] <- "inverse.gamma"
                   names.qmix[i]  <- "inverse.gamma"
                   function(u, nu) 1 / qgamma(1 - u, shape = nu/2, rate = nu/2)
                },
                "pareto" = {
                   special.mix[i] <- "pareto"
                   names.qmix[i]  <- "pareto"
                   function(u, nu) (1-u)^(-1/nu)
                },
                stop("only 'inverse.gamma' and 'pareto' allowed when 'qmix' is a string"))
      } else if(is.function(qmix[[i]])){
         ## Otherwise no further knowledge about 'qmix[[i]]' 
         if(!is.null(names(qmix))) names.qmix[i] <- names(qmix[[i]])
         if(any(fromalArgs(qmix[[i]]) == "nu")){
            function(u, nu) qmix[[i]](u, nu, ...)
         } else {
            ## For instance estimation of 'loc' and 'scale' for multivariate normal
            ## or multivariate-t with known degree-of-freedom parameter; experiment
            ## however focusses on estimation of 'nu' 
            warning("'qmix()' has no mixing parameter, estimation for 'loc' and 'scale' only")
            function(u) qmix[[i]](u, ...)
         }
      } else stop("'qmix' has to be (a vector) of type 'character' or (a list of) 'function'")
   }
   ## Result array
   fitnvmix.results <- 
      array(NA, dim = c( length(d), length(n), length(qmix), 
                         control$ECME.maxiter + 3),
            dimnames = list(d = d, n = n, qmix = names.qmix, 
                            c(1:(control$ECME.maxiter+1), "Analytical", "CPU")))
   ## Set up progress bar
   pb. <- txtProgressBar(max = length(d)*length(n), style = 3)
   for(i in seq_along(d)){ # for each dimenion
      ## Sample 'loc' and 'scale'
      loc       <- rep(0, d[i])
      diag_var  <- diag(runif(d[i], min = 2, max = 5))
      rWish     <- rWishart(1, d[i], diag(d[i]))[,, 1]
      scale     <- diag_var %*% as.matrix(nearPD(cov2cor(rWish))$mat) %*% diag_var
      for(j in seq_along(n)){
         for(k in seq_along(qmix)){
            ## Sample data 
            x <- rnvmix(n[j], qmix = qmix.[[k]], nu = nu, loc = loc, scale = scale)
            ## Call 'fitnvmix()' forcing estimatation of weights/densities
            ## using the function 'qmix.[[i]]'
            CPUused <- 
               system.time(fit <- fitnvmix(x, qmix = qmix.[[k]],
                                           mix.param.bounds = mix.param.bounds, 
                                           nu.init = nu.init, 
                                           size.subsample = size.subsample[j], 
                                           control = control, verbose = verbose))
            ## Re-perform using analytical weights?
            fit.analytical.nu <- if(analytical.weights){
               if(!special.mix[k] == "none"){
                  ## Analytical weights available
                  if(use.EM.mvt && special.mix[k] == "inverse.gamma"){
                     fit.mst(x, method = "BFGS")$df
                  } else {
                     fitnvmix(x, qmix = special.mix[[k]], # "inverse.gamma" or "pareto"
                              mix.param.bounds = mix.param.bounds, nu.init = nu.init, 
                              size.subsample = size.subsample[j], 
                              control = control, verbose = verbose)$nu
                  }
               }
            } else NA
            fitnvmix.results[i, j, k, ] <- 
               c(fit$nu.Ests[, 1], fit.analytical.nu, CPUused[1])
         }
         setTxtProgressBar(pb., (i-1)*length(n)+j) # update progress bar
      }
   }
   ## Save results
   duration <- Sys.time() - start
   attr(fitnvmix.results, 'duration') <- duration 
   attr(fitnvmix.results, 'nu')       <- nu # store true value of 'nu' as well
   attr(fitnvmix.results, 'nu.init')  <- nu.init
   ## Return (Note: 'fitnvmix_testing_plot()' returns input invisibly) 
   if(plot) fitnvmix_testing_plot(fitnvmix.results) else fitnvmix.results
}

#' Title: Plot results of numerical experiments for 'fitnvmix()'
#'
#' @param fitnvmix.results a 4-dimensional result array exactly as produced 
#'    by 'fitnvmix_testing()' 
#' @param index.qmix index of 'qmix' (third dimension 'fitnvmix.results') to be used 
#' @param index.d index of 'd' (first dimension of 'fitnvmix.results') to be used
#' @return plots estimate of 'nu' as a function of the number of iterations 
#'    for different sample sizes in a plot along with result for estimated weights
#'    (if available) for mixing variable in specified 'index.qmix' and dimension
#'    'index.d'.Additionally, 'fitnvmix.results' (the input) is invisbly returned. 
#' @author Erik Hintz
fitnvmix_testing_plot <- function(fitnvmix.results, index.qmix = 1,
                                  index.d = 1)
{
   ## Input checking
   stopifnot(length(dim(fitnvmix.results)) == 4)
   ## Grab sample sizes from dimension names 
   n         <- as.integer(dimnames(fitnvmix.results)$n)
   length.n  <- length(n)
   d         <- as.integer(dimnames(fitnvmix.results)$d)
   
   numiter <- length(dimnames(fitnvmix.results)[[4]]) - 2
   ## If 'nu.init' wasn't provided (=> estimated) start at iteration 0
   iter <- if(is.na(attr(fitnvmix.results, 'nu.init'))) 0:(numiter-1) else 1:numiter
   ## Colors
   pal  <- colorRampPalette(c("#000000", brewer.pal(8, name = "Dark2")[c(7, 3, 5)]))
   cols <- pal(length.n) # colors
   ## 'lwd's for different 'lty's 
   lwds.curr <- (if(length.n <= 6) lwds[1:length.n] else 
      c(lwds, rep(max(lwds), length.n - 5)))
   ## Set up legend
   lgn <- vector("expression", length.n)
   plot(NA, xlim = range(iter), 
        ylim = range(fitnvmix.results[index.d, , index.qmix, 1:(numiter+1)], na.rm = TRUE),
        xlab = "Number of iterations", ylab = expression(widehat(nu)))
   mtext(paste0("d = ", d[index.d]), side = 4)     
   for(j in seq_along(n)){
      ## Grab analytical result (possibly NA) and CPU needed
      analytical.result <- fitnvmix.results[index.d, j, index.qmix, "Analytical"]
      CPU               <- round(fitnvmix.results[index.d, j, index.qmix, "CPU"], 1)
      lines(iter, fitnvmix.results[index.d, j, index.qmix, 1:numiter],
            col = cols[j], lty = j, lwd = lwds.curr[j])
      lgn[[j]] <- bquote(hat(nu)[fitnvmix]~textstyle('for')~textstyle(n == .(n[j]))~
                            textstyle(' (')~textstyle(.(CPU))~textstyle('s)'))
      ## If available, put a mark at the end of the plot with analytical result
      if(!is.na(analytical.result)) points(tail(iter, 1), analytical.result, 
                                           col = cols[j], 
                                           pch = if(index.qmix == 1) 1 else 3)
   }
   legend('topright', legend = lgn, col = cols, lty = 1:length(n), box.lty = 0)
   ## Invisbly return data (eg in case it was called by 'fitnvix_testing()')
   invisible(fitnvmix.results)
}


## 2. Perform numerical experiments for 'pnvmix()' #############################

## 2.1 Estimate absolute error as a function of 'n' ############################
maxiter        <- 11 # maximum number of iterations
max.fun.evals  <- 3 * 2^8 * 2^seq(from = 1, to = maxiter, by = 1)
## Dimension(s) of the underlying normal variance mixture
d              <- c(5, 50, 100, 500) 
## Number of repetitions in each dimension
n              <- 15
## Specification of mixtures
qmix           <- c("inverse.gamma", "pareto")
## Parameter 'nu' of 'qmix' (degree-of-freedom for MVT(), alpha for PNVM())
nu             <- rep(2, 2)

if(doRUN){ # approximately 29 hours 
   set.seed(271)
   pnvmix.abserrors <-  pnvmix_testing_abserr(qmix, nu = nu, d = d, n = n, 
                                              max.fun.evals = max.fun.evals)
}

## 2.2 Estimate variance of the integrand with/without reordering
N        <- 5e4 # number of runs
n        <- 1e4 # number of points to estimate variance in each run
maxdim   <- 500 # dimension will be sampled uniformly from {3,...,maxdim}
mindim   <- 5
nu.lower <- 0.1 # nu is uniformly sampled between 'nu.lower', 'nu,upper'
nu.upper <- 5

if(doRUN){ # approximately 52 hours 
   set.seed(271)
   pnvmix.t.variances <- 
      precond_testing_variance(qmix = "inverse.gamma", N = N, n = n, 
                               mindim = mindim, maxdim = maxdim, 
                               nu.lower = nu.lower, nu.upper = nu.upper)
   
}

## 2.3 Estimate Sobol indices with/without reordering
## Parameters:
d  <- 10
nu <- 1 # df for multivariate t
n  <- 5e5 # Sample size of the 3 design matrices needed in 'sobolowen()'
## Specify seeds to generate 'scale' and 'upper'. The following 'seeds' were
## found by simulation to generate examples of high/moderate/no variance reduction
##    seed = 10575       => high variance reduction
##    seed = 96865       => slight variance reduction
##    seed = 30367       => no variance reduction / increase in variance 
## Those seeds are *specific* to dim = 10 and df = 1 and were found by 
## extensive simulation. 
seeds          <- c(10575, 96865, 30367) 
original.seed  <- 271 # seed to generate design matrices. 

if(doRUN){ # approximately 30 min
   pnvmix.t.sobolind <- pnvmix_estimate_sobolind(
      qmix = "inverse.gamma", d = d, nu = nu, n = n, seeds = seeds, 
      original.seed = original.seed)
   if(FALSE) warnings() # Warning: Conversion of the response to numeric
}

## 2.4 Compare CPU of 'pmvt()' from pkg 'mvtnorm' and 'pStudent()' from pkg 'nvmix'
d        <- c(seq(from = 5, to = 95, by = 10), seq(from = 105, to = 1455, by = 50))
tol      <- 1e-3
n        <- 15 # number of runs per dimension and method
rep      <- 3 # Repetitions for *each* setting (=> eliminate noise from time measurment)
df       <- 2 # degree-of-freedom parameter for the multivariate t dist'n 

if(doRUN){ # approximately 3.5 hours
   set.seed(271)
   pnvmix.t.timing <- pnvmix_timing_mvt(d, n = n, rep = rep, tol = tol, df = df)
}


## 3. Numerical experiments for 'dnvmix()' #####################################
## Specify simulation parameters for MVT() and PNVM()
qmix      <- c("inverse.gamma", "pareto")
nu.sample <- c(1, 0.5)
nu.dens   <- c(4, 2.5)
seed      <- 271

if(doRUN){ # approximately 3 min 
   dnvmix.results <- dnvmix_testing(qmix = qmix, nu.sample = nu.sample,
                                    nu.dens = nu.dens, seed = seed)
}


## 4. Numerical experiments for 'fitnvmix()' ###################################
## Parameter setttings
qmix <- c("inverse.gamma", "pareto")
n    <- c(250, 500, 1000, 2000, 5000) # sample sizes of the input data
d    <- c(10, 50) # dimensions
nu   <- 2.5 # true parameter to be sampled from 

if(doRUN){ # approximately 30 min
   set.seed(271) # for reproducibility
   ## Generate results
   fitnvmix.results <- fitnvmix_testing(qmix = qmix, n = n, d = d, nu = nu)
}


## 5. Data Analysis for 5 Stock portfolio ######################################

## Specifications for the mixing distributions considered 
qmix.strings <- c("constant", "inverse.gamma", "inverse.burr", "pareto")

if(doRUN){ # approx. 5 min 
   ### 5.1. Load data, grab stocks ################################################
   set.seed(123) # for reproducibility (RQMC methods depend on .Random.seed)
   data("SP500_const") # load the constituents data from qrmdata
   time <- c("2007-01-03", "2009-12-31") # time period
   x <- SP500_const[paste0(time, collapse = "/"),] # data
   stocks <- c("AAPL",  # Apple
               "ADBE",  # Adobe
               "INTC",  # Intel
               "ORCL",  # Oracle
               "GOOGL") # Google
   x <- x[, stocks]
   stopifnot(all(!is.na(x))) # no NAs 
   
   ### 5.2. Fitting marginal ARMA(1,1)-GARCH(1,1) models ##########################
   X <- -returns(x) # -log-returns
   n <- nrow(X) # 755
   d <- ncol(X) # 5
   uspec <- rep(list(ugarchspec(distribution.model = "std")), ncol(X)) # GARCH specs
   fit.ARMA.GARCH <- fit_ARMA_GARCH(X, ugarchspec.list = uspec) # fit ARMA-GARCH
   stopifnot(sapply(fit.ARMA.GARCH$error, is.null)) # NULL = no error
   fits <- fit.ARMA.GARCH$fit # fitted models
   resi <- lapply(fits, residuals, standardize = TRUE) # grab out standardized residuals
   Z <- as.matrix(do.call(merge, resi)) # standardized residuals
   stopifnot(is.matrix(Z), nrow(Z) == 755, ncol(Z) == 5) # fail-safe programming
   colnames(Z) <- colnames(x)
   stopifnot(nrow(Z) == n, ncol(Z) == d)
   
   ### 5.3. Fitting normal variance mixture models to the standardized residuals ##
   ## Note: qmix = "inverse.burr" not supported as no analytical density available
   qmix.functions <- list( 
      contant       = function(u, nu) rep(1, length(u)),
      inverse.gamma = function(u, nu) 1 / qgamma(1 - u, shape = nu/2, rate = nu/2),
      inverse.burr  = function(u, nu) (u^(-1/nu[2])-1)^(-1/nu[1]),
      pareto        = function(u, nu) (1-u)^(-1/nu))
   ## Specify parameter bounds for mixing parameters
   mix.bounds.1d <- c(0.5, 12) # parameter bounds for models with 1 param
   mix.bounds.2d <- 
      matrix(c(0.5, 0.5, 8, 8), ncol = 2) # parameter bounds for 'inverse.burr'
   ## Result objects 
   fit.analytical <- vector('list', 4)
   fit.estimated  <- vector('list', 4)
   qqplots <- vector('list', 4) 
   ## Fit using *analytical* densities/weights (not available for 'inverse.burr')
   for(j in 1:4){
      if(j == 3) next 
      t <- system.time(
         fit <- fitnvmix(Z, qmix = qmix.strings[j], mix.param.bounds = mix.bounds.1d, 
                         verbose = FALSE))
      attr(fit, "CPU") <- as.numeric(t[1]) # remove name
      fit.analytical[[j]] <- fit
   }
   ## Fit using *estimated* densities/weights and generate qq-plots
   for(j in 1:4){
      bounds <- if(j != 3) mix.bounds.1d else mix.bounds.2d
      t <- system.time(fit <- fitnvmix(Z, qmix = qmix.functions[[j]], 
                                       mix.param.bounds = bounds, verbose = FALSE))
      attr(fit, "CPU") <- as.numeric(t[1]) # remove name
      fit.estimated[[j]] <- fit
      t <- system.time(
         qq <- qqplot_maha(Z, qmix = qmix.functions[[j]], loc = fit$loc, 
                           scale = fit$scale, nu = fit$nu, verbose = FALSE))
      attr(qq, "CPU") <- as.numeric(t[1])
      qqplots[[j]] <- qq
   }
   
   ### 5.4. Estimating joint quantile shortfall probs for the different models,
   ### ie P(X_i <= F_{X_i}^{-1}(u), i = 1,...,d) ##################################
   n <- 32 
   u <- 1 - seq(0.95, to = 0.9995, length.out = n) # small levels 
   tailprobs <- array(NA, dim = c(length(u), length(qmix.strings), 2), 
                      dimnames = list(u = u, mix = qmix.strings, 
                                      c("univariate quantile", "exceedence prob")))
   for(j in 1:4){
      ## Obtain quantile of 'standard univariate nvmix' (mu = 0, sig = 1)
      tailprobs[, j, 1] <- qnvmix(u, qmix = qmix.functions[[j]],
                                  nu = fit.estimated[[j]]$nu)
      ## Can ignore 'loc' and restrict to correlation matrix of 'scale' 
      ## (as shortfall prob only depends on the underlying copula)
      tailprobs[, j, 2] <-
         pnvmix(matrix(rep(tailprobs[, j, 1], each = d),
                       ncol = d, byrow = TRUE), # matrix of upper limits
                qmix  = qmix.functions[[j]],
                scale = cov2cor(fit.estimated[[j]]$scale),
                nu = fit.estimated[[j]]$nu, 
                control = list(pnvmix.abstol = 1e-6)) # higher accuracy as prob's are small
   }
}

if(doSTORE){
   numerical_experiments_data <- list(fit.estimated = fit.estimated,
                                      fitnvmix.results = fitnvmix.results,
                                      fit.analytical = fit.analytical,
                                      qqplots = qqplots,
                                      pnvmix.t.variances = pnvmix.t.variances,
                                      pnvmix.t.sobolind = pnvmix.t.sobolind,
                                      pnvmix.t.timing = pnvmix.t.timing,
                                      tailprobs = tailprobs,
                                      dnvmix.results = dnvmix.results,
                                      pnvmix.abserrors = pnvmix.abserrors)
   ## 'version = 2' needed to avoid creating dependency on >= R 3.5.0
   ## (see https://github.com/r-lib/devtools/issues/1912#event-1942861670)
   save(numerical_experiments_data, file = "numerical_experiments_data.RData",
        version = 2)
}

## 6. Plot results #############################################################

## 6.1 Plot results for 'pnvmix()' #############################################

d    <- c(5, 50, 100, 500) 
qmix <- c("inverse.gamma", "pareto")
if(doPLOT){ # absolute errors as a fct of 'n'
   height <- 6 # for plotting
   width  <- 6
   ## Generate a plot in each dimension for each mixing variable separately
   for(i in seq_along(d)){
      for(j in seq_along(qmix)){
         if(doPDF){
            curr.qmix <- if(qmix[j] == "inverse.gamma") "t" else qmix[j] # keep fname short
            filename  <- paste0("fig_pnvmix_", curr.qmix, "_abserrors_dim", d[i],".pdf")
            pdf(file = filename, width = width, height = height)
         }
         pnvmix_testing_abserr_plot(pnvmix.abserrors, index.qmix = j, index.dim = i)
         if(!doPDF) mtext(paste0(qmix[j], "-mixture"), line = 1.5)
         if(doPDF) dev.off()
      }
   }
}

if(doPLOT){ # variances with/without reordering
   height <- 6 # for plotting
   width  <- 9
   ## Produce scatterplot
   if(doPDF) pdf(file="fig_pnvmix_t_variances_scatterplot.pdf", width = height, 
                 height = height)
   precond_testing_variance_plot(pnvmix.t.variances, scatterplot = TRUE)
   if(doPDF) dev.off()
   ## Produce histogram of variance ratios
   if(doPDF) pdf(file="fig_pnvmix_t_variances_histogram.pdf", width = height, 
                 height = height)
   precond_testing_variance_plot(pnvmix.t.variances, scatterplot = FALSE,
                                 density = TRUE)
   if(doPDF) dev.off()
}

if(doPLOT){ # sobol indices with/without re-ordering 
   height <- 10 # for plotting
   width  <- 7
   for(i in seq_along(seeds)){
      if(doPDF){
         fname <- paste("fig_pnvmix_t_sobolind_seed", seeds[i],".pdf", sep = "")
         pdf(file = fname, height = height, width = width)
      }
      pnvmix_estimate_sobolind_plot(pnvmix.t.sobolind, index.seed = i)
      if(doPDF) dev.off() 
   }
}

if(doPLOT){ # timing experiment
   height <- 6 # for plotting
   width  <- 6
   if(doPDF) pdf(file = "fig_pnvmix_t_timing.pdf", width = width, height = height)
   pnvmix_timing_mvt_plot(pnvmix.t.timing)
   if(doPDF) dev.off()
   if(doPDF) pdf(file = "fig_pnvmix_t_timing_ratio.pdf", width = width, height = height)
   pnvmix_timing_mvt_plot(pnvmix.t.timing, ratio = TRUE)
   if(doPDF) dev.off()
}

## 6.2 Plot results for 'dnvmix()' #############################################

if(doPLOT){
   height <- 6 # for plotting
   width  <- 6 
   for(i in seq_along(qmix)){
      if(doPDF){
         curr.qmix <- if(qmix[i] == "inverse.gamma") "t" else qmix[i] # keep fname short
         filename  <- paste0("fig_dnvmix_", curr.qmix, ".pdf")
         pdf(file = filename, width = width, height = height)
      }
      dnvmix_testing_plot(dnvmix.results, index.qmix = i)
      if(!doPDF) mtext(paste0(qmix[j], "-mixture"), line = 1.5)
      if(doPDF) dev.off()
   }
}

## 6.3 Plot results for 'fitnvmix()' ###########################################

d <- c(10, 50) # dimensions
if(doPLOT){
   height <- 6 # for pdf(..., height, width)
   width  <- 6
   names.qmix <- dimnames(fitnvmix.results)$qmix
   ## For each mixing quantile fct 'qmix' and each dimension a separate plot
   for(i in seq_along(qmix)){
      ## Get name of 'qmix' if possible 
      curr.qmix <-  if(names.qmix[i] == "inverse.gamma") "t" else 
         names.qmix[i] # keep fname short
      for(j in seq_along(d)){
         ## Generate .pdf?
         if(doPDF){
            fname <- paste("fig_fitnvmix_", curr.qmix, "_dim", d[j],".pdf", sep = "")
            pdf(file = fname, height = height, width = width)
         }
         fitnvmix_testing_plot(fitnvmix.results, index.qmix = i, index.d = j)
         if(!doPDF) mtext(paste0(qmix[i], "-mixture"), line = 1.5)
         if(doPDF) dev.off()
      }
   }
}

## 6.4 Plot results for the stock data analysis ################################

if(doPLOT){
   size <- 6 # width and height for doPDF()
   ## Result matrix for paper 
   res.matrix <- matrix(0, ncol = 4, nrow = 8)
   colnames(res.matrix) <- c("analytical", "", "estimated", "")
   rownames(res.matrix) <- rep(qmix.strings, each = 2)
   rownames(res.matrix)[c(2, 4, 6, 8)] <- "CPU"
   for(j in 2:4){ # nothing really estimated for "qmix = 'constant'"
      if(j == 3) next 
      ## Analytical results
      fit.nu.time <- round(c(fit.estimated[[j]]$nu,
                             attr(fit.estimated[[j]], 'CPU')), 2)
      res.matrix[(j-1)*2+(1:2), 1] <- fit.nu.time
      ## Estimated results
      fit.nu.time <- round(c(fit.analytical[[j]]$nu,
                             attr(fit.analytical[[j]], 'CPU')), 2)
      res.matrix[(j-1)*2+(1:2), 3] <- fit.nu.time
   }
   res.matrix[5, 3:4] <- round(fit.estimated[[3]]$nu, 2)
   res.matrix[6, 3] <- round(attr(fit.estimated[[3]], 'CPU'), 2)
   print(res.matrix)
   
   ## Produce QQ-Plots (results already obtained above)
   size <- 6 # width and height for doPDF()
   qmix.strings. <- 
      c("Constant", "Inverse-gamma", "Inverse-Burr", "Pareto") # capitalized 
   for(j in 1:4){
      if(doPDF) pdf(file = (file <- paste0("fig_qqplot_", qmix.strings[j], ".pdf")),
                    width = size, height = size)
      plot(qqplots[[j]]$q, qqplots[[j]]$maha2, 
           xlab = "Theoretical quantiles", ylab = "Sample quantiles", main = "",
           pch = 4)
      lines(qqplots[[j]]$q, qqplots[[j]]$q, lty = 2, lwd = 0.5) # diagonal
      mtext(qmix.strings.[j],  cex = 1.1, line = 1)
      if(doPDF) dev.off()
   }
   
   ## Plot shortfall probabilities 
   u    <- as.numeric(dimnames(tailprobs)$u)
   pal  <- colorRampPalette(c("#000000", brewer.pal(8, name = "Dark2")[c(7, 3, 5)]))
   cols <- pal(4) # colors
   
   if(doPDF) pdf(file = "fig_shortfallprob.pdf", width = size, height = size)
   plot(NA, xlim = range(u), ylim = range(tailprobs[,,2]), xlab = "u",
        ylab = expression(P(X[1]<=q[u],...,X[30]<=q[u])),
        log = "y")
   for(j in 1:4){
      lines(u, tailprobs[, j, 2], type = 'l', col = cols[j], lty = j, 
            lwd = lwds[j])
   }
   legend("bottomright", c("Pareto mixture", "Inverse-Burr mixture",
                           "Inverse-gamma mixture", "Multiv. normal"),
          col = rev(cols), lty = 4:1, lwd = lwds[4:1], box.lty = 0)
   if(doPDF) dev.off()
   ## Plot shortfall probabilities (standardized by normal)
   if(doPDF) pdf(file = "fig_shortfallprob_standardized.pdf", width = size, height = size)
   plot(NA, xlim = range(u), ylim = c(1, 1e2), xlab = "u",
        ylab = expression(P(X[1]<=q[u],...,X[30]<=q[u])~"standardized by the normal"),
        log = "y")
   for(j in 1:4){ # omit normal case 
      lines(u, tailprobs[, j, 2]/tailprobs[, 1, 2], lty = j,
            col = cols[j], lwd = lwds[j])
   }
   legend("topright", c("Pareto mixture", "Inverse-Burr mixture",
                        "Inverse-gamma mixture", "Multiv. normal"),
          col = rev(cols), lty = 4:1, lwd = lwds[4:1], box.lty = 0)
   if(doPDF) dev.off() 
}

##### END DEMO #######