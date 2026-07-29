### get_set_qqplot_param( ) #######################################################

##' @title  Retrieve specific default parameters for plotting and overwrite them
##' @return list with default (or changed) values used in 'plot.qqplot_maha()'
##' @author Erik Hintz and Marius Hofert
get_set_qqplot_param <- function(plot.pars = list(log = "")){
   ## Set up default pars:
   pars <- 
      list(log = "", xlim = NULL, ylim = NULL, xlab = "Theoretical quantiles", 
           ylab = "Sample quantiles", sub = "", main = "", plot_legend = TRUE, 
           plot_test = TRUE, plot_line = TRUE, pch = 1,  lty = 1:3, 
           col = c("black", "red", "azure4", "chocolate4"))
   if(length(plot.pars) > 0){
      ## If input provided, grab input controls and overwrite:
      names.default <- names(pars)
      pars[(names.provided <- names(plot.pars))] <- plot.pars
      if (length(unmatched <- names.provided[!names.provided %in% names.default]))
         warning("Unknown names in control: ", paste(unmatched, collapse = ", "))
   }
   ## *Some* checking
   stopifnot(pars$log %in% c("", "x", "y", "xy"),
             is.logical(pars$plot_legend),
             is.logical(pars$plot_test),
             is.logical(pars$plot_line),
             length(pars$pch) == 1)
   ## Return
   pars
}


### get_set_param() #######################################################

##' @title  Retrieve algorithm specific default parameters and overwrite them
##' @return list with default values for all functions in the 'nvmix' package
##' @author Erik Hintz and Marius Hofert
get_set_param <- function(control = list())
{
   ## Set up default controls:
   ctrl <- list(
      ## For pnvmix():
      mean.sqrt.mix = NULL,
      precond = TRUE,
      pnvmix.abstol = 1e-3,
      pnvmix.reltol = NA,
      cholesky.tol = 1e-9,
      pnvmix.do.ant = TRUE,
      maxiter.stored = 4, 
      ## For dnvmix():
      dnvmix.abstol = 1e-3,
      dnvmix.reltol = 1e-2, # If !NA, 'reltol' is used instead of 'abstol'
      dnvmix.max.iter.rqmc.pilot = 6,
      dnvmix.doAdapt = TRUE,
      dnvmix.tol.int.lower = 1e-100,
      dnvmix.order.lower = 5,
      dnvmix.tol.bisec = c(1e-6, 1e-1, 1e-1),
      dnvmix.tol.stratlength = 1e-50,
      dnvmix.max.iter.bisec = 15,
      ## For pgammamix():
      pgammamix.reltol = NA,
      pgammamix.abstol = 1e-3,
      ## For qnvmix():
      max.iter.newton = 45,
      newton.conv.abstol = 5e-4,
      newton.df.reltol = 2.5e-4,
      newton.logdens.abstol = 1e-2,
      newton.df.max.iter.rqmc = 350, # 'doubling' used here!
      qqplot.df.reltol = 5e-3,
      ## For fitnvmix():
      ## Algorithm specifications:
      ECMEstep = TRUE,
      ECMEstep.do.nu = TRUE,
      laststep.do.nu = FALSE,
      resample = FALSE,
      ## Tolerances:
      ECME.maxiter = 200,
      ECME.miniter = 5,
      max.iter.locscaleupdate = 50,
      weights.abstol = NA, # currently not used
      weights.reltol = 1e-2,
      weights.interpol.reltol = 1e-2,
      ECME.rel.conv.tol = c(1e-2, 1e-2, 1e-3), # [1] => 'loc'; [2] => 'scale'; [3] => 'nu'
      ## For the underlying 'optim':
      control.optim = list(maxit = 75),
      control.optim.laststep = list(),
      ## For riskmeasures:
      riskmeasures.abstol = NA,
      riskmeasures.reltol = 5e-2,
      ## For dependence measures:
      dependencemeasures.abstol = 1e-3,
      dependencemeasures.reltol = NA,
      ## For all (randomized) algorithms:
      method = "sobol",
      increment = "num.init", # "doubling" or "num.init"
      max.iter.rqmc = NA, # defined below, depending on 'increment'
      CI.factor = 3.5,
      fun.eval = c(2^7, 1e12),
      B = 15)
   if(length(control) > 0) {
      ## If input provided, grab input controls and overwrite:
      names.control <- names(ctrl)
      ctrl[(names.provided <- names(control))] <- control
      ## Did the user provide something that is not used?
      if (length(unmatched <- names.provided[!names.provided %in% names.control]))
         warning("Unknown names in control: ", paste(unmatched, collapse = ", "))
      ## If 'pnvmix.reltol' was provided, set 'pnvmix.abstol' to NA:
      if(any(names.provided == 'pnvmix.reltol') && !is.na(control$pnvmix.reltol))
         ctrl$pnvmix.abstol <- NA
      ## Check if 'method' and 'increment' were provided correctly
      ctrl$method     <- match.arg(ctrl$method,
                                   choices = c("sobol", "ghalton", "PRNG"))
      ctrl$increment  <- match.arg(ctrl$increment,
                                   choices = c("doubling", "num.init"))
      ## Now some more checkings: ('max.iter.rqmc' checked at the end)
      stopifnot(is.logical(ctrl$precond),
                ctrl$cholesky.tol > 0,
                is.logical(ctrl$pnvmix.do.ant), 
                ctrl$dnvmix.max.iter.rqmc.pilot >= 1,
                is.logical(ctrl$dnvmix.doAdapt),
                ctrl$dnvmix.tol.int.lower > 0,
                ctrl$dnvmix.order.lower > 0,
                ctrl$dnvmix.tol.bisec >0,
                ctrl$dnvmix.tol.stratlength > 0,
                ctrl$dnvmix.max.iter.bisec > 0,
                ctrl$max.iter.newton >= 0,
                ctrl$newton.conv.abstol >= 0,
                ctrl$newton.df.reltol >= 0,
                ctrl$newton.logdens.abstol >= 0,
                ctrl$newton.df.max.iter.rqmc > 0,
                ctrl$qqplot.df.reltol > 0,
                is.logical(ctrl$ECMEstep),
                is.logical(ctrl$ECMEstep.do.nu),
                is.logical(ctrl$laststep.do.nu),
                is.logical(ctrl$resample),
                ctrl$ECME.maxiter > 0,
                ctrl$ECME.miniter >= 0,
                ctrl$ECME.maxiter >= ctrl$ECME.miniter,
                ctrl$max.iter.locscaleupdate > 0,
                ctrl$weights.abstol > 0 | ctrl$weights.reltol > 0, 
                ctrl$weights.interpol.reltol > 0,
                ctrl$riskmeasures.reltol > 0 | ctrl$riskmeasures.abstol > 0,
                length(ctrl$ECME.rel.conv.tol) == 3, ctrl$ECME.rel.conv.tol >= 0,
                ctrl$CI.factor >= 0,
                length(ctrl$fun.eval) == 2, ctrl$fun.eval >= 0,
                ctrl$B > 1) # If B = 1 error estimates are NA => need B > 1
   }
   ## Define 'max.iter.rqmc': If it was not provided (=> NA), set defaults
   if(is.na(ctrl$max.iter.rqmc)) {
      ctrl$max.iter.rqmc <- if(ctrl$increment == "doubling") 15 else 1250
   } else {
      ## If it was provided (=> not NA), check if it's reasonable
      stopifnot(ctrl$max.iter.rqmc > 1)
   }
   ## Return
   ctrl
}


##' @title Define quantile function from use provided input
##' @param qmix see calling function (string, list or function)
##' @param rmix see calling function (string, list or function)
##' @param callingfun string, name of calling function (eg "pnvmix")
##' @param groupings d-vector giving grouping in case of generalized NVM dist'n 
##' @return list with three elements. First element is a 'function<args> 'where 
##'         <args> is
##             * (u) for qmix + callingfun =
##                - pnvmix/pgammamix 
##                - rnvmix 
##                - dnvmix 
##                - qnvmix/qgammamix
##                - ESnvmix
##             * (u, nu) for qmix + callingfun =
##                - fitnvmix
##             * (n) for rmix + callingfun =
##                - pnvmix/pgammamix
##                - rnvmix 
##          in the other cases (rmix +  callingfun in ("dnvmix", "fitnvmix", 
##          "qnvmix", "qgammamix")), get_mix_ throws an error (for adaptive 
##          procedures 'qmix' must be provided). 
##          The 2nd element is 'special.mix', a character string  (one of 
##          "constant", "inverse.gamma", "pareto") for special mixing
##          distributions.
##          The third element is 'mean.sqrt.mix' (numeric or NULL), the
##          fourth element 'use.q' (logical if returned function is a quantile fun)
##          The fourth element 'param' is either special mixing parameter ('df' or
##          'alpha' for mix = "inverse.gamma" or "pareto") or NA
##' @author Erik Hintz
get_mix_ <- function(qmix = NULL, rmix = NULL, callingfun = NULL, 
                     groupings = NULL, ...)
{
   ## Grab ellipsis argument
   ell.args   <- list(...)
   n.ell.args <- names(ell.args)
   ## Determine if generalized or non-generalized normal variance mixture 
   is.grpd <- if(!is.null(groupings)){
      d <- length(groupings <- as.vector(groupings))
      stopifnot(all(groupings %in% 1:d)) 
      num.groups <- length(unique(groupings)) # number of different groups 
      (num.groups > 1) 
   } else {
      ## Groupings is null => non-generalized normal variance mixture
      num.groups <- 1
      FALSE
   }
   ## Determine if 'qmix' or 'rmix' used; 'qmix' always has priority over 'rmix'
   use.q <- !is.null(qmix) # logical if quantile function to be used
   ## If 'qmix' not provided, check if 'rmix' provided
   if(!use.q & is.null(rmix)) stop("Either 'qmix' or 'rmix' must be provided")
   mix_usr <- if(use.q){ # relevant 'qmix' or 'rmix'
      mix.prov <- "qmix" # for potential error/warning messages further below
      qmix
   } else {
      mix.prov <- "rmix" # for potential error/warning messages further below
      rmix 
   }
   ## In the grouped case, "qmix" must be provided (=> o/w cannot draw comonotone mixing rvs)
   if(is.grpd & !use.q)
      stop("In the grouped case, argument 'qmix' must be provided")
   ## Calling functions for which 'qmix' *must* be provided (due to adaptiveness)
   need.qmix <- c("dnvmix", "dgammamix", "fitnvmix", "qnvmix", "qgammamix", "ESnvmix",
                  "dgnvmix", "pgnvmix", "dependencemeasures")
   if(any(callingfun == need.qmix) & !use.q) 
      stop(paste(callingfun, "()", " needs argument 'qmix' to be provided", sep = ""))
   special.mix   <- NA # special mixing variable? TBD below
   param         <- NA # special parameter? TBD below
   mean.sqrt.mix <- NULL # E(sqrt(W)) for 'pnvmix()'
   ## 'mix_' is one of {function(u, nu), function(u), function(n)} depending
   ## on provided 'qmix'/'rmix' and 'callingfun'
   mix_ <-  if(callingfun == "fitnvmix"){ # 'mix_usr' is definitely a 'qmix' 
      if(is.character(mix_usr)){
         mix_usr <- match.arg(mix_usr, choices = c("constant", "inverse.gamma", "pareto"))
         special.mix <- mix_usr # for later
         switch(mix_usr,
                "constant" = {
                   function(u, nu) rep(1, length(u))},
                "inverse.gamma" = {
                   function(u, nu) 1 / qgamma(1 - u, shape = nu/2, rate = nu/2)},
                "pareto" = {
                   function(u, nu) (1-u)^(-1/nu)
                },
                stop(paste0("Currently unsupported '", mix.prov,"'")))
      } else if(is.function(mix_usr)){  # 'mix_usr is the quantile function of F_W
         function(u, nu) mix_usr(u, nu)
      } else {
         stop("fitnvmix() needs 'qmix' provided as a string or function")
      }
   } else if(!is.grpd) { # case of a normal variance mixture (not generalized)
      if(is.character(mix_usr)){
         ## 'qmix' is a character vector specifying supported mixture distributions (utilizing '...')
         mix_usr <- match.arg(mix_usr, 
                              choices = c("constant", "inverse.gamma", "pareto"))
         special.mix <- mix_usr # for later
         switch(mix_usr,
                "constant" = {
                   mean.sqrt.mix <- 1
                   function(u) rep(1, length(u)) # also works for 'rmix' 
                },
                "inverse.gamma" = {
                   if(hasArg(df)) {
                      df <- list(...)$df
                   } else if(hasArg(nu)) {
                      nu <- list(...)$nu
                      df <- nu
                   } else { 
                      stop(paste(mix.prov, " = \"inverse.gamma\"' requires 'df' to be provided.", sep = ""))
                   }
                   ## Still allow df = Inf (normal distribution)
                   stopifnot(is.numeric(df), length(df) == 1, df > 0)
                   if(is.finite(df)){
                      param <- df
                      df2 <- df / 2
                      mean.sqrt.mix <- sqrt(df) * gamma(df2) / (sqrt(2) * gamma((df+1)/2)) 
                      if(use.q) function(u) 1 / qgamma(1 - u, shape = df2, rate = df2) else
                         function(n) 1 / rgamma(n, shape = df2, rate = df2)
                   } else {
                      special.mix <- "constant"
                      mean.sqrt.mix <- 1
                      function(u) rep(1, length(u)) # also works for 'rmix' 
                   }
                },
                "pareto" = {
                   if(hasArg(alpha)) {
                      alpha <- list(...)$alpha
                   } else if(hasArg(nu)){
                      nu <- list(...)$nu
                      alpha <- nu
                   } else { 
                      stop(paste(mix.prov, " = \"pareto\"' requires 'alpha' to be provided.", sep = ""))
                   }
                   param <- alpha
                   mean.sqrt.mix <- if(alpha > 0.5) alpha/(alpha-0.5) else NULL
                   if(use.q) function(u) (1-u)^(-1/alpha) else 
                      function(n) (1 - runif(n))^(-1/alpha)
                },
                stop(paste0("Currently unsupported '", mix.prov,"'")))
      } else if(is.list(mix_usr)) {
         ## 'mix_usr' is a list of the form (<character string>, <parameters>)
         stopifnot(length(mix_usr) >= 1, is.character(distr <- mix_usr[[1]]))
         mix_usr_ <- if(use.q) paste0("q", distr) else paste0("r", distr)
         if(!existsFunction(mix_usr_))
            stop("No function named '", mix_usr_, "'.")
         function(u) 
            do.call(mix_usr_, append(list(u), mix_usr[-1]))
      } else if(is.function(mix_usr)){
         ## 'mix_usr' is the quantile function or a RNG for W 
         if(use.q) function(u) mix_usr(u, ...) else function(n) mix_usr(n, ...)
      } else stop(paste(mix.prov, "must be a character string, list or quantile function."))
   } else if(is.grpd){
      ## Case 1: 'mix_usr' is one string => all mixing rv's of same type
      if(is.character(mix_usr)){
         mix_usr <- match.arg(mix_usr, 
                              choices = c("inverse.gamma", "pareto"))
         special.mix <- mix_usr # for later
         switch(mix_usr, # note: no "constant" here! 
                "inverse.gamma" = {
                   if(hasArg(df)) {
                      df <- ell.args$df
                   } else if(hasArg(nu)) {
                      nu <- ell.args$nu
                      df <- nu
                   } else { 
                      stop(paste(mix.prov, " = \"inverse.gamma\"' requires 'df' to be provided.", sep = ""))
                   }
                   param <- df
                   ## Still allow df = Inf (normal distribution)
                   stopifnot(is.numeric(df), length(df) == num.groups, all(df > 0))
                   fin.df <- is.finite(df)
                   ## Construct 'mean.sqrt.mix_' (length num.groups)
                   mean.sqrt.mix_ <- rep(1, num.groups)
                   ## If X ~ InvGam(nu/2, nu/2) and df > 1
                   ## => E(sqrt(X)) = sqrt(df/2) * gamma((df-1)/2) / gamma(df/2);
                   ## If df <= 1 => E(sqrt(X)) DNE => take some proxy 
                   dffingr <- pmax(df[fin.df], 1.0001)
                   mean.sqrt.mix_[fin.df] <- sqrt(dffingr/2) * 
                      gamma((dffingr-1)/2) / gamma(dffingr/2)
                   ## Construct 'mean.sqrt.mix' (length d, i'th element = E(sqrt(W_i))
                   mean.sqrt.mix <- mean.sqrt.mix_[groupings]
                   ## Quantile fun or rng as function of (u, df) or (n, df)
                   mixfun <- if(use.q){
                      function(u, df){
                         if(is.finite(df)){
                            1/qgamma(1 - u, shape = df/2, rate = df/2)
                         } else rep(1, length(u))} 
                   } else {
                      function(n, df){
                         if(is.finite(df)){
                            1/rgamma(n, shape = df/2, rate = df/2)
                         } else rep(1, n)} 
                   }
                   ## Construct grouped quantile fun/RNG
                   function(u){
                      res <- sapply(1:num.groups, function(i) 
                         mixfun(u, df = df[i])) # also works with 'rmix'
                      if(is.matrix(res)) res else rbind(res, deparse.level = 0)
                   } 
                },
                "pareto" = {
                   if(hasArg(alpha)) {
                      alpha <- ell.args$alpha
                   } else if(hasArg(nu)){
                      nu <- ell.args$nu
                      alpha <- nu
                   } else { 
                      stop(paste(mix.prov, " = \"pareto\"' requires 'alpha' to be provided.", sep = ""))
                   }
                   param <- alpha
                   stopifnot(is.numeric(alpha), length(alpha) == num.groups, 
                             all(alpha > 0))
                   alpha.gr <- (alpha > 0.5)
                   ## Construct 'mean.sqrt.mix_' (length num.groups)
                   mean.sqrt.mix_ <- rep(1, num.groups)
                   mean.sqrt.mix_[alpha.gr] <- alpha[alpha.gr]/(alpha[alpha.gr]-0.5)
                   ## Construct 'mean.sqrt.mix' (length d, i'th element = E(sqrt(W_i))
                   mean.sqrt.mix <- mean.sqrt.mix_[groupings]
                   ## Quantile fun or rng as function of (u, df) or (n, df)
                   mixfun <- if(use.q){
                      function(u, alpha) (1-u)^(-1/alpha)
                   } else {
                      function(n, alpha) (1 - runif(n))^(-1/alpha)
                   }
                   ## Construct grouped quantile fun/RNG
                   function(u){
                      res <- sapply(1:num.groups, function(i) 
                         mixfun(u, df = alpha[i])) # also works with 'rmix'
                      if(is.matrix(res)) res else rbind(res, deparse.level = 0)
                   } 
                },
                stop(paste0("Currently unsupported '", mix.prov,"'")))
      } else {
         ## 'mix_usr' must be a list with 'num.groups' elements
         stopifnot(is.list(mix_usr), length(mix_usr) == num.groups)
         ## Element of 'mix_usr' can be functions or lists;
         addArgs   <- vector("list", num.groups) # for handling ellipsis (...) etc
         hasaddArg <- logical(num.groups)
         ## Go through the list and check arguments
         for(i in 1:num.groups){
            if(is.list(mix_usr[[i]])){ # element i of 'qmix' is a list
               stopifnot(length(mix_usr[[i]]) >= 1, 
                         is.character(distr <- mix_usr[[i]][[1]]))
               ## Name of the function to be called 
               mix_usr_i <- if(use.q) paste0("q", distr) else paste0("r", distr)
               if(!existsFunction(mix_usr_i))
                  stop("No function named '", mix_usr_i, "'.")
               hasaddArg[i] <- if(length(mix_usr[[i]][-1]) > 0){
                  addArgs[i] <- list(mix_usr[[i]][-1]) # => additional arguments
                  TRUE
               } else FALSE # => no additional arguments
               mix_usr[i] <- mix_usr_i # called below via 'do.call(..)' 
            } else if (is.function(mix_usr[[i]])) { # element i of 'qmix' is a fun
               l.addargs_i <- length(addargs_i <- names(
                  formals(mix_usr[[i]]))[-1]) # first one is 'u'
               ## Match 'addargs_i' with arguments provided via (...) 
               if(l.addargs_i > 0){
                  whichmatch_i <- sapply(seq_along(addargs_i), function(ii) 
                     which(addargs_i[ii] == n.ell.args)[1])
                  ## All arguments provided?
                  if(length(whichmatch_i) != length(addargs_i))
                     stop(paste("Function specified in element", i, "of", mix.prov, "is missing at least one argument."))
                  hasaddArg[i] <- TRUE
                  addArgs[i] <- list(ell.args[whichmatch_i]) # store additional args 
               } else hasaddArg[i] <- FALSE
            } else {
               stop(mix.prov, " must be either a character string, a list of functions or a list of lists.")
            }
         } # for()
         ## Build final return function 
         function(u){
            ## Call all quantile functions
            W. <- sapply(1:num.groups, function(i){ 
               if(hasaddArg[i]) do.call(mix_usr[[i]], append(list(u), addArgs[[i]])) else
                  do.call(mix_usr[[i]], list(u))})
            if(is.matrix(W.)) W. else rbind(W., deparse.level = 0)
         }
      }
   } 
   ## Return
   list(mix_ = mix_, special.mix = special.mix, mean.sqrt.mix = mean.sqrt.mix,
        use.q = use.q, param = param)
}
