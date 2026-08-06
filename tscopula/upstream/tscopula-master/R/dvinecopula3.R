#' D-vine copula processes of type 3
#'
#' Class of objects for d-vine copula processes. See \link{dvinecopula3} for more details.
#'
#' @slot name name of the d-vine copula process.
#' @slot modelspec list containing the family, rotation, and name of KPACF
#' @slot pars list comprising of the parameters.
#'
#' @export
#'
setClass("dvinecopula3", contains = "tscopula", slots = list(
  name = "character",
  modelspec = "list",
  pars = "list"
))

#' Constructor function for dvinecopula3 process
#'
#' This function sets up a stationary d-vine process of finite or infinite order based on a
#' sequence of Gaussian copulas with a finite number of non-Gaussian substitutions at specified lags.
#' The substituted families can be Gumbel, Clayton, Joe, Frank, t and BB1 copulas as implemented by the
#' \code{\link[rvinecopulib]{bicop_dist}} in the \code{rvinecopulib} package.
#'
#' For the substituted copulas (other than t and Frank) the user must specify the rotation that should be used for
#' positive dependencies (0 or 180) and the rotation that should be used for negative dependencies (90 or 270).
#'
#' The copulas are parameterized using the Kendall partial autocorrelation function (kpacf) specified
#' by the \code{kpacf} argument. The default choice is the kpacf of a standard ARMA process which is
#' implemented in the function \code{\link{kpacf_arma}}. The parameters
#' of the kpacf should be set as a list using the \code{pars} argument; the required parameters should usually
#' be clear from the documentation of the chosen kpacf function and must be correctly named.
#'
#' The \code{maxlag} parameter specifies the maximum lag of the process; a finite number gives a finite-order
#' stationary d-vine process, but the parameter may also be set to \code{Inf} for an infinite-order process.
#'
#' If one or more of the substituted copulas are t or BB1 copulas the argument \code{auxpar} should be used to
#' specify the additional parameters. These are the degree-of-freedom parameter for t and the delta parameter for BB1;
#' the former must be greater or equal 2 and the latter greater or equal 1.
#'
#' @param location vector of locations of non-Gaussian copula substitutions
#' @param family vector of family names for non-Gaussian copula substitutions
#' @param posrot vector of rotations for substituted families under positive dependence (default is 0)
#' @param negrot vector of rotations for substituted families under negative dependence (default is 90)
#' @param kpacf a character string giving the name of the Kendall pacf
#' @param pars a list containing the parameters of the model
#' @param auxpar vector of additional parameters for two-parameter copulas
#' @param maxlag a scalar specifying the maximum lag
#'
#' @return An object of class \linkS4class{dvinecopula3}.
#' @export
#'
#' @examples
#' dvinecopula3(location = c(1,4), family = c("Gumbel", "clayton"),
#' posrot = c(0, 180), negrot = c(90, 270), kpacf = "kpacf_arma",
#' pars = list(ar = 0.95, ma = 0.85), maxlag = 20)
dvinecopula3 <- function(location = 1,
                         family = "gumbel",
                         posrot = 0,
                         negrot = 90,
                         kpacf = "kpacf_arma",
                         pars = list(ar = 0.1, ma = 0.1),
                         auxpar = NA,
                         maxlag = Inf) {
  if (!(is(family, "character")))
    stop("copula family must be specified by name")
  if (is.null(names(pars)))
    stop("parameters should be named (p1 and p2 for exp/power)")
  fam <- tolower(family)
  nsub <- length(location)
  if (nsub > 1)
  {
    if (length(fam) == 1)
      fam <- rep(fam, nsub)
    if (length(posrot) == 1)
      posrot <- rep(posrot, nsub)
    if (length(negrot) == 1)
      negrot <- rep(negrot, nsub)
  }
  if ((length(posrot) != nsub) | (length(negrot) != nsub) | (length(fam) != nsub))
    stop("Length of family and rotations must
         be equal to length of location or 1")
  posrot[(fam == "frank") | (fam == "t")] <- 0 # radial symmetry
  negrot[(fam == "frank") | (fam == "t")] <- 0 # radial symmetry
  twoparfamily <- fam[(fam == "t") | (fam == "bb1")]
  ntwopar <- length(twoparfamily)
  if (ntwopar > 0){
    if (length(auxpar) != ntwopar)
      stop("Require second parameter for two-parameter copulas")
    pars$auxpar <- auxpar
  }
  modelspec <- list(location = location,
                    family = fam,
                    posrot = posrot,
                    negrot = negrot,
                    kpacf = kpacf,
                    maxlag = maxlag,
                    npar = length(unlist(pars)))
  output <- new("dvinecopula3",
      name = paste("type3-d-vine"),
      modelspec = modelspec,
      pars = pars
  )
  check <- mklist_dvine3(output, max(location), truncate = TRUE)
  output
}

#' @describeIn dvinecopula3 Coef Method for dvinecopula3 class
#'
#' @param object an object of the class.
#'
#' @export
#'
setMethod("coef", c(object = "dvinecopula3"), function(object) {
  if (length(object@pars) == 1) {
    return(object@pars[[1]])
  } else {
    return(unlist(object@pars))
  }
})

#' @describeIn dvinecopula3 Calculate Kendall's tau values for pair copulas in type 3 d-vine copula
#'
#' @param object an object of the class.
#' @param lagmax maximum value of lag.
#'
#' @export
#'
setMethod("kendall", c(object = "dvinecopula3"), function(object, lagmax = 20) {
  kpacf <- eval(parse(text = object@modelspec$kpacf))
  kpacf(lagmax, object@pars)
}
)


#' @describeIn dvinecopula3 Show method for dvinecopula3 class
#'
#' @param object an object of the class.
#'
#' @export
#'
setMethod("show", c(object = "dvinecopula3"), function(object) {
  cat("object class: ", is(object)[[1]], "\n", sep = "")
  cat("name: ", object@name, "\n", sep = "")
  cat("non-Gaussian substitutions:\n")
  locs <- object@modelspec$location
  cat(" - locations:", locs, "\n", sep = " ")
  cat(" - families:", object@modelspec$family, "\n", sep = " ")
  tau <- kendall(object, lagmax = max(locs))[locs]
  rot <- ifelse(tau >= 0, object@modelspec$posrot,object@modelspec$negrot)
  cat(" - rotations:", rot, "\n", sep = " ")
  cat(" - Kendall's tau:", round(tau,3), "\n", sep = " ")
  kpacf  <- object@modelspec$kpacf
  if (object@modelspec$maxlag != Inf)
    kpacf <- paste(kpacf, "with max lag", object@modelspec$maxlag)
  cat("KPACF: ", kpacf,"\n", sep = "")
  cat("parameters: \n")
  print(coef(object))
})

#' Make list of pair copulas for dvinecopula3 object
#'
#' @param x an object of class dvinecopula3
#' @param maxlag maximum possible lag to consider
#' @param truncate logical variable stating whether to truncate copulas with negligible dependence
#' @param tol tolerance for truncation
#'
#' @return a list of pair copulas
#' @keywords internal
#'
mklist_dvine3 <- function(x, maxlag, truncate, tol = 1){
  k <- min(maxlag, x@modelspec$maxlag)
  kpacf <- eval(parse(text = x@modelspec$kpacf))
  tauvals <- kpacf(k, x@pars)
  if (truncate)
    k <- max(c(1, which(abs(tauvals) > .Machine$double.eps^tol)))
  pc_list <- vector("list", k)
  auxpar <- 1
  for (i in 1:k) {
    fam <- "gauss"
    rot <- 0
    for (j in 1:length(x@modelspec$location)){
      if (i == x@modelspec$location[j]){
        fam <- tolower(x@modelspec$family)[j]
        if ((tauvals[i] >= 0) &
            (fam %in% c("gumbel", "clayton", "joe", "bb1")))
          rot <- x@modelspec$posrot[j]
        if ((tauvals[i] < 0) &
            (fam %in% c("gumbel", "clayton", "joe", "bb1")))
          rot <- x@modelspec$negrot[j]
      }
    }
    coppars <- ktau_to_par(
      family = fam,
      tau = tauvals[i]
    )
    if (fam == "t"){
      coppars <- c(coppars, x@pars$auxpar[auxpar])
      auxpar <- auxpar + 1
    }
    if (fam == "bb1"){
      par2 <- x@pars$auxpar[auxpar]
      par1 <- (coppars+2)/par2 -2
      coppars <- c(par1, par2)
      auxpar <- auxpar + 1
    }
    pc_list[[i]] <- rvinecopulib::bicop_dist(
      family = fam,
      rotation = rot,
      parameters = coppars)
  }
  pc_list
}

#' @describeIn dvinecopula3 Simulation method for dvinecopula3 class
#'
#' @param object an object of the class.
#' @param n length of realization.
#'
#' @export
#'
setMethod("sim", c(object = "dvinecopula3"), function(object, n = 1000) {
  pc_list <- mklist_dvine3(object, n-1, truncate = TRUE, tol = 1/3)
  simdvine(pc_list, n, innov = NA, start = NA)
})


#' Objective function for dvinecopula3 process
#'
#' @param theta parameters of kpacf
#' @param modelspec list specifying model
#' @param u data
#'
#' @return Value of objective function at parameters.
#'
#' @keywords internal
#'
dvinecopula3_objective <- function(theta, modelspec, u) {
  n <- length(u)
  kpacf <- eval(parse(text = modelspec$kpacf))
  tauvals <- kpacf((n-1), theta)
  if (is.na(sum(tauvals)))
    return(NA)
  k <- 1
  largetau <- (abs(tauvals) > .Machine$double.eps)
  if (sum(largetau) > 0)
    k <- max((1:(n-1))[largetau])
  k <- min(k, modelspec$maxlag)
  pc_list <- vector("list", k)
  tpar <- 1
  for (i in 1:k) {
    fam <- "gauss"
    rot <- 0
    for (j in 1:length(modelspec$location)){
      if (i == modelspec$location[j]){
        fam <- tolower(modelspec$family)[j]
        if ((tauvals[i] >= 0) &
            (fam %in% c("gumbel", "clayton", "joe", "bb1")))
          rot <- modelspec$posrot[j]
        if ((tauvals[i] < 0) &
            (fam %in% c("gumbel", "clayton", "joe", "bb1")))
          rot <- modelspec$negrot[j]
      }
    }
    coppars <- ktau_to_par(
      family = fam,
      tau = tauvals[i]
    )
    if (fam == "t"){
      auxpar <- theta[substring(names(theta),1,2) == "au"]
      coppars <- c(coppars, auxpar[tpar])
      tpar <- tpar + 1
    }
    if (fam == "bb1"){
      auxpar <- theta[substring(names(theta),1,2) == "au"]
      par2 <- auxpar[tpar]
      par1 <- (coppars+2)/par2 -2
      coppars <- c(par1, par2)
      tpar <- tpar + 1
    }
    pc_list[[i]] <- tryCatch(rvinecopulib::bicop_dist(
      family = fam,
      rotation = rot,
      parameters = coppars
    ),
    error = function(e) {
      return(NA)
    }
    )
    if (is.na(pc_list[[i]][[1]])) {
      return(NA)
    }
  }
  v <- cbind(u[1:(n - 1)], u[2:n])
  LL <- 0
  for (j in 1:k) {
    LL <- LL + sum(log(rvinecopulib::dbicop(u = v, family = pc_list[[j]])))
    if (j == k) {
      return(-LL)
    }
    n <- dim(v)[1]
    v <- cbind(
      rvinecopulib::hbicop(v[(1:(n - 1)), ], cond_var = 2, family = pc_list[[j]]),
      rvinecopulib::hbicop(v[(2:n), ], cond_var = 1, family = pc_list[[j]])
    )
  }
}

#' @describeIn dvinecopula3 Prediction method for dvinecopula2 class
#'
#' @param object an object of the class.
#' @param data vector of past data values.
#' @param x vector of arguments of prediction function.
#' @param type type of prediction function ("df" for density, "qf" for quantile function
#' or "dens" for density).
#'
#' @export
#'
setMethod("predict", c(object = "dvinecopula3"), function(object, data, x, type = "df") {
  pc_list <- mklist_dvine3(object, length(data)-1, truncate = TRUE, tol = 1/3)
  switch(type,
         "df" = Rblatt(pc_list, data, x),
         "qf" = IRblatt(pc_list, data, x),
         "dens" = Rblattdens(pc_list, data, x))

})

#' Residual function for dvinecopula3 object
#'
#' @param object a fitted dvinecopula3 object.
#' @param data the data to which copula is fitted.
#' @param trace extract trace instead of residuals.
#'
#' @return vector of model residuals
#' @keywords internal
#'
resid_dvinecopula3 <- function(object, data = NA, trace = FALSE){
  n <- length(data)
  pc_list <- mklist_dvine3(object, n-1, truncate = TRUE, tol = 1/3)
  k <- length(pc_list)
  for (i in 1:k)
    if (pc_list[[i]]$rotation %in% c(90,270))
      pc_list[[i]]$rotation <- 360 - pc_list[[i]]$rotation
  if (trace)
    target <- rep(0.5, n)
  else
    target <- data
  res <- rep(NA, n)
  res[1] <- target[1]
  if (k >1){
    for (i in 2:k){
      pcs <- lapply(1:(i-1), function(j) {
        replicate(i - j, pc_list[[j]], simplify = FALSE)
      })
      vc_short <- rvinecopulib::vinecop_dist(pcs, rvinecopulib::dvine_structure(i:1))
      vals <- c(data[1:(i-1)], target[i])
      res[i] <- rvinecopulib::rosenblatt(t(vals), vc_short)[i]
    }
  }
  pcs <- lapply(1:k, function(j) {
    replicate(k + 1 - j, pc_list[[j]], simplify = FALSE)
  })
  vc_short <- rvinecopulib::vinecop_dist(pcs, rvinecopulib::dvine_structure((k + 1):1))
  for (i in ((k+1):n)){
    vals <- c(data[(i-k):(i-1)], target[i])
    res[i] <- rvinecopulib::rosenblatt(t(vals), vc_short)[k+1]
  }
  qnorm(res)
}

#' Generalized lagging for fitted dvinecopula3 objects
#'
#' @param copula a dvinecopula3 object
#' @param data the data to which copula is fitted
#' @param lagmax the maximum lag value.
#' @param glagplot logical value indicating generalized lag plot.
#'
#' @return If \code{glagplot} is \code{TRUE} a list of generalized lagged datasets
#' of maximum length 9 is returned to facilitate a generalized lagplot.
#' If \code{glagplot} is \code{FALSE} a vector of length \code{lagmax} containing
#' the Kendall rank correlations for the generalized lagged datasets is returned.
#' @keywords internal
glag_for_dvinecopula3 <- function(copula, data, lagmax, glagplot = FALSE) {
  if (glagplot)
    lagmax <- min(lagmax, 9)
  pc_list <- mklist_dvine3(copula, lagmax, truncate = FALSE)
  k <- length(pc_list)
  n <- length(data)
  data <- cbind(as.numeric(data[1:(n - 1)]), as.numeric(data[2:n]))
  if (glagplot){
    output <- vector(mode = "list", length = k)
    output[[1]] <- data
  }
  else{
    output <- rep(NA, k)
    output[1] <- cor(data, method = "kendall")[1, 2]
  }
  if (k >1){
    for (i in 1:(k - 1)) {
      n <- dim(data)[1]
      data <-
        cbind(rvinecopulib::hbicop(data[(1:(n - 1)), ], pc_list[[i]], cond_var = 2),
              rvinecopulib::hbicop(data[(2:n), ], pc_list[[i]], cond_var = 1))
      if (glagplot)
        output[[i+1]] <- data
      else
        output[i+1] <- cor(data, method = "kendall")[1, 2]
    }
  }
  output
}
