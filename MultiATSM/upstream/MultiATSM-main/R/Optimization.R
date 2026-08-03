#' Perform the minimization of ML function
#'
#' @param ML_fun vector-valued objective ML function
#' @param ListInputSet list containing :
#'    \enumerate{
#'        \item a starting value for K1XQ and/or SSZ
#'        \item a variable label among the following:
#'                  \itemize{
#'                      \item 'Jordan' or 'Jordan; stationary' for single countries setups or
#'                      'Jordan MultiCountry' or 'Jordan MultiCountry; stationary'.  All cases related to the computation of a K1XQ parameter;
#'                      \item 'psd': PSD matrix, used for JPS-based models. It relates to the SSZ parameter;
#'                      \item  'BlockDiag': block diagonal matrix, used for JPS-based models. It relates to the SSZ parameter.
#'                      \item 'JLLstructure': to impose the zero-restrictions on the SSZ term along
#'                              the lines of the JLL models
#'                          }
#'   }
#' @param FactorLabels list. Labels for all variables present in the model, as returned by \code{\link{LabFac}}.
#' @param Economies character vector. Names of the \code{C} economies included in the system.
#' @param ModelType character. Model type to be estimated. Permissible choices: "JPS original", "JPS global", "GVAR single", "JPS multi", "GVAR multi", "JLL original", "JLL No DomUnit", "JLL joint Sigma".
#' @param JLLinputs List. Inputs for JLL model estimation (see \code{\link{JLL}}). Default is NULL.
#' @param GVARinputs List. Inputs for GVAR model estimation (see \code{\link{GVAR}}). Default is NULL.
#' @param tol convergence tolerance (scalar). Default value is 1e-4.
#' @param EstType Available options are"BFGS" and/or "Nelder-Mead".
#' @param TimeCount computes the required time for estimation of the model. Default is TRUE.
#' @param verbose Logical flag controlling function messaging.
#'
#' @keywords internal

Optimization_PE <- function(ML_fun, ListInputSet, FactorLabels, Economies, ModelType, JLLinputs = NULL,
                            GVARinputs = NULL, tol = 1e-4, EstType, TimeCount = TRUE, verbose) {
  # 1) Vectorize initial guesses from K1XQ and SSZ
  x0 <- c()
  for (i in seq_len(length(ListInputSet))) {
    Temp_x <- GetAuxPara(ListInputSet[[i]]$Value, ListInputSet[[i]]$Label, Economies, FactorLabels, JLLinputs)
    x0 <- rbind(x0, t(t(as.vector(Temp_x))))
  }

  # 2) Likelihood function:
  MLvec <- function(...) {
    FunctionML_vec(..., ML_fun, ListInputSet, ModelType, FactorLabels,
      Economies, JLLinputs, GVARinputs,
      WithEstimation = TRUE
    )
  }

  # 3) Optimization of the llk
  if (TimeCount) {
    start_time <- Sys.time()
  }
  x_opt <- OptimizationSetup_ATSM(x0, MLvec, EstType, tol, verbose)
  if (TimeCount) {
    Optimization_Time(start_time, verbose)
  }

  # 4) Build the full auxiliary vector, including concentrated parameters
  Up_Temp <- Update_ParaList(x_opt, ModelType, FactorLabels, Economies, JLLinputs, GVARinputs, ListInputSet)

  ML_opt <- function(...) {
    FunctionML_vec(..., ML_fun,
      ListInputSet = Up_Temp, ModelType, FactorLabels,
      Economies, JLLinputs, GVARinputs, WithEstimation = FALSE
    )
  }
  out <- ML_opt(x = x_opt)

  return(out)
}

#############################################################################################################
#' Perform the optimization of the log-likelihood function of the chosen ATSM
#'
#' @param MLEinputs list. Contains the inputs for building the log-likelihood function (see \code{\link{InputsForOpt}}).
#' @param StatQ A logical value indicating whether to impose that the largest eigenvalue under Q is strictly
#'                          smaller than 1. Set TRUE to impose this restriction.
#' @param DataFreq  character vector specifying the data frequency. Available options: "Daily All Days", "Daily Business Days",
#'                      "Weekly", "Monthly", "Quarterly", "Annually".
#' @param FactorLabels list. Labels for all variables present in the model, as returned by \code{\link{LabFac}}.
#' @param Economies character vector. Names of the \code{C} economies included in the system.
#' @param ModelType character. Model type to be estimated. Permissible choices: "JPS original", "JPS global", "GVAR single", "JPS multi", "GVAR multi", "JLL original", "JLL No DomUnit", "JLL joint Sigma".
#' @param tol numeric. Convergence tolerance. The default is 1e-4.
#' @param EstType Available options are"BFGS" and/or "Nelder-Mead".
#' @param TimeCount Logical. If TRUE, computes the time required for model estimation. Default is TRUE.
#' @param BS_outputs Logical. If TRUE, generates a simplified output list in the bootstrap setting. Default is FALSE.
#' @param verbose Logical flag controlling function messaging. Default is TRUE.
#'
#' @examples
#' LoadData("CM_2024")
#' ModelType <- "JPS original"
#' Economy <- "Brazil"
#' t0 <- "01-05-2007" # Initial Sample Date (Format: "dd-mm-yyyy")
#' tF <- "01-12-2018" # Final Sample Date (Format: "dd-mm-yyyy")
#' N <- 1
#' GlobalVar <- "Gl_Eco_Act" # Global Variables
#' DomVar <- "Eco_Act" # Domestic Variables
#' DataFreq <- "Monthly"
#' StatQ <- FALSE
#'
#' FacLab <- LabFac(N, DomVar, GlobalVar, Economy, ModelType)
#' ATSMInputs <- InputsForOpt(t0, tF, ModelType, Yields, GlobalMacro, DomMacro,
#'   FacLab, Economy, DataFreq,
#'   CheckInputs = FALSE, verbose = FALSE
#' )
#'
#' OptPara <- Optimization(ATSMInputs, StatQ, DataFreq, FacLab, Economy, ModelType, verbose = FALSE)
#'
#' @return
#' An object of class 'ATSMModelOutputs' containing model outputs after the optimization of the chosen ATSM specification.
#'
#' @section Available Methods:
#' - `summary(object)`
#'
#' @references
#'
#' \itemize{
#'  \item Candelon, C. and Moura, R. (2024). “A Multicountry Model of the Term Structures of Interest Rates with a GVAR.”
#'  Journal of Financial Econometrics 22 (5): 1558–87.
#'  \item Jotikasthira, C; Le, A. and Lundblad, C (2015). “Why Do Term Structures in Different Currencies Co-Move?”
#'  Journal of Financial Economics 115: 58–83.
#'  \item Joslin, S,; Priebsch, M. and Singleton, K. (2014). “Risk Premiums in Dynamic Term Structure Models with Unspanned Macro Risks.”
#'  Journal of Finance 69 (3): 1197–1233.
#'  \item Joslin, S., Singleton, K. and Zhu, H. (2011). "A new perspective on Gaussian dynamic term structure models".
#'  The Review of Financial Studies.
#'  \item Le, A. and Singleton, K. (2018). "A Small Package of Matlab Routines for the Estimation of Some Term Structure Models."
#'  Euro Area Business Cycle Network Training School - Term Structure Modelling.
#' }
#' @export

Optimization <- function(MLEinputs, StatQ, DataFreq, FactorLabels, Economies, ModelType, tol = 1e-4,
                         EstType = c("BFGS", "Nelder-Mead"), TimeCount = TRUE, BS_outputs = FALSE, verbose = TRUE) {
  if (verbose) {
    message("2) ATSM ESTIMATION : POINT ESTIMATE ANALYSIS")
    message("2.1) Estimating ATSM ...")
  }
  GVARinputs <- MLEinputs$GVARinputs
  JLLinputs <- MLEinputs$JLLinputs

  # Prepare optimization
  ModelParaList <- list()

  if (any(ModelType == c("JPS original", "JPS global", "GVAR single"))) {
    for (i in 1:length(Economies)) {
      MLEinputsCS <- MLEinputs[[Economies[i]]]
      Economy <- Economies[i]
      # 1) Build the objective function
      ML_fun <- MLFunction(MLEinputsCS, Economy, DataFreq, FactorLabels, ModelType, BS_outputs)
      # 2) Set the optimization settings
      VarLab <- ParaLabelsOpt(ModelType, StatQ, MLEinputsCS, BS_outputs)

      if (verbose) message(paste(" ... for country:", Economies[i], ". This may take several minutes."))
      ModelParaList[[ModelType]][[Economies[i]]] <- Optimization_PE(ML_fun, VarLab, FactorLabels, Economies, ModelType,
        JLLinputs, GVARinputs, tol, EstType,
        TimeCount = TimeCount, verbose
      )
    }
  } else {
    # 1) Build the objective function
    ML_fun <- MLFunction(MLEinputs, Economies, DataFreq, FactorLabels, ModelType, BS_outputs)
    # 2) Set the optimization settings
    VarLab <- ParaLabelsOpt(ModelType, StatQ, MLEinputs, BS_outputs)

    if (verbose) message("... This may take several minutes.")
    ModelParaList[[ModelType]] <- Optimization_PE(ML_fun, VarLab, FactorLabels, Economies, ModelType, JLLinputs,
      GVARinputs, tol, EstType,
      TimeCount = TimeCount, verbose
    )
  }

  # Store metadata inside the class without explicitly exporting it
  attr(ModelParaList, "ModelOutInfo") <- list(
    Outs = ModelParaList[[ModelType]], Economies = Economies,
    ModelType = ModelType
  )

  return(structure(ModelParaList, class = "ATSMModelOutputs"))
}


##############################################################################################
#' Create the variable labels used in the estimation
#'
#' @param ModelType a string-vector containing the label of the model to be estimated
#' @param WishStationarityQ User must set TRUE is she wishes to impose the largest eigenvalue under the Q to be strictly
#'                       smaller than 1. Otherwise set FALSE
#' @param MLEinputs Set of inputs that are necessary to the log-likelihood function
#' @param BS_outputs Generates simplified output list in the bootstrap setting. Default is set to FALSE.
#'
#' @keywords internal

ParaLabelsOpt <- function(ModelType, WishStationarityQ, MLEinputs, BS_outputs = FALSE) {
  ParaLabelsList <- list()

  # a) K1XQ
  K1XQType <- K1XQStationary(WishStationarityQ)$SepQ

  if (ModelType %in% c("JPS original", "JPS global", "GVAR single")) {
    ParaLabelsList[[ModelType]]$K1XQ <- K1XQStationary(WishStationarityQ)$SepQ
  } else {
    ParaLabelsList[[ModelType]]$K1XQ <- K1XQStationary(WishStationarityQ)$JointQ
  }

  # b) SSZ
  if (ModelType %in% c("JPS original", "JPS global", "JPS multi")) {
    ParaLabelsList[[ModelType]]$SSZ <- "SSZ: psd"
  } else if (ModelType %in% c("GVAR single", "GVAR multi")) {
    ParaLabelsList[[ModelType]]$SSZ <- "SSZ: BlockDiag"
  } else if (ModelType == "JLL joint Sigma") {
    ParaLabelsList[[ModelType]]$SSZ <- "SSZ: JLLstructure"
  }
  # Ensures that the structure of the Variance-covariance matrix of the JLL is preserved

  # 3.2) Initial guesses for Variables that will be concentrated out of from the log-likelihood function
  K1XQ <- MLEinputs$K1XQ
  if (ModelType %in% c("JLL original", "JLL No DomUnit")) {
    SSZ <- NULL
  } else {
    SSZ <- MLEinputs$SSZ
  }

  # Prepare list for optimization
  VarArgList <- list()
  VarArgList$K1XQ <- list(K1XQ, ParaLabelsList[[ModelType]]$K1XQ)
  if (!(ModelType %in% c("JLL original", "JLL No DomUnit"))) {
    VarArgList$SSZ <- list(SSZ, ParaLabelsList[[ModelType]]$SSZ)
  }

  LabelVar <- c("Value", "Label") # Elements of each parameter
  for (d in 1:(length(VarArgList))) {
    names(VarArgList[[d]]) <- LabelVar
  }

  return(VarArgList)
}

################################################################################################################
#' Impose stationarity under the Q-measure
#'
#' @param StationaryEigenvalues Binary variable: set "1" if the user wishes the largest eigenvalue
#'                            to be strictly smaller than 1. Set "0", otherwise
#'
#' @keywords internal

K1XQStationary <- function(StationaryEigenvalues) {
  K1Type <- list()

  if (StationaryEigenvalues == 1) {
    K1Type$SepQ <- paste("K1XQ: ", "Jordan", "; stationary", sep = "")
    K1Type$JointQ <- paste("K1XQ: ", "Jordan", " MultiCountry", "; stationary", sep = "")
  } else {
    K1Type$SepQ <- paste("K1XQ: ", "Jordan", sep = "")
    K1Type$JointQ <- paste("K1XQ: ", "Jordan", " MultiCountry", sep = "")
  }

  return(K1Type)
}

#########################################################################################################
#' Optimization routine for the entire selected ATSM
#'
#' @param x0 List containing features for estimation of the risk-neutral parameters.
#' @param MLvec Log-likelihood function
#' @param EstType Estimation type. Available options: "BFGS" and "Nelder-Mead".
#' @param tol convergence tolerance (scalar). Default value is set as 1e-4.
#' @param verbose Logical flag controlling function messaging.
#'
#' @keywords internal

OptimizationSetup_ATSM <- function(x0, MLvec, EstType, tol = 1e-4, verbose = FALSE) {
  # --- 1) Optimization Settings ---
  max_global_iter <- 1e4
  prev_best <- -1e20
  max_iter <- 2000
  converged <- FALSE
  old_val <- mean(MLvec(x0))
  count <- 0

  # --- 2) Gradient-based optimization (with scaling if required) ---
  while (!converged) {
    if ("BFGS" %in% EstType) {
      if (count == 0) {
        D_MLvec <- Jac_approx(MLvec, x0)
        scaling_factor <- scaling_from_jacobian(D_MLvec)
      }

      MLtemp <- function(par) MLtemporary(par, scaling_factor, MLvec)

      res <- stats::optim(
        par     = x0 / scaling_factor,
        fn      = MLtemp,
        method  = "L-BFGS-B",
        control = list(maxit = max_iter, factr = 1e7) # factr ≈ reltol
      )

      x1 <- res$par * scaling_factor

      if (mean(MLvec(x1)) < mean(MLvec(x0))) x0 <- x1
    }

    # --- 3) Nelder–Mead fallback ---
    if ("Nelder-Mead" %in% EstType) {
      res <- stats::optim(
        par     = x0,
        fn      = function(par) ML_stable(par, MLvec),
        method  = "Nelder-Mead",
        control = list(maxit = 1000, reltol = 1e-8)
      )

      if (mean(MLvec(res$par)) < mean(MLvec(x0))) x0 <- res$par
    }

    # --- 4) Convergence checks ---
    new_val <- mean(MLvec(x0))
    if (verbose) {
      message(paste("*** Estimation round", count + 1, "***"))
    }

    converged <- (abs(old_val - new_val) < tol) ||
      (count > max_global_iter && new_val > prev_best)

    old_val <- new_val
    count <- count + 1
  }

  if (verbose) message("-- Done!")

  return(x0)
}

#######################################################################################################
#' Set up the vector-valued objective function (Point estimate)
#'
#' @param MLEinputs Set of inputs that are necessary to the log-likelihood function
#' @param Economies string-vector containing the names of the economies which are part of the economic system
#' @param DataFrequency  character-based vector: "Daily All Days", "Daily Business Days", "Weekly", "Monthly", "Quarterly", "Annually"
#' @param FactorLabels string-list based which contains the labels of all the variables present in the model
#' @param ModelType string-vector containing the label of the model to be estimated
#' @param BS_outType Generates simplified output list in the bootstrap setting. Default is set to FALSE.
#'
#' @returns
#' objective function
#'
#' @keywords internal

MLFunction <- function(MLEinputs, Economies, DataFrequency, FactorLabels, ModelType, BS_outType = FALSE) {
  dt <- Getdt(DataFrequency)
  mat <- MLEinputs$mat

  # 1) If one choose models in which the estimation is done country-by-country
  if (any(ModelType == c("JPS original", "JPS global", "JPS multi", "GVAR single", "GVAR multi", "JLL joint Sigma"))) {
    ML_Func <- function(...) {
      MLEdensity(...,
        r0 = NULL, MLEinputs$K0Z, MLEinputs$K1Z, se = NULL, MLEinputs$Gy.0, mat,
        MLEinputs$Y, MLEinputs$RiskFactors, MLEinputs$SpaFact, MLEinputs$Wpca, MLEinputs$We,
        MLEinputs$WpcaFull, dt, Economies, FactorLabels, ModelType, MLEinputs$GVARinputs,
        MLEinputs$JLLinputs, BS_outType
      )
    }
  } else {
    # 2) If one choose models in which the estimation is done jointly for all countries of the system with
    # the Sigma matrix estimated exclusively under the P-dynamics

    ML_Func <- function(...) {
      MLEdensity(...,
        r0 = NULL, MLEinputs$SSZ, MLEinputs$K0Z, MLEinputs$K1Z, se = NULL,
        MLEinputs$Gy.0, mat, MLEinputs$Y, MLEinputs$RiskFactors, MLEinputs$SpaFact,
        MLEinputs$Wpca, MLEinputs$We, MLEinputs$WpcaFull, dt, Economies, FactorLabels,
        ModelType, MLEinputs$GVARinputs, MLEinputs$JLLinputs, BS_outType
      )
    }
  }

  return(ML_Func)
}
###################################################################################################### "
#' Mean of the llk function used in the estimation of the selected ATSM
#'
#' @param xtemp temporary vector of parameters to be estimated numerically
#' @param scaling_vector scaling factor
#' @param MLvec log-likelihood function
#'
#' @keywords internal

MLtemporary <- function(xtemp, scaling_vector, MLvec) {
  scaled_xtemp <- scaling_vector * xtemp
  mean(MLvec(scaled_xtemp))
}


################################################################################################
#' Prevents algorithm to end up in ill-defined likelihood
#'
#' @param par vector of parameters
#' @param MLvec objective function
#'
#' @keywords internal

ML_stable <- function(par, MLvec) {
  val <- tryCatch(mean(MLvec(par)), error = function(e) NA)
  if (!is.finite(val)) {
    return(1e6) # big penalty
  }
  return(val)
}

############################################################################################################
#' Compute the time elapsed in the numerical optimization
#'
#' @param start_time Starting time
#' @param verbose Logical flag controlling function messaging.
#'
#' @keywords internal

Optimization_Time <- function(start_time, verbose) {
  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  units <- c("seconds", "minutes", "hours")
  scale <- c(1, 60, 3600)
  idx <- max(which(elapsed >= scale))
  if (verbose) message(sprintf("Elapsed time: %.2f %s", elapsed / scale[idx], units[idx]))
}


##########################################################################################################
#' Richardson extrapolation
#'
#' @param f function to differentiate
#' @param x evaluation point
#' @param i coordinate index
#' @param h base step size
#' @param fx0 f(x), precomputed
#' @param nlevels depth of extrapolation. Default is 5
#'
#' @keywords internal

richardson_diff <- function(f, x, i, h, fx0, nlevels = 5L) {
  nlevels <- as.integer(nlevels)

  # Pre-allocate as a nested list
  D <- vector("list", nlevels)
  for (n in seq_len(nlevels)) {
    D[[n]] <- vector("list", n)
  }

  for (n in seq_len(nlevels)) {
    step <- h / (2^(n - 1))
    dx <- numeric(length(x))
    dx[i] <- step

    # vector-valued central difference
    D[[n]][[1]] <- (f(x + dx) - f(x - dx)) / (2 * step)

    if (n > 1) {
      for (k in 2:n) {
        D[[n]][[k]] <- ((2^(k - 1)) * D[[n]][[k - 1]] - D[[n - 1]][[k - 1]]) / (2^(k - 1) - 1)
      }
    }
  }

  D[[nlevels]][[nlevels]] # return best Richardson estimate (vector)
}
########################################################################################## "
#'  Main Jacobian approximation
#'
#' @param f function to differentiate
#' @param x evaluation point
#' @param nlevels depth of extrapolation. Default is 5.
#'
#' @references
#' Le, A., & Singleton, K. J. (2018). Small Package of Matlab Routines for
#' Estimation of Some Term Structure Models. EABCN Training School.\cr
#' This function offers an independent R implementation that is informed
#' by the conceptual framework outlined in Le and Singleton (2018), but adapted to the
#' present modeling context.
#'
#' @keywords internal

Jac_approx <- function(f, x, nlevels = 5L) {
  h0 <- pmax(pmin(abs(x) * 1e-3, 1e-3), 1e-6) # step size per coordinate
  fx0 <- f(x)
  T_dim <- length(fx0)

  # Jacobian: nrow = #params, ncol = #outputs
  Jacobian_apx <- matrix(NA_real_, nrow = length(x), ncol = T_dim)

  for (i in seq_along(x)) {
    h <- h0[i]
    # one row is a vector of partial derivatives wrt x[i]
    Jacobian_apx[i, ] <- richardson_diff(f, x, i, h, fx0, nlevels = nlevels)
  }

  return(Jacobian_apx)
}

########################################################################################## "
#' Scaling vector computation
#'
#' @param J jacobian approximation term
#'
#' @keywords internal

scaling_from_jacobian <- function(J) {
  vec <- 1 / rowMeans(abs(J))

  # Handle infinities and zeros
  if (any(!is.finite(vec))) {
    vec[!is.finite(vec)] <- max(vec[is.finite(vec)])
  }
  if (any(vec == 0)) {
    vec[vec == 0] <- min(vec[vec > 0])
  }

  return(vec)
}
############################################################################################################
#' Adjust delta for numerical differentiation
#'
#' @param f function which contains vector (J x T) valued function handle
#' @param x parameter values
#' @param delta initial delta value
#' @param i index of the parameter being adjusted
#' @param fx0 initial function value
#' @param direction direction of adjustment (1 for positive, -1 for negative)
#'
#' @return adjusted delta value
#' @keywords internal

adjust_delta <- function(f, x, delta, i, fx0, direction) {
  fixed <- FALSE
  count <- 1
  while (!fixed && count < 10) {
    x_temp <- x
    x_temp[i] <- x[i] + direction * delta
    fx_temp <- f(x = x_temp)

    if (abs(mean((fx_temp - fx0) / delta)) > 1e8) {
      delta <- delta / 2
      count <- count + 1
    } else {
      fixed <- TRUE
    }
  }

  if (!fixed) {
    delta <- 0
  }

  return(delta)
}
