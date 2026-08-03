#' Generates the bootstrap-related outputs
#'
#' @param ModelType character. Model type to be estimated. Permissible choices: "JPS original", "JPS global", "GVAR single", "JPS multi", "GVAR multi", "JLL original", "JLL No DomUnit", "JLL joint Sigma".
#' @param ModelParaPE list. Point estimates of the model parameters. See outputs from \code{\link{Optimization}}.
#' @param NumOutPE list. Point estimates from numerical outputs. See outputs from \code{\link{NumOutputs}}.
#' @param Economies character vector. Names of the \code{C} economies included in the system.
#' @param InputsForOutputs list. Inputs for generating IRFs, GIRFs, FEVDs, GFEVDs, and Term Premia.
#' @param FactorLabels list. Labels for all variables present in the model, as returned by \code{\link{LabFac}}.
#' @param JLLlist list. Inputs for JLL model estimation (see \code{\link{JLL}}). Default is NULL.
#' @param GVARlist list. Inputs for GVAR model estimation (see \code{\link{GVAR}}). Default is NULL.
#' @param WishBC logical. Whether to estimate the physical parameter model with bias correction (see \code{\link{Bias_Correc_VAR}}). Default is FALSE.
#' @param BRWlist list. Inputs for bias-corrected estimation (see \code{\link{Bias_Correc_VAR}}).
#' @param Folder2save character. Folder path where outputs will be stored. Default saves outputs in a temporary directory.
#' @param verbose logical. Print progress messages. Default is TRUE.
#'
#' @section Permissible options - Bootstrap list in \code{InputsForOutputs}:
#' \itemize{
#'    \item \code{methodBS} : \code{"bs"} (standard bootstrap), \code{"wild"} (wild bootstrap), \code{"block"} (block bootstrap)
#'    \item \code{BlockLength} : required input for the block bootstrap method. Block length must be larger than 0 and smallar than the model time series dimension (\code{Td}).
#'    \item \code{ndraws}: number of draws. Must be a positive integer.
#'    \item \code{pctg} : confidence level. Must be a positive integer. Common choices are: 68, 90 and  95.
#' }
#'
#' @examples
#' \donttest{
#' data("ParaSetEx")
#' data("InpForOutEx")
#' data("NumOutEx")
#' ModelType <- "JPS original"
#' Economy <- "Brazil"
#' FacLab <- LabFac(N = 1, DomVar = "Eco_Act", GlobalVar = "Gl_Eco_Act", Economy, ModelType)
#'
#' # Adjust Forecasting setting
#' InpForOutEx[[ModelType]]$Bootstrap <- list(
#'   WishBootstrap = 1, methodBS = "bs", BlockLength = 4,
#'   ndraws = 5, pctg = 95
#' )
#'
#' Boot <- Bootstrap(ModelType, ParaSetEx, NumOutEx, Economy, InpForOutEx, FacLab,
#'   JLLlist = NULL,
#'   GVARlist = NULL, WishBC = FALSE, BRWlist = NULL, Folder2save = NULL, verbose = FALSE
#' )
#' }
#'
#' @section Available methods:
#' - \code{autoplot(object, NumOutPE, type)}
#'
#' @returns
#' An object of class 'ATSMModelBoot' containing:
#' \itemize{
#'   \item List of model parameters for each draw
#'   \item List of numerical outputs (IRFs, GIRFs, FEVDs, GFEVDs) for each draw
#'   \item Confidence bounds for the chosen level of significance
#' }
#'
#' @export

Bootstrap <- function(ModelType, ModelParaPE, NumOutPE, Economies, InputsForOutputs, FactorLabels, JLLlist,
                      GVARlist, WishBC = FALSE, BRWlist = NULL, Folder2save = NULL, verbose = TRUE) {
  if (verbose) message("3) BOOTSTRAP ANALYSIS")
  WishBoot <- InputsForOutputs[[ModelType]]$Bootstrap$WishBoot

  if (!WishBoot) {
    if (verbose) message("No Bootstrap analysis was generated \n")
    return(NULL)
  }

  if (verbose) message("3.1) Estimating bootstrap setup. This may take several hours.")

  StatQ <- InputsForOutputs$StationaryQ
  UMatY <- InputsForOutputs$UnitMatYields

  # 1) Pre-allocation of list of outputs
  ndraws <- InputsForOutputs[[ModelType]]$Bootstrap$ndraws
  DataFreq <- InputsForOutputs$DataFreq
  N <- length(FactorLabels$Spanned)

  dt <- Getdt(DataFreq)
  mat <- if (any(ModelType == c("JPS original", "JPS global", "GVAR single"))) {
    ModelParaPE[[ModelType]][[Economies[1]]]$Inputs$mat
  } else {
    ModelParaPE[[ModelType]]$Inputs$mat
  }

  ModelBootstrap <- list(GeneralInputs = list(mat = mat, dt = dt, N = N))

  # 2) Obtain the residuals from the original model
  # a) P-dynamics residuals
  residPdynOriginal <- PdynResid_BS(ModelType, Economies, ModelParaPE)

  # b) Yield residuals
  BFull_Original <- Get_BFull(ModelParaPE, FactorLabels, mat, Economies, ModelType)
  residYieOriginal <- residY_original(ModelParaPE, BFull_Original, ModelType, Economies)

  # 3) Bootstrap samples
  tt <- 1 # numbers of accepted draws
  ww <- 1 # index for printing on screen

  start_time <- Sys.time()

  while (tt <= ndraws) {
    if (tt == 10 * ww) {
      if (verbose) message(paste("Loop ", tt, " / ", ndraws, " draws"))
      ww <- ww + 1
    }

    # a) generate the artificial data
    invisible(utils::capture.output(Series_artificial <- Gen_Artificial_Series(ModelParaPE, residPdynOriginal,
      residYieOriginal, ModelType, BFull_Original,
      InputsForOutputs, Economies, FactorLabels,
      GVARlist, JLLlist, WishBC, BRWlist,
      verbose = F
    )))

    # b) Prepare the inputs for model estimation
    Y_BS <- t(Series_artificial$Y_BS)
    Global_BS <- t(Series_artificial$GlobalMacro_BS)
    Dom_BS <- t(Series_artificial$DomesticMacro_BS)
    t0_BS <- colnames(Global_BS)[1]
    tF_BS <- utils::tail(colnames(Global_BS), 1)

    invisible(utils::capture.output(ATSMInputs_BS <- InputsForOpt(t0_BS, tF_BS, ModelType, Y_BS, Global_BS, Dom_BS,
      FactorLabels, Economies, DataFreq, GVARlist, JLLlist,
      WishBC, BRWlist, UMatY,
      CheckInputs = F,
      BS_Adj = T, verbose = F
    )))

    # c) Run the optimization
    invisible(utils::capture.output(Draw_Opt <- Optimization(ATSMInputs_BS, StatQ, DataFreq, FactorLabels,
      Economies, ModelType,
      tol = 1e-1, EstType = "Nelder-Mead",
      TimeCount = FALSE, BS_outputs = TRUE, verbose = FALSE
    )))

    ModelBootstrap <- AdjustOptm_BS(ModelType, ModelBootstrap, Draw_Opt, Economies, tt)

    # if the optimization crashes after a particular draws, we can still keep the outputs of draws before
    FolderPath <- if (is.null(Folder2save)) tempdir() else Folder2save
    saveRDS(ModelBootstrap, paste(FolderPath, "/Bootstrap_", InputsForOutputs$"Label Outputs", ".rds", sep = ""))
    tt <- tt + 1
  }

  if (verbose) message("-- Done!")
  Optimization_Time(start_time, verbose)

  # 4) Compute the numerical outputs from the bootstrap samples
  if (verbose) message("3.2) Computing numerical outputs.")
  ModelBootstrap$NumOutDraws <- NumOutputs_Bootstrap(
    ModelType, ModelBootstrap, InputsForOutputs,
    FactorLabels, Economies
  )

  # 5) Compute confidence bounds
  if (verbose) message("3.3) Computing confidence bounds and producing graphical outputs.")
  ModelBootstrap$ConfBounds <- BootstrapBoundsSet(
    ModelType, ModelBootstrap, NumOutPE, InputsForOutputs,
    Economies, FolderPath, verbose
  )

  # To save space, clean the repeated outputs from the JLL outputs
  if (any(ModelType == c("JLL original", "JLL No DomUnit", "JLL joint Sigma"))) {
    ModelBootstrap <- CleanOrthoJLL_Boot(ModelBootstrap, ndraws, ModelType)
  }

  saveRDS(ModelBootstrap, paste(FolderPath, "/Bootstrap_", InputsForOutputs$"Label Outputs", ".rds", sep = ""))

  return(structure(ModelBootstrap, class = "ATSMModelBoot"))
}

################################################################################################################
################################################################################################################
#' Compute some key parameters from the P-dynamics (Bootstrap set)
#'
#' @param ModelType A character vector containing the label of the model to be estimated
#' @param Economies A character vector containing the names of the economies included in the system.
#' @param ModelPara_PE Point estimate from the model parameters
#'
#' @keywords internal

PdynResid_BS <- function(ModelType, Economies, ModelPara_PE) {
  Get_residuals <- function(ZZ, K0Z, K1Z) {
    T_dim <- ncol(ZZ)
    t(ZZ[, 2:T_dim] - matrix(K0Z, nrow = nrow(K0Z), ncol = T_dim - 1) - K1Z %*% ZZ[, 1:(T_dim - 1)])
  }

  # SepQ models
  if (ModelType %in% c("JPS original", "JPS global", "GVAR single")) {
    eZ <- stats::setNames(lapply(Economies, function(economy) {
      Get_residuals(
        ModelPara_PE[[ModelType]][[economy]]$Inputs$AllFactors,
        ModelPara_PE[[ModelType]][[economy]]$ModEst$P$K0Z,
        ModelPara_PE[[ModelType]][[economy]]$ModEst$P$K1Z
      )
    }), Economies)
  } else {
    # JointQ models
    eZ <- Get_residuals(
      ModelPara_PE[[ModelType]]$Inputs$AllFactors,
      ModelPara_PE[[ModelType]]$ModEst$P$K0Z,
      ModelPara_PE[[ModelType]]$ModEst$P$K1Z
    )
  }

  return(eZ)
}
################################################################################################################
#' Compute the residuals from the original model
#'
#' @param residPdynOriginal Time-series of the residuals from the P-dynamics equation (T x F)
#' @param residYieOriginal Time-series of the residuals from the observational equation (T x J or T x CJ)
#' @param InputsForOutputs List containing the desired inputs for the construction of the numerical outputs.
#' @param ModelType A character vector indicating the model type to be estimated
#' @param nlag Number of lags in the P-dynamics. Default is set to 1.
#'
#' @keywords internal

ResampleResiduals_BS <- function(residPdynOriginal, residYieOriginal, InputsForOutputs, ModelType, nlag = 1) {
  methodBS <- InputsForOutputs[[ModelType]]$Bootstrap$methodBS
  BlockLength <- InputsForOutputs[[ModelType]]$Bootstrap$BlockLength

  T_dim <- nrow(residYieOriginal)
  K <- ncol(residPdynOriginal)

  # Define a function to apply bootstrap resampling
  # a) Residuals to bootstrap:
  if (methodBS == "bs") {
    Rand <- matrix(stats::runif(T_dim, min = 0, max = 1), ncol = 1)
    rr <- ceiling((T_dim - nlag) * Rand)
    uPdyn <- residPdynOriginal[rr[1:(T_dim - nlag)], ]
    uYiel <- residYieOriginal[rr, ]
  } else if (methodBS == "wild") {
    # b) Wild bootstrap based on simple distribution (~Rademacher)
    Rand <- matrix(stats::runif(T_dim, min = 0, max = 1), ncol = 1)
    rr <- 1 - 2 * (Rand > 0.5)
    uPdyn <- residPdynOriginal * (rr[1:(T_dim - nlag)] %*% matrix(1, nrow = 1, ncol = K))
    uYiel <- residYieOriginal * (rr %*% matrix(1, nrow = 1, ncol = ncol(residYieOriginal)))
  } else if (methodBS == "block") {
    # c) Blocks overlap and are drawn with replacement
    FullBlocksSet <- dim(residPdynOriginal)[1] - BlockLength + 1 # all possible blocks that can be drawn
    SampleBlock <- ceiling((T_dim - nlag) / BlockLength)

    Rand <- matrix(stats::runif(SampleBlock, min = 0, max = 1), ncol = 1)
    bb <- ceiling(SampleBlock * Rand)
    IdxBlocks <- matrix(NA, nrow = BlockLength, ncol = FullBlocksSet)
    for (mm in 1:FullBlocksSet) {
      IdxBlocks[, mm] <- mm:(mm + BlockLength - 1)
    }
    rr <- as.vector(IdxBlocks[, bb])[1:T_dim]
    uPdyn <- residPdynOriginal[rr[1:(T_dim - nlag)], ]
    uYiel <- residYieOriginal[rr, ]
  } else {
    stop(paste("The method ", methodBS, " is not available"))
  }

  return(list(residFact = uPdyn, residYields = uYiel))
}

##############################################################################################################
#' Compute the residuals from the observational equation
#'
#' @param ModelParaPE List of point estimates of the model parameter
#' @param BFull Matrix B of loadings (CJ x F or J x F)
#' @param ModelType A character vector indicating the model type to be estimated
#' @param Economies String-vector containing the names of the economies which are part of the economic system
#'
#' @keywords internal

residY_original <- function(ModelParaPE, BFull, ModelType, Economies) {
  sepQ_Labels <- ModelType %in% c("JPS original", "JPS global", "GVAR single")

  compute_residuals <- function(YY, ZZ, A, B) {
    t(YY - (matrix(A, nrow = nrow(YY), ncol = ncol(YY)) + B %*% ZZ))
  }

  if (sepQ_Labels) {
    residYie <- stats::setNames(lapply(Economies, function(economy) {
      params <- ModelParaPE[[ModelType]][[economy]]
      compute_residuals(
        params$Inputs$Y, params$Inputs$AllFactors,
        params$ModEst$Q$Load$P$A, BFull[[economy]]
      )
    }), Economies)
  } else {
    params <- ModelParaPE[[ModelType]]
    residYie <- compute_residuals(
      params$Inputs$Y, params$Inputs$AllFactors,
      params$ModEst$Q$Load$P$A, BFull
    )
  }

  return(residYie)
}

######################################################################################################### ""
#' Compute the B matrix of loadings
#'
#' @param ModelParaPE List of point estimates of the model parameter
#' @param FactorLabels String-list based which contains the labels of all the variables present in the model
#' @param mat Vector of bond yield maturities
#' @param Economies String-vector containing the names of the economies which are part of the economic system
#' @param ModelType A character vector indicating the model type to be estimated
#'
#' @keywords internal

Get_BFull <- function(ModelParaPE, FactorLabels, mat, Economies, ModelType) {
  J <- length(mat)
  G <- length(FactorLabels$Global)
  N <- length(FactorLabels$Spanned)
  M <- length(FactorLabels$Domestic) - N
  C <- length(Economies)

  # For models estimated separately
  if (any(ModelType == c("JPS original", "JPS global", "GVAR single"))) {
    K <- nrow(ModelParaPE[[ModelType]][[Economies[1]]]$Inputs$AllFactors)
    BFull <- lapply(Economies, function(economy) {
      B_CS <- matrix(0, nrow = J, ncol = K)
      LabelSpannedCS <- c(FactorLabels$Tables[[economy]][-(1:M)])
      RiskFactorLabels <- rownames(ModelParaPE[[ModelType]][[economy]]$Inputs$AllFactors)
      idxSpanned <- match(LabelSpannedCS, RiskFactorLabels)
      B <- ModelParaPE[[ModelType]][[economy]]$ModEst$Q$Load$P$B
      B_CS[, idxSpanned] <- B
      B_CS
    })
    names(BFull) <- Economies
  } else {
    # For models estimated jointly
    K <- nrow(ModelParaPE[[ModelType]]$Inputs$AllFactors)
    B <- ModelParaPE[[ModelType]]$ModEst$Q$Load$P$B
    BFull <- BUnspannedAdapJoint(G, M, N, C, J, B)
  }

  return(BFull)
}
###################################################################################################################
#' Build the time-series of the risk factors in each bootstrap draw
#'
#' @param ModelParaPE List of point estimates of the model parameter
#' @param residPdynOriginal Time-series of the residuals from the P-dynamics equation (T x F)
#' @param residYieOriginal Time-series of the residuals from the observational equation (T x J or T x CJ)
#' @param InputsForOutputs List containing the desired inputs for the construction of the numerical outputs.
#' @param Economies String-vector containing the names of the economies which are part of the economic system
#' @param ModelType Desired model to be estimated
#' @param FactorLabels String-list based which contains the labels of all the variables present in the model
#' @param GVARlist List of necessary inputs for the estimation of GVAR-based models
#' @param JLLlist List of necessary inputs for the estimation of JLL-based models
#' @param WishBRW Whether the user wishes to estimate the physical parameter model with the Bias correction model from BRW (2012) (see "Bias_Correc_VAR" function). Default is set to 0.
#' @param BRWlist List of necessary inputs for performing the bias-corrected estimation (see "Bias_Correc_VAR" function)
#' @param nlag Number of lags in the P-dynamics. Default is set to 1.
#' @param verbose Logical flag controlling function messaging.
#'
#' @keywords internal

BuildRiskFactors_BS <- function(ModelParaPE, residPdynOriginal, residYieOriginal, InputsForOutputs, Economies,
                                ModelType, FactorLabels, GVARlist, JLLlist, WishBRW, BRWlist, nlag = 1, verbose) {
  sepQ_Labels <- ModelType %in% c("JPS original", "JPS global", "GVAR single")
  multiQ_Labels <- ModelType %in% c("JPS multi", "GVAR multi", "JLL original", "JLL No DomUnit", "JLL joint Sigma")

  # 1) Initialization
  if (sepQ_Labels) {
    ZZ_list <- list()
    resid_list <- list()
  }


  for (i in 1:length(Economies)) {
    if (multiQ_Labels & i > 1) break

    MaxEigen <- 1.1 # Initialization Max eigenvalues
    while (MaxEigen > 1) { # Test whether the VAR is stationary (if not, drop the draw)

      # Extract model parameters
      if (sepQ_Labels) {
        params <- ModelParaPE[[ModelType]][[Economies[i]]]
        resids_BS <- ResampleResiduals_BS(
          residPdynOriginal[[Economies[i]]], residYieOriginal[[Economies[i]]],
          InputsForOutputs, ModelType
        )
      } else {
        params <- ModelParaPE[[ModelType]]
        resids_BS <- ResampleResiduals_BS(residPdynOriginal, residYieOriginal, InputsForOutputs, ModelType)
      }

      RiskFactors <- params$Inputs$AllFactors
      T_dim <- ncol(RiskFactors)
      K <- nrow(RiskFactors)

      ZZ_Boot <- matrix(NA,
        nrow = T_dim - 1 + nlag, ncol = K,
        dimnames = list(rownames(t(RiskFactors)), colnames(t(RiskFactors)))
      )

      # 2)  Compute artificial time-series
      # 2.1) initial values for the artificial data
      ZZ_Boot[1:nlag, ] <- t(RiskFactors)[1:nlag, ] # Initial values
      LAG <- ZZ_Boot[1:nlag, , drop = FALSE]
      LAGplus <- cbind(1, LAG)

      # 2.2) generate artificial series
      Ft <- rbind(t(params$ModEst$P$K0Z), t(params$ModEst$P$K1Z))
      # From observation nlag+1 to nobs, compute the artificial data
      for (jj in (nlag + 1):(T_dim - 1 + nlag)) {
        for (mm in 1:K) {
          ZZ_Boot[jj, mm] <- LAGplus %*% as.matrix(Ft[, mm]) + resids_BS$residFact[jj - nlag, mm]
        }
        # Update the LAG matrix
        if (jj < T_dim - 1 + nlag) {
          LAG <- rbind(ZZ_Boot[jj, ], LAG[1, seq_len((nlag - 1) * K)])
          LAGplus <- cbind(1, LAG)
        }
      }

      # 3) Test whether the VAR is stationary (if not, drop the draw)
      if (multiQ_Labels) {
        K1Z_BS <- FeedbackMat_BS(ModelType, t(ZZ_Boot), FactorLabels, Economies, GVARlist, JLLlist, WishBRW, BRWlist, verbose)
      } else {
        Economies_temp <- if (ModelType %in% c("JPS original", "JPS global")) Economies[i] else Economies
        RiskFact_Temp <- stats::setNames(list(t(ZZ_Boot)), Economies[i])
        K1Z_BS <- FeedbackMat_BS(
          ModelType, RiskFact_Temp, FactorLabels, Economies_temp,
          GVARlist, JLLlist, WishBRW, BRWlist, verbose
        )
      }
      MaxEigen <- max(abs(eigen(K1Z_BS)$values))
    }

    # 4) Store outputs to export for sepQ models
    if (sepQ_Labels) {
      ZZ_list[[Economies[i]]] <- ZZ_Boot
      resid_list[[Economies[i]]] <- resids_BS
    }
  }

  return(if (sepQ_Labels) list(ZZ_BS = ZZ_list, resids_BS = resid_list) else list(ZZ_BS = ZZ_Boot, resids_BS = resids_BS))
}
###################################################################################################################
#' Generate artificial time-series in the bootstrap setup
#'
#' @param ModelParaPE List of point estimates of the model parameter
#' @param residPdynOriginal Time-series of the residuals from the P-dynamics equation (T x F)
#' @param residYieOriginal Time-series of the residuals from the observational equation (T x J or T x CJ)
#' @param ModelType Desired model to be estimated
#' @param BFull Matrix B of loadings (CJ x F or J x F)
#' @param InputsForOutputs List containing the desired inputs for the construction
#' @param Economies String-vector containing the names of the economies which are part of the economic system
#' @param FactorLabels String-list based which contains the labels of all the variables present in the model
#' @param GVARlist List of necessary inputs for the estimation of GVAR-based models
#' @param JLLlist List of necessary inputs for the estimation of JLL-based models
#' @param WishBRW Whether the user wishes to estimate the physical parameter model with the Bias correction model from BRW (2012) (see "Bias_Correc_VAR" function). Default is set to 0.
#' @param BRWlist List of necessary inputs for performing the bias-corrected estimation (see "Bias_Correc_VAR" function)
#' @param verbose Logical flag controlling function messaging.
#' @param nlag Number of lags in the P-dynamics. Default is set to 1.
#'
#' @keywords internal

Gen_Artificial_Series <- function(ModelParaPE, residPdynOriginal, residYieOriginal, ModelType, BFull,
                                  InputsForOutputs, Economies, FactorLabels, GVARlist, JLLlist, WishBRW, BRWlist,
                                  verbose, nlag = 1) {
  # 1) Artificial time-series of the risk factors
  BS_Set <- BuildRiskFactors_BS(
    ModelParaPE, residPdynOriginal, residYieOriginal, InputsForOutputs, Economies,
    ModelType, FactorLabels, GVARlist, JLLlist, WishBRW, BRWlist, nlag, verbose
  )

  # Extract Unspanned factors from the full risk factor set
  ZZ_BS <- BS_Set$ZZ_BS

  if (any(ModelType == c("JPS original", "JPS global", "GVAR single"))) {
    G <- length(FactorLabels$Global)
    N <- length(FactorLabels$Spanned)

    UnspannedFactors_CS_BS <- function(Economy, G, N) {
      AllFactors <- BS_Set$ZZ_BS[[Economy]]
      AllFactors[, (G + 1):(ncol(AllFactors) - N), drop = FALSE]
    }

    DomesticMacro_BS <- do.call(cbind, lapply(Economies, UnspannedFactors_CS_BS, G, N))
    GlobalMacro_BS <- do.call(cbind, lapply(Economies, function(country) ZZ_BS[[country]][, seq_len(G), drop = FALSE]))
  } else {
    Idxs <- Idx_UnspanFact(t(ZZ_BS), FactorLabels, Economies)
    GlobalMacro_BS <- ZZ_BS[, Idxs$IdxGlobal, drop = FALSE]
    DomesticMacro_BS <- ZZ_BS[, Idxs$IdxunSpa, drop = FALSE]
  }

  Y_BS <- BuildYields_BS(ModelParaPE, ModelType, ZZ_BS, BFull, BS_Set, Economies)

  return(list(
    ZZ_BS = ZZ_BS, GlobalMacro_BS = GlobalMacro_BS, DomesticMacro_BS = DomesticMacro_BS,
    Y_BS = Y_BS
  ))
}
#################################################################################################################
#' Clean unnecessary outputs of JLL models in the bootstrap setup
#'
#' @param ModelBootstrap List of outputs to store bootstrap draws
#' @param ndraws Total number of bootstrap draws
#' @param ModelType A character vector indicating the model type to be estimated
#'
#' @keywords internal

CleanOrthoJLL_Boot <- function(ModelBootstrap, ndraws, ModelType) {
  # a) Bootstrap
  # IRF
  ModelBootstrap$NumOutDraws$IRF[[ModelType]] <- lapply(
    ModelBootstrap$NumOutDraws$IRF[[ModelType]], function(draw) {
      if (!is.null(draw)) draw$Yields$Ortho <- NULL
      return(draw)
    }
  )

  # FEVD
  ModelBootstrap$NumOutDraws$FEVD[[ModelType]] <- lapply(
    ModelBootstrap$NumOutDraws$FEVD[[ModelType]], function(draw) {
      if (!is.null(draw)) draw$Yields$Ortho <- NULL
      return(draw)
    }
  )

  # b) Point Estimates
  ModelBootstrap$ConfBounds$IRF[[ModelType]]$Yields$Ortho <- NULL
  ModelBootstrap$ConfBounds$FEVD[[ModelType]]$Yields$Ortho <- NULL

  return(ModelBootstrap)
}

##############################################################################################################
#' Compute the Feedback matrix of each bootstrap draw
#'
#' @param ModelType String-vector containing the label of the model to be estimated
#' @param RiskFactors_TS Time-series of risk factors of the bootstrap (F x T)
#' @param FactorLabels String-list based which contains the labels of all the variables present in the model
#' @param Economies String-vector containing the names of the economies which are part of the economic system
#' @param GVARlist List of necessary inputs for the estimation of GVAR-based models
#' @param JLLlist List of necessary inputs for the estimation of JLL-based models
#' @param WishBRW Whether the user wishes to estimate the physical parameter model with the Bias correction model from BRW (2012) (see "Bias_Correc_VAR" function). Default is set to 0.
#' @param BRWlist List of necessary inputs for performing the bias-corrected estimation (see \code{\link{Bias_Correc_VAR}} function)
#' @param verbose Logical flag controlling function messaging.
#'
#' @keywords internal

FeedbackMat_BS <- function(ModelType, RiskFactors_TS, FactorLabels, Economies, GVARlist, JLLlist,
                           WishBRW, BRWlist, verbose) {
  # Model-specific inputs
  SpeInputs <- SpecificMLEInputs(
    ModelType, Economies, RiskFactors_TS, FactorLabels, GVARlist, JLLlist,
    WishBRW, BRWlist
  )

  if (WishBRW) SpeInputs$BRWinputs[c("checkBRW", "checkSigma")] <- 0

  # 1) Two special cases
  # a) JLL models without bias correction (avoid the unnecessary numerical optimization from the Sigma matrix)
  if (ModelType %in% c("JLL original", "JLL No DomUnit", "JLL joint Sigma") & !WishBRW) {
    SpeInputs$JLLinputs$WishSigmas <- 0
    N <- length(FactorLabels$Spanned)
    PdynPara <- JLL(RiskFactors_TS, N, SpeInputs$JLLinputs)
    K1Z_BS <- PdynPara$k1

    # b) GVAR single model and  bias correction
  } else if (ModelType == "GVAR single" & WishBRW) {
    N <- length(FactorLabels$Spanned)
    PdynPara <- Bias_Correc_VAR(
      ModelType, SpeInputs$BRWinputs, t(RiskFactors_TS[[1]]), Economies, FactorLabels,
      SpeInputs$GVARinputs
    )

    K1Z_BS <- PdynPara$K1Z_BC
  } else {
    # 2) All other specifications
    PdynPara <- GetPdynPara(RiskFactors_TS, FactorLabels, Economies, ModelType, SpeInputs$BRWinputs,
      SpeInputs$GVARinputs, SpeInputs$JLLinputs,
      verbose = verbose
    )

    if (ModelType %in% c("JPS original", "JPS global", "GVAR single")) {
      K1Z_BS <- PdynPara[[Economies[1]]]$K1Z
    } else {
      K1Z_BS <- PdynPara$K1Z
    }
  }

  return(K1Z_BS)
}

################################################################################################################
#' Build the time-series of bond yields for each bootstrap draw
#'
#' @param ModelParaPE List of point estimates of the model parameter
#' @param ModelType String-vector containing the label of the model to be estimated
#' @param RiskFactors_BS Time-series of the risk factors (F x T)
#' @param BFull B matrix of loadings
#' @param BS_Set Set of bootstrap inputs
#' @param Economies String-vector containing the names of the economies which are part of the economic system
#'
#' @keywords internal

BuildYields_BS <- function(ModelParaPE, ModelType, RiskFactors_BS, BFull, BS_Set, Economies) {
  # Models estimated jointly
  if (any(ModelType == c("JPS original", "JPS global", "GVAR single"))) {
    T_dim <- nrow(RiskFactors_BS[[Economies[1]]])
    TS_Labels <- rownames(RiskFactors_BS[[Economies[1]]])
    Y_listBS <- list()

    for (i in seq_along(Economies)) {
      YieldsLabels <- rownames(ModelParaPE[[ModelType]][[Economies[i]]]$Inputs$Y)

      Y_CS <- matrix(NA, nrow = T_dim, ncol = length(YieldsLabels))
      dimnames(Y_CS) <- list(TS_Labels, YieldsLabels)

      A <- ModelParaPE[[ModelType]][[Economies[i]]]$ModEst$Q$Load$P$A
      ZZ_BS <- RiskFactors_BS[[Economies[i]]]
      Y_CS <- BS_Set$resids_BS[[Economies[i]]]$residYields + matrix(A, nrow = T_dim, ncol = length(YieldsLabels), byrow = TRUE) + ZZ_BS %*% t(BFull[[Economies[i]]])
      rownames(Y_CS) <- TS_Labels

      Y_listBS[[Economies[i]]] <- Y_CS
    }

    Y_BS <- do.call(cbind, lapply(seq_along(Economies), function(i) {
      Y_listBS[[Economies[i]]]
    }))

    # Models estimated separately
  } else {
    T_dim <- nrow(RiskFactors_BS)
    TS_Labels <- rownames(RiskFactors_BS)
    YieldsLabels <- rownames(ModelParaPE[[ModelType]]$Inputs$Y)

    Y_BS <- matrix(NA, nrow = T_dim, ncol = length(YieldsLabels))
    dimnames(Y_BS) <- list(TS_Labels, YieldsLabels)

    A <- ModelParaPE[[ModelType]]$ModEst$Q$Load$P$A
    Y_BS <- BS_Set$resids_BS$residYields + matrix(A, nrow = T_dim, ncol = length(YieldsLabels), byrow = T_dim) + RiskFactors_BS %*% t(BFull)
    rownames(Y_BS) <- TS_Labels
  }

  return(Y_BS)
}

##############################################################################################################
#' Gathers the estimate of the bootstrap draws
#'
#' @param ModelType String-vector containing the label of the model to be estimated
#' @param ModelBootstrap List to store the bootstrap set
#' @param Draw_Opt List of model estimated parameters
#' @param Economies String-vector containing the names of the economies which are part of the economic system
#' @param tt Number of the bootstrap draw
#'
#' @keywords internal

AdjustOptm_BS <- function(ModelType, ModelBootstrap, Draw_Opt, Economies, tt) {
  if (ModelType %in% c("JPS original", "JPS global", "GVAR single")) {
    for (i in seq_along(Economies)) {
      ModelBootstrap$ParaDraws[[ModelType]][[Economies[i]]][[tt]] <- Draw_Opt[[ModelType]][[Economies[i]]]
    }
  } else {
    ModelBootstrap$ParaDraws[[ModelType]][[tt]] <- Draw_Opt[[ModelType]]
  }

  return(ModelBootstrap)
}

#############################################################################################################
#' Prepare the factor set for GVAR models (Bootstrap version)
#'
#' @param ModelType A character vector containing the label of the model to be estimated.
#' @param RiskFactors A matrix of the complete set of risk factors (K x T).
#' @param Wgvar A transition matrix from GVAR models (C x C).
#' @param Economies A character vector containing the names of the economies included in the system.
#' @param FactorLabels A list of character vectors with labels for all variables in the model.
#'
#' @return A list containing the factor set for GVAR models.
#' @keywords internal

DataSet_BS <- function(ModelType, RiskFactors, Wgvar, Economies, FactorLabels) {
  if (!any(ModelType %in% c("GVAR single", "GVAR multi"))) {
    return(NULL)
  }

  # 1) Extract dimensions
  T_dim <- ncol(RiskFactors) # Time dimension
  C <- length(Economies) # Number of economies
  N <- length(FactorLabels$Spanned)
  M <- length(FactorLabels$Domestic) - N
  M_star <- length(FactorLabels$Star) - N
  G <- length(FactorLabels$Global)

  # 2) Initialize lists
  ListFactors <- stats::setNames(vector("list", C + 1), c(Economies, "Global"))
  CSF <- stats::setNames(vector("list", length(FactorLabels$Domestic)), FactorLabels$Domestic)
  SF <- stats::setNames(vector("list", length(FactorLabels$Star)), FactorLabels$Star)

  # 3) Merge country-specific lists
  ListFactors[Economies] <- lapply(Economies, function(e) list(Factors = append(CSF, SF)))

  # 4) Assign global factors
  ListFactors$Global <- stats::setNames(lapply(FactorLabels$Global, function(f) as.matrix(RiskFactors[f, ])), FactorLabels$Global)

  # 5) Fill country-specific factors
  for (economy in Economies) {
    ListFactors[[economy]] <- list(
      Factors = stats::setNames(
        Map(function(j) as.matrix(RiskFactors[FactorLabels$Tables[[economy]][j], ]), seq_len(M + N)),
        FactorLabels$Domestic
      )
    )
  }

  # 6) Compute star variables
  Z <- stats::setNames(lapply(seq_len(M + N), function(j) {
    matrix(
      unlist(lapply(Economies, function(e) ListFactors[[e]]$Factors[[j]])),
      nrow = length(Economies), byrow = TRUE
    )
  }), FactorLabels$Domestic)

  # 6.1) If star variables are computed with time-varying weights
  if (is.list(Wgvar)) {
    ListFactors <- TimeVarWeights_GVAR(RiskFactors, Economies, Z, ListFactors, Wgvar, FactorLabels)
  } else {
    # 6.2) Compute star variables with fixed weights
    idx1 <- M + N
    for (i in seq_len(C)) {
      for (j in seq_len(M + N)) {
        ListFactors[[Economies[i]]]$Factors[[idx1 + j]] <- t(Wgvar[i, ] %*% Z[[j]])
      }
      names(ListFactors[[Economies[i]]]$Factors) <- c(FactorLabels$Domestic, FactorLabels$Star)
    }
  }
  return(ListFactors)
}

#################################################################################################
#' Compute the star variables with time-varying weights
#'
#' @param RiskFactors A matrix of the complete set of risk factors (F x T).
#' @param Economies A character vector containing the names of the economies included in the system.
#' @param RiskFactors_List List of domestic risk factors (both spanned and unspanned)
#' @param ListFactors List of risk factors
#' @param Wgvar List of transition matrices
#' @param FactorLabels A list of character vectors with labels for all variables in the model.
#'
#' @keywords internal

TimeVarWeights_GVAR <- function(RiskFactors, Economies, RiskFactors_List, ListFactors, Wgvar, FactorLabels) {
  T_dim <- ncol(RiskFactors)
  K <- length(RiskFactors_List)

  # Use only the transition matrices that are included in the sample span
  t_First <- as.Date(colnames(RiskFactors)[1], format = "%d-%m-%Y")
  t_Last <- as.Date(colnames(RiskFactors)[T_dim], format = "%d-%m-%Y")

  t_First_Wgvar <- format(t_First, "%Y")
  t_Last_Wgvar <- format(t_Last, "%Y")

  Wgvar_subset <- Wgvar[names(Wgvar) >= t_First_Wgvar & names(Wgvar) <= t_Last_Wgvar]

  # Add common column label (i.e. the year of the observation) to all variables
  Dates <- as.Date(colnames(RiskFactors), format = "%d-%m-%Y")
  YearLabels <- substr(Dates, 1, 4)
  Z <- lapply(RiskFactors_List, function(x) {
    colnames(x) <- paste0(YearLabels, seq_along(colnames(x)))
    return(x)
  })


  # Compute the star variables with time-varying weights
  for (i in 1:length(Economies)) {
    for (j in 1:K) {
      StarTimeVarTemp <- matrix(NA, nrow = T_dim, ncol = 1)

      for (k in 1:length(Wgvar_subset)) {
        YearRef <- names(Wgvar_subset)[k] # year of reference
        IdxYear <- grep(YearRef, colnames(Z[[j]]))
        WgvarYear <- Wgvar_subset[[k]]
        StarTimeVarTemp[IdxYear] <- t(WgvarYear[i, ] %*% Z[[j]][, IdxYear])
      }
      # If the last year of the transition matrix happens earlier than the year of the last observation from the sample,
      # then use the last transition matrix for the remaining observations
      if (anyNA(StarTimeVarTemp)) {
        LenlastYear <- length(IdxYear)
        IdxLastObs <- IdxYear[LenlastYear] + 1
        StarTimeVarTemp[IdxLastObs:T_dim] <- t(WgvarYear[i, ] %*% Z[[j]][, (IdxLastObs):T_dim])
      }

      ListFactors[[Economies[i]]]$Factors[[K + j]] <- StarTimeVarTemp
    }
    names(ListFactors[[Economies[i]]]$Factors) <- c(FactorLabels$Domestic, FactorLabels$Star)
  }

  return(ListFactors)
}
