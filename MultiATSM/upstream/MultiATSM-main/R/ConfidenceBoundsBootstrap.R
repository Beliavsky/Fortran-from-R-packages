#' Builds the confidence bounds and graphs (Bootstrap set)
#'
#' @param ModelType string-vector containing the label of the model to be estimated
#' @param ModelBootstrap list containing the complete set of model parameters after the bootstrap estimation procedure
#' @param NumOutPE point estimate from the numerical outputs (see the outputs of the "NumOutputs" function)
#' @param InputsForOutputs list containing the desired inputs for the construction of IRFs, GIRFs, FEVDs, and GFEVDs
#' @param Economies string-vector containing the names of the economies which are part of the economic system
#' @param Folder2save Folder path where the outputs will be stored.
#' @param verbose Logical flag controlling function messaging.
#'
#' @keywords internal

BootstrapBoundsSet <- function(ModelType, ModelBootstrap, NumOutPE, InputsForOutputs, Economies, Folder2save, verbose) {
  # Generate the graph paths and the graph folders
  FolderPath <- if (is.null(Folder2save)) tempdir() else Folder2save
  dir.create(paste(FolderPath, "/Outputs/", ModelType, "/Bootstrap", sep = ""))
  if (any(ModelType %in% c("JPS original", "JPS global", "GVAR single"))) {
    for (i in seq_along(Economies)) {
      dir.create(paste(FolderPath, "/Outputs/", ModelType, "/Bootstrap/Model ", Economies[i], sep = ""))
    }
  }

  PathsGraphs <- FolderCreationBoot(ModelType, Economies, FolderPath)

  AllNumOutputs <- list()

  IRFandGIRF <- IRFandGIRFbs(
    ModelType, ModelBootstrap, NumOutPE, InputsForOutputs, Economies, PathsGraphs,
    Folder2save, verbose
  )
  FEVDandGFEVD <- FEVDandGFEVDbs(
    ModelType, ModelBootstrap, NumOutPE, InputsForOutputs, Economies, PathsGraphs,
    Folder2save, verbose
  )

  AllNumOutputs$Factors <- append(IRFandGIRF$Factors, FEVDandGFEVD$Factors)
  AllNumOutputs$Yields <- append(IRFandGIRF$Yields, FEVDandGFEVD$Yields)

  if (verbose) message(paste("Desired graphs are saved in your chosen directory. Please, check:", FolderPath, "\n"))
  return(AllNumOutputs)
}

######################################################################################################
######################################################################################################
#' Creates the confidence bounds and the graphs of IRFs and GIRFs after bootstrap
#'
#' @param ModelType string-vector containing the label of the model to be estimated
#' @param ModelBootstrap list containing the complete set of model parameters after bootstrap estimation procedure
#' @param NumOutPE list of model parameter point estimates
#' @param InputsForOutputs list containing the desired inputs for the construction of the outputs of interest
#' @param Economies string-vector containing the names of the economies which are part of the economic system
#' @param PathsGraphs path of the folder in which the graphs will be saved
#' @param Folder2save Folder path where the outputs will be stored.
#' @param verbose Logical flag controlling function messaging.
#'
#' @keywords internal

IRFandGIRFbs <- function(ModelType, ModelBootstrap, NumOutPE, InputsForOutputs, Economies, PathsGraphs,
                         Folder2save, verbose) {
  ndraws <- InputsForOutputs[[ModelType]]$Bootstrap$ndraws
  pctg <- InputsForOutputs[[ModelType]]$Bootstrap$pctg

  # 1) Define the percentiles
  pctg_inf <- (100 - pctg) / 2
  pctg_sup <- 100 - (100 - pctg) / 2
  quants <- c(pctg_inf, 50, pctg_sup) / 100 # Desired quantiles

  # 2) Define some elements of interest
  J <- length(ModelBootstrap$GeneralInputs$mat)
  C <- length(Economies)
  Horiz <- InputsForOutputs[[ModelType]]$IRF$horiz

  # 3) Compute confidence bounds
  if (ModelType %in% c("JPS original", "JPS global", "GVAR single")) {
    K <- nrow(ModelBootstrap$ParaDraws[[ModelType]][[Economies[1]]][[1]]$ModEst$P$K1Z)
    NumOutBounds <- ComputeBounds_IRFandGIRF(ModelBootstrap, quants, K, J, ModelType, Economies, ndraws, Horiz)
  } else {
    K <- nrow(ModelBootstrap$ParaDraws[[ModelType]][[1]]$ModEst$P$K1Z)
    NumOutBounds <- ComputeBounds_IRFandGIRF(ModelBootstrap, quants, K, C * J, ModelType, Economies, ndraws, Horiz)
  }

  # 4) PLOTS
  WG <- WishGraphs_IRFandGIRF_Boot(InputsForOutputs, ModelType)

  if (any(c(WG$Fac, WG$Yields))) {
    if (verbose) message(" ** IRFs/GIRFs (Bootstrap)")

    # a) Factors
    if (any(WG$Fac)) {
      Boot_Fac_Graphs(
        NumOutBounds, NumOutPE, ModelType, K, Horiz, Economies, PathsGraphs, "IRF", Folder2save,
        WG$Fac, WG$Fac_Ortho
      )
    }

    # b) Yields
    if (any(WG$Yields)) {
      Boot_Yields_Graphs(
        NumOutBounds, NumOutPE, ModelType, K, J, Horiz, Economies, PathsGraphs, "IRF", Folder2save,
        WG$Yields, WG$Yields_Ortho
      )
    }
  }

  return(NumOutBounds)
}

##############################################################################################################
#' Creates the confidence bounds and the graphs of FEVDs and GFEVDs after bootstrap (all models)
#'
#' @param ModelType string-vector containing the label of the model to be estimated
#' @param ModelBootstrap list containing the complete set of model parameters after bootstrap estimation procedure
#' @param NumOutPE list of model parameter point estimates
#' @param InputsForOutputs list containing the desired inputs for the construction of the outputs of interest
#' @param Economies string-vector containing the names of the economies which are part of the economic system
#' @param PathsGraphs path of the folder in which the graphs will be saved
#' @param Folder2save Folder path where the outputs will be stored.
#' @param verbose Logical flag controlling function messaging.
#'
#' @keywords internal

FEVDandGFEVDbs <- function(ModelType, ModelBootstrap, NumOutPE, InputsForOutputs, Economies,
                           PathsGraphs, Folder2save, verbose) {
  ndraws <- InputsForOutputs[[ModelType]]$Bootstrap$ndraws
  pctg <- InputsForOutputs[[ModelType]]$Bootstrap$pctg

  # 1) Define the percentiles
  pctg_inf <- (100 - pctg) / 2
  pctg_sup <- 100 - (100 - pctg) / 2
  quants <- c(pctg_inf, 50, pctg_sup) / 100 # Desired quantiles

  # 2) Define some elements of interest
  J <- length(ModelBootstrap$GeneralInputs$mat)
  C <- length(Economies)
  Horiz <- InputsForOutputs[[ModelType]]$IRF$horiz

  # 3) Compute confidence bounds
  if (ModelType %in% c("JPS original", "JPS global", "GVAR single")) {
    K <- nrow(ModelBootstrap$ParaDraws[[ModelType]][[Economies[1]]][[1]]$ModEst$P$K1Z)
    NumOutBounds <- ComputeBounds_FEVDandGFEVD(ModelBootstrap, quants, K, J, ModelType, Economies, ndraws, Horiz)
  } else {
    K <- nrow(ModelBootstrap$ParaDraws[[ModelType]][[1]]$ModEst$P$K1Z)
    NumOutBounds <- ComputeBounds_FEVDandGFEVD(ModelBootstrap, quants, K, C * J, ModelType, Economies, ndraws, Horiz)
  }

  # 4) PLOTS
  WG <- WishGraphs_FEVDandGFEVD_Boot(InputsForOutputs, ModelType)

  if (any(c(WG$Fac, WG$Yields) == 1)) {
    if (verbose) message(" ** FEVDs/GFEVDs (Bootstrap)")

    # a) Factors
    if (any(WG$Fac)) {
      Boot_Fac_Graphs(
        NumOutBounds, NumOutPE, ModelType, K, Horiz, Economies, PathsGraphs, "FEVD",
        Folder2save, WG$Fac, WG$Fac_Ortho
      )
    }

    # b) Yields
    if (any(WG$Yields)) {
      Boot_Yields_Graphs(
        NumOutBounds, NumOutPE, ModelType, K, J, Horiz, Economies, PathsGraphs, "FEVD",
        Folder2save, WG$Yields, WG$Yields_Ortho
      )
    }
  }

  return(NumOutBounds)
}

##############################################################################################################
#' Compute the confidence bounds from the model's numerical outputs
#'
#' @param ModelBootstrap numerical output set from the bootstrap analysis
#' @param quants quantile of the confidence bounds
#' @param FacDim dimension of the risk factor set
#' @param YieDim dimension of the bond yield set
#' @param ModelType Desired model type
#' @param Economies Economies that are part of the economic system
#' @param ndraws number of draws selected
#' @param Horiz horizon of numerical outputs
#'
#' @keywords internal

ComputeBounds_IRFandGIRF <- function(ModelBootstrap, quants, FacDim, YieDim, ModelType, Economies, ndraws, Horiz) {
  LabIRF <- c("IRF", "GIRF")

  # 1) Factors
  NumOutBounds_Fac <- FactorBounds_IRFandGIRF(
    ModelBootstrap, quants, ModelType, ndraws, Horiz, FacDim, LabIRF,
    Economies
  )
  # 2) Yields
  NumOutBounds_Yie <- YieldBounds_IRFandGIRF(
    ModelBootstrap, quants, ModelType, ndraws, Horiz, FacDim, YieDim,
    LabIRF, Economies
  )
  # Export output
  Out <- list(Factors = NumOutBounds_Fac, Yields = NumOutBounds_Yie)

  return(Out)
}

###########################################################################################################
#' Compute the confidence bounds for the model P-dynamics
#'
#' @param ModelBootstrap numerical output set from the bootstrap analysis
#' @param quants quantile of the confidence bounds
#' @param ModelType desired model type
#' @param ndraws number of draws
#' @param Horiz horizon of numerical outputs
#' @param FacDim dimension of the risk factor set
#' @param LabIRF vector containing the labels "IRF" and "GIRF"
#' @param Economies Economies that are part of the economic system
#'
#' @keywords internal

FactorBounds_IRFandGIRF <- function(ModelBootstrap, quants, ModelType, ndraws, Horiz, FacDim, LabIRF, Economies) {
  NumOutBounds_Fac <- list()

  # 1) For models estimated on a country-by-country basis
  if (ModelType %in% c("JPS original", "JPS global", "GVAR single")) {
    NumOutBounds_Fac <- stats::setNames(
      lapply(LabIRF, function(lab) {
        stats::setNames(
          lapply(Economies, function(econ) {
            DrawSet <- ModelBootstrap$NumOutDraws[[lab]][[ModelType]][[econ]]
            DimLabelsFac <- dimnames(DrawSet[[1]]$Factors)
            FacQuantile_bs(DrawSet, lab, ndraws, quants, Horiz, FacDim, DimLabelsFac, ModelType)
          }),
          Economies
        )
      }),
      LabIRF
    )

    # Adjust lists to export
    NumOutBounds_Fac <- lapply(NumOutBounds_Fac, function(country_list) {
      stats::setNames(list(country_list), ModelType)
    })
  } else {
    # 2) For models estimated for countries jointly
    JLL_ModLabel <- c("JLL original", "JLL No DomUnit", "JLL joint Sigma")

    GetBounds_joint_country <- function(lab) {
      DrawSet <- ModelBootstrap$NumOutDraws[[lab]][[ModelType]]
      DimLabelsFac <- if (ModelType %in% JLL_ModLabel) dimnames(DrawSet[[1]]$Factors$NonOrtho) else dimnames(DrawSet[[1]]$Factors)
      FacQuantile_bs(DrawSet, lab, ndraws, quants, Horiz, FacDim, DimLabelsFac, ModelType)
    }

    # Adjust lists to export
    NumOutBounds_Fac <- stats::setNames(
      lapply(LabIRF, function(lab) {
        result <- GetBounds_joint_country(lab)
        stats::setNames(list(result), ModelType)
      }),
      LabIRF
    )

    # 3) Orthogonalized version (JLL-based models)
    if (ModelType %in% JLL_ModLabel) {
      LabIRF_Ortho <- c("IRF_Ortho", "GIRF_Ortho")

      GetBounds_joint_country_ortho <- function(nn) {
        DrawSet <- ModelBootstrap$NumOutDraws[[LabIRF[nn]]][[ModelType]]
        DimLabelsFac <- dimnames(DrawSet[[1]]$Factors$Ortho)
        FacQuantile_bs(DrawSet, LabIRF_Ortho[nn], ndraws, quants, Horiz, FacDim, DimLabelsFac, ModelType, Ortho = TRUE)
      }

      NumOutBounds_Ortho <- stats::setNames(
        lapply(seq_along(LabIRF_Ortho), GetBounds_joint_country_ortho),
        LabIRF
      )

      # Adjust lists to export
      for (nn in seq_along(LabIRF)) {
        NumOutBounds_Fac[[LabIRF[nn]]][[ModelType]]$Ortho <- NumOutBounds_Ortho[[LabIRF[nn]]]
      }
    }
  }

  return(NumOutBounds_Fac)
}

###############################################################################################
#' Compute quantiles for model P-dynamics
#'
#' @param DrawSet Draw-specific set
#' @param LabIRF vector containing the labels "IRF" and "GIRF"
#' @param ndraws number of draws
#' @param quants quantile of the confidence bounds
#' @param Horiz horizon of numerical outputs
#' @param FacDim dimension of the risk factor set
#' @param DimLabelsFac labels of the factor set
#' @param ModelType desired model type
#' @param Ortho Orthogonolized version for the JLL models. Default is FALSE.
#'
#' @keywords internal

FacQuantile_bs <- function(DrawSet, LabIRF, ndraws, quants, Horiz, FacDim, DimLabelsFac, ModelType, Ortho = FALSE) {
  # Initialization
  INFfacs <- array(NA, c(Horiz, FacDim, FacDim)) # Lower bound
  MEDfacs <- array(NA, c(Horiz, FacDim, FacDim)) # Median
  SUPfacs <- array(NA, c(Horiz, FacDim, FacDim)) # Upper bound

  dimnames(INFfacs) <- DimLabelsFac
  dimnames(MEDfacs) <- DimLabelsFac
  dimnames(SUPfacs) <- DimLabelsFac

  # Extract the correct data based on the model type
  if (ModelType %in% c("JLL original", "JLL No DomUnit", "JLL joint Sigma")) {
    DataArray <- array(NA, dim = c(ndraws, FacDim, Horiz, FacDim))

    for (g in seq_len(ndraws)) {
      DataArray[g, , , ] <- if (Ortho) {
        DrawSet[[g]]$Factors$Ortho
      } else {
        DrawSet[[g]]$Factors$NonOrtho
      }
    }
  } else {
    DataArray <- array(NA, dim = c(ndraws, FacDim, Horiz, FacDim))

    for (g in seq_len(ndraws)) {
      DataArray[g, , , ] <- DrawSet[[g]]$Factors
    }
  }

  # Calculate quantiles without looping through shocks
  INFfacs[] <- apply(DataArray, c(2, 3, 4), stats::quantile, probs = quants[1])
  MEDfacs[] <- apply(DataArray, c(2, 3, 4), stats::quantile, probs = quants[2])
  SUPfacs[] <- apply(DataArray, c(2, 3, 4), stats::quantile, probs = quants[3])

  # Store results
  NumOutBounds_Fac <- list(
    INF = INFfacs,
    MED = MEDfacs,
    SUP = SUPfacs
  )

  return(NumOutBounds_Fac)
}

################################################################################################
#' Compute the confidence bounds for the model bond yield-related outputs
#'
#' @param ModelBootstrap numerical output set from the bootstrap analysis
#' @param quants quantile of the confidence bounds
#' @param ModelType desired model type
#' @param ndraws number of draws
#' @param Horiz horizon of numerical outputs
#' @param FacDim dimension of the risk factor set
#' @param YieDim dimension of the bond yield set
#' @param LabIRF vector containing the labels "IRF" and "GIRF"
#' @param Economies Economies that are part of the economic system
#'
#' @keywords internal

YieldBounds_IRFandGIRF <- function(ModelBootstrap, quants, ModelType, ndraws, Horiz, FacDim, YieDim, LabIRF, Economies) {
  NumOutBounds_Yields <- list()

  # 1) For models estimated on a country-by-country basis
  if (ModelType %in% c("JPS original", "JPS global", "GVAR single")) {
    C <- length(Economies)

    for (nn in 1:length(LabIRF)) {
      for (i in 1:C) {
        DrawSet <- ModelBootstrap$NumOutDraws[[LabIRF[nn]]][[ModelType]][[Economies[i]]]
        DimLabelsYields <- dimnames(DrawSet[[1]]$Yields)

        NumOutBounds_CS <- YieldQuantile_bs(
          DrawSet, LabIRF[nn], ndraws, quants, Horiz, FacDim, YieDim,
          DimLabelsYields, ModelType
        )
        NumOutBounds_Yields[[LabIRF[nn]]][[ModelType]][[Economies[i]]] <- NumOutBounds_CS
      }
    }
  } else {
    # 2) For models estimated for countries jointly
    JLL_ModLabel <- c("JLL original", "JLL No DomUnit", "JLL joint Sigma")

    for (nn in 1:length(LabIRF)) {
      DrawSet <- ModelBootstrap$NumOutDraws[[LabIRF[nn]]][[ModelType]]

      if (ModelType %in% JLL_ModLabel) {
        DimLabelsYields <- dimnames(DrawSet[[1]]$Yields$NonOrtho)
      } else {
        DimLabelsYields <- dimnames(DrawSet[[1]]$Yields)
      }

      NumOutBounds_AllEco <- YieldQuantile_bs(
        DrawSet, LabIRF[nn], ndraws, quants, Horiz, FacDim, YieDim,
        DimLabelsYields, ModelType
      )
      NumOutBounds_Yields[[LabIRF[nn]]][[ModelType]] <- NumOutBounds_AllEco
    }

    # Orthogonalized version (JLL-based models)
    if (ModelType %in% JLL_ModLabel) {
      LabIRF_Ortho <- c("IRF_Ortho", "GIRF_Ortho")

      for (nn in 1:length(LabIRF_Ortho)) {
        DrawSet <- ModelBootstrap$NumOutDraws[[LabIRF[nn]]][[ModelType]]
        DimLabelsYields <- dimnames(DrawSet[[1]]$Yields$Ortho)
        NumOutBounds_AllEco_Ortho <- YieldQuantile_bs(DrawSet, LabIRF_Ortho[nn], ndraws, quants, Horiz, FacDim,
          YieDim, DimLabelsYields, ModelType,
          Ortho = TRUE
        )

        NumOutBounds_Yields[[LabIRF[nn]]][[ModelType]]$Ortho <- NumOutBounds_AllEco_Ortho
      }
    }
  }

  return(NumOutBounds_Yields)
}

###############################################################################################
#' Compute quantiles for model bond yield-related outputs
#'
#' @param DrawSet Draw-specific set
#' @param LabIRF vector containing the labels "IRF" and "GIRF"
#' @param ndraws number of draws
#' @param quants quantile of the confidence bounds
#' @param Horiz horizon of numerical outputs
#' @param FacDim dimension of the risk factor set
#' @param YieDim dimension of the bond yield set
#' @param LabelsYies labels of the factor set
#' @param ModelType desired model type
#' @param Ortho Orthogonolized version for the JLL models. Default is FALSE.
#'
#' @keywords internal

YieldQuantile_bs <- function(DrawSet, LabIRF, ndraws, quants, Horiz, FacDim, YieDim, LabelsYies, ModelType, Ortho = FALSE) {
  # Initialization
  INFyields <- array(NA, c(Horiz, YieDim, FacDim)) # Lower bound
  MEDyields <- array(NA, c(Horiz, YieDim, FacDim)) # Median
  SUPyields <- array(NA, c(Horiz, YieDim, FacDim)) # Upper bound

  dimnames(INFyields) <- LabelsYies
  dimnames(MEDyields) <- LabelsYies
  dimnames(SUPyields) <- LabelsYies

  # Efficiently extract all data at once
  if (ModelType %in% c("JLL original", "JLL No DomUnit", "JLL joint Sigma")) {
    DataArray <- array(NA, dim = c(ndraws, YieDim, Horiz, FacDim))

    for (g in seq_len(ndraws)) {
      DataArray[g, , , ] <- if (Ortho) {
        DrawSet[[g]]$Yields$Ortho
      } else {
        DrawSet[[g]]$Yields$NonOrtho
      }
    }
  } else {
    DataArray <- array(NA, dim = c(ndraws, YieDim, Horiz, FacDim))

    for (g in seq_len(ndraws)) {
      DataArray[g, , , ] <- DrawSet[[g]]$Yields
    }
  }

  # Directly compute quantiles without loops
  INFyields[] <- apply(DataArray, c(2, 3, 4), stats::quantile, probs = quants[1])
  MEDyields[] <- apply(DataArray, c(2, 3, 4), stats::quantile, probs = quants[2])
  SUPyields[] <- apply(DataArray, c(2, 3, 4), stats::quantile, probs = quants[3])

  # Store results
  NumOutBounds_Yie <- list(
    INF = INFyields,
    MED = MEDyields,
    SUP = SUPyields
  )

  return(NumOutBounds_Yie)
}
###################################################################################################
#' Extract graphs of interest (bootstrap version)
#'
#' @param InputsForOutputs list containing the desired inputs for the construction of IRFs, GIRFs, FEVDs, and GFEVDs
#' @param ModelType desired model type
#'
#' @keywords internal

WishGraphs_IRFandGIRF_Boot <- function(InputsForOutputs, ModelType) {
  # Factors
  WishGraphFac <- c(
    InputsForOutputs[[ModelType]]$IRF$WishGraphs$RiskFactorsBootstrap,
    InputsForOutputs[[ModelType]]$GIRF$WishGraphs$RiskFactorsBootstrap
  )

  if (ModelType %in% c("JLL original", "JLL No DomUnit", "JLL joint Sigma")) {
    WishGraphFac_Ortho <- c(
      InputsForOutputs[[ModelType]]$IRF$WishGraphsOrtho$RiskFactorsBootstrap,
      InputsForOutputs[[ModelType]]$GIRF$WishGraphsOrtho$RiskFactorsBootstrap
    )
  } else {
    WishGraphFac_Ortho <- NULL
  }


  # Yields
  WishGraphYields <- c(
    InputsForOutputs[[ModelType]]$IRF$WishGraphs$YieldsBootstrap,
    InputsForOutputs[[ModelType]]$GIRF$WishGraphs$YieldsBootstrap
  )

  if (ModelType %in% c("JLL original", "JLL No DomUnit", "JLL joint Sigma")) {
    WishGraphiYields_Ortho <- c(
      InputsForOutputs[[ModelType]]$IRF$WishGraphsOrtho$YieldsBootstrap,
      InputsForOutputs[[ModelType]]$GIRF$WishGraphsOrtho$YieldsBootstrap
    )
  } else {
    WishGraphiYields_Ortho <- NULL
  }


  Out <- list(
    Fac = WishGraphFac, Fac_Ortho = WishGraphFac_Ortho, Yields = WishGraphYields,
    Yields_Ortho = WishGraphiYields_Ortho
  )

  return(Out)
}
####################################################################################################
#' Build P-dynamic graphs after the bootstrap implementation
#'
#' @param NumOutBounds numerical output set from the bootstrap analysis
#' @param NumOutPE numerical output set from the point estimate analysis
#' @param ModelType desired model type
#' @param FacDim dimension of the risk factor set
#' @param Horiz horizon of numerical outputs
#' @param Economies Economies that are part of the economic system
#' @param PathsGraphs Path to save the desired graphs
#' @param OutInt available options are "IRF", "FEVD", "GIRF" or "GFEVD"
#' @param Folder2save Folder path where the outputs will be stored.
#' @param WishFacGraphs Binary variable reflecting the graphs of interest
#' @param WishFacGraphsOrtho Binary variable reflecting the graphs of interest (orthogonalized version). Default is NULL
#'
#' @keywords internal

Boot_Fac_Graphs <- function(NumOutBounds, NumOutPE, ModelType, FacDim, Horiz, Economies, PathsGraphs, OutInt,
                            Folder2save, WishFacGraphs, WishFacGraphsOrtho = NULL) {
  C <- length(Economies)

  if (OutInt == "IRF") {
    Lab_Int <- c("IRF", "GIRF")
    Graph_Lab <- "Factors_shock_to_"
  } else {
    Lab_Int <- c("FEVD", "GFEVD")
    Graph_Lab <- "Factors_"
  }

  JLL_ModLabels <- c("JLL original", "JLL No DomUnit", "JLL joint Sigma")

  IdxFacGraphs <- which(WishFacGraphs == 1)

  Autoplot_List <- list()

  ################ 1) Estimation done for countries individually ################
  if (any(ModelType %in% c("JPS original", "JPS global", "GVAR single"))) {
    for (f in 1:C) {
      for (d in IdxFacGraphs) {
        # Adjust the graph path
        if (!is.null(Folder2save)) {
          if (OutInt == "IRF") {
            PathAdj <- AdjustPathIRFs(Lab_Int[d], "Factors", PathsGraphs, Economies[f], ModelType)
          } else {
            PathAdj <- AdjustPathFEVDs(Lab_Int[d], "Factors", PathsGraphs, Economies[f], ModelType)
          }
        }

        # Labels of shocks and response variables
        nmResponse <- dimnames(NumOutPE[[Lab_Int[d]]][[ModelType]][[Economies[f]]]$Factors)[[2]]
        nmShock <- dimnames(NumOutPE[[Lab_Int[d]]][[ModelType]][[Economies[f]]]$Factors)[[3]]

        # Folder Creation
        if (!is.null(Folder2save)) {
          FolderCreation_Boot(ModelType, Lab_Int[d], Economies[f], "Factors", Folder2save)
        }

        # Create plots
        for (b in 1:FacDim) {
          plotlist_OneShock <- Boot_DataGraphFact_perShock(
            NumOutBounds, NumOutPE, b, nmResponse, Lab_Int[d],
            ModelType, FacDim, Horiz, Economies[f]
          )

          subplots <- plot_grid(plotlist = plotlist_OneShock, ncol = 3)
          if (!is.null(Folder2save)) {
            suppressMessages(ggplot2::ggsave(subplots,
              file = paste0(Lab_Int[d], Graph_Lab, nmShock[b], ".png"),
              path = PathAdj
            ))
            print(subplots)
          }
          Autoplot_List[[Economies[f]]][[nmShock[b]]] <- subplots
        }
      }
    }
  } else {
    ################ 2) Estimation done for countries jointly ###############################
    for (d in IdxFacGraphs) {
      # Adjust the graph path
      if (!is.null(Folder2save)) {
        if (OutInt == "IRF") {
          PathAdj <- AdjustPathIRFs(Lab_Int[d], "Factors", PathsGraphs, Economies, ModelType)
        } else {
          PathAdj <- AdjustPathFEVDs(Lab_Int[d], "Factors", PathsGraphs, Economies, ModelType)
        }
      }

      # Labels of shocks and response variables
      if (ModelType %in% JLL_ModLabels) {
        nmResponse <- dimnames(NumOutPE[[Lab_Int[d]]][[ModelType]]$Factors$NonOrtho)[[2]]
        nmShock <- dimnames(NumOutPE[[Lab_Int[d]]][[ModelType]]$Factors$NonOrtho)[[3]]
      } else {
        nmResponse <- dimnames(NumOutPE[[Lab_Int[d]]][[ModelType]]$Factors)[[2]]
        nmShock <- dimnames(NumOutPE[[Lab_Int[d]]][[ModelType]]$Factors)[[3]]
      }
      # Folder Creation
      if (!is.null(Folder2save)) {
        FolderCreation_Boot(ModelType, Lab_Int[d], Economies, "Factors", Folder2save)
      }

      # Create plots
      for (b in 1:FacDim) {
        plotlist_OneShock <- Boot_DataGraphFact_perShock(
          NumOutBounds, NumOutPE, b, nmResponse, Lab_Int[d],
          ModelType, FacDim, Horiz
        )

        subplots <- plot_grid(plotlist = plotlist_OneShock, ncol = 3)
        if (!is.null(Folder2save)) {
          suppressMessages(ggplot2::ggsave(subplots,
            file = paste0(Lab_Int[d], Graph_Lab, nmShock[b], ".png"),
            path = PathAdj
          ))
          print(subplots)
        }
        Autoplot_List[[nmShock[b]]] <- subplots
      }
    }
  }

  ################ 3) Orthogonalized version for JLL-based models ################################
  if (ModelType %in% JLL_ModLabels) {
    if (OutInt == "IRF") {
      Lab_Int_Ortho <- c("IRF Ortho", "GIRF Ortho")
    } else if (OutInt == "FEVD") {
      Lab_Int_Ortho <- c("FEVD Ortho", "GFEVD Ortho")
    }

    IdxFacGraphs_Ortho <- which(WishFacGraphsOrtho == 1)

    for (d in IdxFacGraphs_Ortho) {
      # Adjust the graph path
      if (!is.null(Folder2save)) {
        if (OutInt == "IRF") {
          PathAdj <- AdjustPathIRFs(Lab_Int_Ortho[d], "Factors", PathsGraphs, Economies, ModelType)
        } else {
          PathAdj <- AdjustPathFEVDs(Lab_Int_Ortho[d], "Factors", PathsGraphs, Economies, ModelType)
        }
      }
      # Labels of shocks and response variables
      nmResponse <- dimnames(NumOutPE[[Lab_Int[d]]][[ModelType]]$Factors$Ortho)[[2]]
      nmShock <- dimnames(NumOutPE[[Lab_Int[d]]][[ModelType]]$Factors$Ortho)[[3]]

      # Folder Creation
      if (!is.null(Folder2save)) {
        FolderCreation_Boot(ModelType, Lab_Int[d], Economies, "Factors", Folder2save, Ortho = TRUE)
      }

      # Create plots
      for (b in 1:FacDim) {
        plotlist_OneShock <- Boot_DataGraphFact_perShock(NumOutBounds, NumOutPE, b, nmResponse, Lab_Int[d],
          ModelType, FacDim, Horiz,
          Ortho = TRUE
        )

        subplots <- plot_grid(plotlist = plotlist_OneShock, ncol = 3)
        if (!is.null(Folder2save)) {
          suppressMessages(ggplot2::ggsave(subplots,
            file = paste0(Lab_Int[d], Graph_Lab, nmShock[b], "_ORTHO", ".png"),
            path = PathAdj
          ))
          print(subplots)
        }
        Autoplot_List[[nmShock[b]]] <- subplots
      }
    }
  }

  return(Autoplot_List)
}

###############################################################################################
#' Creates folder to store graphs generated from the bootstrap analysis
#'
#' @param ModelType Desired model type
#' @param LabIRF Output types "IRF", "GIRF" and "IRF Ortho"
#' @param Economies economies of the economic system
#' @param OutType Available option "Factors" or "Yields
#' @param Folder2save Folder path where the outputs will be stored.
#' @param Ortho Option for orthogonal outputs, for JLL models. Default is FALSE.
#'
#' @keywords internal

FolderCreation_Boot <- function(ModelType, LabIRF, Economies, OutType, Folder2save, Ortho = FALSE) {
  # A) Folder for factors
  if (OutType == "Factors") {
    # 1) Models estimated on a country-by-country basis
    if (any(ModelType %in% c("JPS original", "JPS global", "GVAR single"))) {
      dir.create(paste(Folder2save, "/Outputs/", ModelType, "/Bootstrap/Model ", Economies, "/", LabIRF, sep = ""), showWarnings = FALSE)
      dir.create(paste(Folder2save, "/Outputs/", ModelType, "/Bootstrap/Model ", Economies, "/", LabIRF, "/Factors", sep = ""))
    } else {
      # 2) Country estimated jointly
      dir.create(paste(Folder2save, "/Outputs/", ModelType, "/Bootstrap/", LabIRF, sep = ""), showWarnings = FALSE)
      dir.create(paste(Folder2save, "/Outputs/", ModelType, "/Bootstrap/", LabIRF, "/Factors", sep = ""))
      # 3) Orthogonalized outputs (for JLL models)
      if (Ortho) {
        dir.create(paste(Folder2save, "/Outputs/", ModelType, "/Bootstrap/", LabIRF, "/Factors/Ortho", sep = ""))
      }
    }
  } else {
    # B) Folder for yields
    # 1) Models estimated on a country-by-country basis
    if (any(ModelType %in% c("JPS original", "JPS global", "GVAR single"))) {
      dir.create(paste(Folder2save, "/Outputs/", ModelType, "/Bootstrap/Model ", Economies, "/", LabIRF, sep = ""), showWarnings = FALSE)
      dir.create(paste(Folder2save, "/Outputs/", ModelType, "/Bootstrap/Model ", Economies, "/", LabIRF, "/Yields", sep = ""))
    } else {
      # 2) Country estimated jointly
      dir.create(paste(Folder2save, "/Outputs/", ModelType, "/Bootstrap/", LabIRF, sep = ""), showWarnings = FALSE)
      dir.create(paste(Folder2save, "/Outputs/", ModelType, "/Bootstrap/", LabIRF, "/Yields", sep = ""))
      # 3) Orthogonalized outputs (for JLL models)
      if (Ortho) {
        dir.create(paste(Folder2save, "/Outputs/", ModelType, "/Bootstrap/", LabIRF, "/Yields/Ortho", sep = ""))
      }
    }
  }
}

####################################################################################################
#' Generates the desired bootstrap graphs
#'
#' @param NumOutBounds numerical output set from the bootstrap analysis
#' @param NumOutPE numerical output set from the point estimate analysis
#' @param IdxShock index associated with the shock variable
#' @param nmResponse Label of the response variable
#' @param Lab_Int Output types "IRF", "GIRF" and "IRF Ortho"
#' @param ModelType desired model type
#' @param FacDim dimension from the P-dynamics
#' @param Horiz horizon of analysis
#' @param Economies name of economies forming the economic system
#' @param Ortho Option for orthogonal outputs, for JLL models. Default is FALSE.
#'
#' @keywords internal

Boot_DataGraphFact_perShock <- function(NumOutBounds, NumOutPE, IdxShock, nmResponse, Lab_Int, ModelType,
                                        FacDim, Horiz, Economies = NULL, Ortho = FALSE) {
  plot_list <- list()

  for (i in 1:FacDim) {
    # Confidence Bounds
    if (any(ModelType %in% c("JPS original", "JPS global", "GVAR single"))) {
      Paras <- NumOutBounds$Factors[[Lab_Int]][[ModelType]][[Economies]]
      LL <- Paras$INF[, i, IdxShock]
      UU <- Paras$SUP[, i, IdxShock]
      MM <- Paras$MED[, i, IdxShock]
    } else if (any(ModelType %in% c("JLL original", "JLL No DomUnit", "JLL joint Sigma")) & Ortho == 1) {
      Paras_Ortho <- NumOutBounds$Factors[[Lab_Int]][[ModelType]]$Ortho
      LL <- Paras_Ortho$INF[, i, IdxShock]
      UU <- Paras_Ortho$SUP[, i, IdxShock]
      MM <- Paras_Ortho$MED[, i, IdxShock]
    } else {
      Paras <- NumOutBounds$Factors[[Lab_Int]][[ModelType]]
      LL <- Paras$INF[, i, IdxShock]
      UU <- Paras$SUP[, i, IdxShock]
      MM <- Paras$MED[, i, IdxShock]
    }

    # Point estimate
    if (any(ModelType %in% c("JPS original", "JPS global", "GVAR single"))) {
      PP <- NumOutPE[[Lab_Int]][[ModelType]][[Economies]]$Factors[, i, IdxShock]
    } else if (any(ModelType %in% c("JLL original", "JLL No DomUnit", "JLL joint Sigma"))) {
      Para_PE <- NumOutPE[[Lab_Int]][[ModelType]]$Factors
      if (Ortho) {
        PP <- Para_PE$Ortho[, i, IdxShock]
      } else {
        PP <- Para_PE$NonOrtho[, i, IdxShock]
      }
    } else {
      PP <- NumOutPE[[Lab_Int]][[ModelType]]$Factors[, i, IdxShock]
    }

    # Add time-span
    ALLdata <- data.frame(cbind(LL, MM, PP, UU))
    TimeSpan <- 1:Horiz
    ALLdata$TimeSpan <- TimeSpan

    p <- Boot_graph_template(ALLdata, nmResponse[i])

    plot_list[[i]] <- p
  }

  return(plot_list)
}

####################################################################################################
#' Builds template from bootstrap-related graphs
#'
#' @param ALLdata data-frame containing the necessary data for building the graphs
#' @param nmResponse string containing the name of the response variable
#'
#' @keywords internal

Boot_graph_template <- function(ALLdata, nmResponse) {
  LL <- ALLdata$LL
  MM <- ALLdata$MM
  UU <- ALLdata$UU
  PP <- ALLdata$PP
  TimeSpan <- ALLdata$TimeSpan

  p <- ggplot(data = ALLdata, aes(x = TimeSpan)) +
    geom_line(aes(y = LL), color = "gray50", linetype = "dashed") +
    geom_line(aes(y = MM), color = "#009E73", linewidth = 1) +
    geom_line(aes(y = UU), color = "gray50", linetype = "dashed") +
    geom_line(aes(y = PP), linewidth = 1) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 6, face = "bold", hjust = 0.5),
      axis.title.x = element_blank(), axis.title.y = element_blank()
    ) +
    ggtitle(nmResponse) +
    geom_hline(yintercept = 0)

  return(p)
}

######################################################################################################
#' Build P-dynamic graphs after the bootstrap implementation
#'
#' @param NumOutBounds numerical output set from the bootstrap analysis
#' @param NumOutPE numerical output set from the point estimate analysis
#' @param ModelType desired model type
#' @param FacDim dimension of the risk factor set
#' @param Horiz horizon of numerical outputs
#' @param Economies Economies that are part of the economic system
#' @param PathsGraphs Path to save the desired graphs
#' @param OutInt Available option are "IRF", "FEVD", "GIRF" or "GFEVD"
#' @param Folder2save Folder path where the outputs will be stored.
#' @param WishYieldGraphs Binary variable reflecting the graphs of interest
#' @param WishYieldGraphsOrtho Binary variable reflecting the graphs of interest (orthogonalized version). Default is NULL
#'
#' @keywords internal

Boot_Yields_Graphs <- function(NumOutBounds, NumOutPE, ModelType, FacDim, YielDim, Horiz, Economies, PathsGraphs,
                               OutInt, Folder2save, WishYieldGraphs, WishYieldGraphsOrtho = NULL) {
  C <- length(Economies)
  if (OutInt == "IRF") {
    Lab_Int <- c("IRF", "GIRF")
    Graph_Lab <- "Yields_shock_to_"
  } else {
    Lab_Int <- c("FEVD", "GFEVD")
    Graph_Lab <- "Yields_"
  }
  JLL_ModLabels <- c("JLL original", "JLL No DomUnit", "JLL joint Sigma")

  IdxYielddGraphs <- which(WishYieldGraphs == 1)

  Autoplot_List <- list()

  ################ 1) Estimation done for countries individually ################
  if (ModelType %in% c("JPS original", "JPS global", "GVAR single")) {
    for (f in 1:C) {
      for (d in IdxYielddGraphs) {
        # Adjust the graph path
        if (!is.null(Folder2save)) {
          if (OutInt == "IRF") {
            PathAdj <- AdjustPathIRFs(Lab_Int[d], "Yields", PathsGraphs, Economies[f], ModelType)
          } else {
            PathAdj <- AdjustPathFEVDs(Lab_Int[d], "Yields", PathsGraphs, Economies[f], ModelType)
          }
        }

        # Labels of shocks and response variables
        nmResponse <- dimnames(NumOutPE[[Lab_Int[d]]][[ModelType]][[Economies[f]]]$Yields)[[2]]
        nmShock <- dimnames(NumOutPE[[Lab_Int[d]]][[ModelType]][[Economies[f]]]$Yields)[[3]]

        # Folder Creation
        if (!is.null(Folder2save)) {
          FolderCreation_Boot(ModelType, Lab_Int[d], Economies[f], "Yields", Folder2save)
        }

        # Create plots
        if (Lab_Int[d] %in% c("IRF", "GIRF")) {
          DimInt <- FacDim
        } else if (Lab_Int[d] %in% c("FEVD", "GFEVD")) {
          DimInt <- YielDim
        }

        for (b in 1:DimInt) {
          plot_list <- Boot_DataGraphYield_perShock(
            NumOutBounds, NumOutPE, b, nmResponse, Lab_Int[d],
            ModelType, FacDim, YielDim, Horiz, Economies[f]
          )

          subplots <- plot_grid(plotlist = plot_list, ncol = 3)
          if (!is.null(Folder2save)) {
            suppressMessages(ggplot2::ggsave(subplots,
              file = paste0(Lab_Int[d], Graph_Lab, nmShock[b], ".png"),
              path = PathAdj
            ))
            print(subplots)
          }

          Autoplot_List[[Economies[f]]][[nmShock[b]]] <- subplots
        }
      }
    }
  } else {
    ################ 2) Estimation done for countries jointly ###############################
    for (d in IdxYielddGraphs) {
      # Adjust the graph path
      if (!is.null(Folder2save)) {
        if (OutInt == "IRF") {
          PathAdj <- AdjustPathIRFs(Lab_Int[d], "Yields", PathsGraphs, Economies, ModelType)
        } else {
          PathAdj <- AdjustPathFEVDs(Lab_Int[d], "Yields", PathsGraphs, Economies, ModelType)
        }
      }

      # Labels of shocks and response variables
      if (ModelType %in% JLL_ModLabels) {
        nmResponse <- dimnames(NumOutPE[[Lab_Int[d]]][[ModelType]]$Yields$NonOrtho)[[2]]
        nmShock <- dimnames(NumOutPE[[Lab_Int[d]]][[ModelType]]$Yields$NonOrtho)[[3]]
      } else {
        nmResponse <- dimnames(NumOutPE[[Lab_Int[d]]][[ModelType]]$Yields)[[2]]
        nmShock <- dimnames(NumOutPE[[Lab_Int[d]]][[ModelType]]$Yields)[[3]]
      }

      # Folder Creation
      if (!is.null(Folder2save)) {
        FolderCreation_Boot(ModelType, Lab_Int[d], Economies, "Yields", Folder2save)
      }

      # Create plots
      if (Lab_Int[d] %in% c("IRF", "GIRF")) {
        DimInt <- FacDim
      } else if (Lab_Int[d] %in% c("FEVD", "GFEVD")) {
        DimInt <- C * YielDim
      }

      for (b in 1:DimInt) {
        plot_list <- Boot_DataGraphYield_perShock(
          NumOutBounds, NumOutPE, b, nmResponse, Lab_Int[d],
          ModelType, FacDim, C * YielDim, Horiz, Economies
        )

        subplots <- plot_grid(plotlist = plot_list, ncol = 3)
        if (!is.null(Folder2save)) {
          suppressMessages(ggplot2::ggsave(subplots,
            file = paste0(Lab_Int[d], Graph_Lab, nmShock[b], ".png"),
            path = PathAdj
          ))
          print(subplots)
        }

        Autoplot_List[[nmShock[b]]] <- subplots
      }
    }
  }

  ################ 3) Orthogonalized version for JLL-based models ################################
  if (ModelType %in% JLL_ModLabels) {
    if (OutInt == "IRF") {
      Lab_Int_Ortho <- c("IRF Ortho", "GIRF Ortho")
    } else {
      Lab_Int_Ortho <- c("FEVD Ortho", "GFEVD Ortho")
    }
    IdxYieldGraphs_Ortho <- which(WishYieldGraphsOrtho == 1)

    for (d in IdxYieldGraphs_Ortho) {
      # Adjust the graph path
      if (!is.null(Folder2save)) {
        if (OutInt == "IRF") {
          PathAdj <- AdjustPathIRFs(Lab_Int_Ortho[d], "Yields", PathsGraphs, Economies, ModelType)
        } else {
          PathAdj <- AdjustPathFEVDs(Lab_Int_Ortho[d], "Yields", PathsGraphs, Economies, ModelType)
        }
      }

      # Labels of shocks and response variables
      nmResponse <- dimnames(NumOutPE[[Lab_Int[d]]][[ModelType]]$Yields$Ortho)[[2]]
      nmShock <- dimnames(NumOutPE[[Lab_Int[d]]][[ModelType]]$Yields$Ortho)[[3]]

      # Folder Creation
      if (!is.null(Folder2save)) {
        FolderCreation_Boot(ModelType, Lab_Int[d], Economies, "Yields", Folder2save, Ortho = TRUE)
      }

      # Create plots
      if (Lab_Int[d] %in% c("IRF", "GIRF")) {
        DimInt <- FacDim
      } else if (Lab_Int[d] %in% c("FEVD", "GFEVD")) {
        DimInt <- C * YielDim
      }

      for (b in 1:DimInt) {
        plot_list <- Boot_DataGraphYield_perShock(NumOutBounds, NumOutPE, b, nmResponse, Lab_Int[d],
          ModelType, FacDim, C * YielDim, Horiz, Economies,
          Ortho = TRUE
        )

        subplots <- plot_grid(plotlist = plot_list, ncol = 3)
        if (!is.null(Folder2save)) {
          suppressMessages(ggplot2::ggsave(subplots,
            file = paste0(Lab_Int[d], Graph_Lab, nmShock[b], ".png"),
            path = PathAdj
          ))
          print(subplots)
        }
        Autoplot_List[[nmShock[b]]] <- subplots
      }
    }
  }
  if (is.null(Folder2save)) {
    return(Autoplot_List)
  }
}

####################################################################################################
#' Generates the desired bootstrap graphs
#'
#' @param NumOutBounds numerical output set from the bootstrap analysis
#' @param NumOutPE numerical output set from the point estimate analysis
#' @param IdxShock index associated with the shock variable
#' @param nmResponse Label of the response variable
#' @param Lab_Int Output types "IRF" or "FEVD"
#' @param ModelType desired model type
#' @param FacDim dimension from bond yield set
#' @param YieldDim dimension from the P-dynamics
#' @param Horiz horizon of analysis
#' @param Economies name of economies forming the economic system
#' @param Ortho Option for orthogonal outputs, for JLL models. Default is FALSE.
#'
#' @keywords internal

Boot_DataGraphYield_perShock <- function(
    NumOutBounds, NumOutPE, IdxShock, nmResponse, Lab_Int, ModelType,
    FacDim, YieldDim, Horiz, Economies = NULL, Ortho = FALSE) {
  if (Lab_Int %in% c("IRF", "GIRF")) {
    DimInt <- YieldDim
  } else if (Lab_Int %in% c("FEVD", "GFEVD")) {
    DimInt <- FacDim
  }

  plot_list <- list()

  for (i in 1:DimInt) {
    # Get Confidence interval set
    CI_data <- BuildCI_Yields(NumOutBounds, NumOutPE, Lab_Int, ModelType, Economies, i, IdxShock, Ortho)

    # Build data-frame and add time-span
    TimeSpan <- 1:Horiz
    ALLdata <- data.frame(LL = CI_data$LL, MM = CI_data$MM, PP = CI_data$PP, UU = CI_data$UU, TimeSpan = TimeSpan)

    p <- Boot_graph_template(ALLdata, nmResponse[i])

    plot_list[[i]] <- p
  }

  return(plot_list)
}
###########################################################################################
#' Build Confidence intervals for yield-related outputs
#'
#' @param NumOutBounds numerical output set from the bootstrap analysis
#' @param NumOutPE numerical output set from the point estimate analysis
#' @param Lab_Int Label of interest. available options are "IRF" and "FEVD"
#' @param ModelType desired model type
#' @param Economies name of the economies forming the economic system
#' @param IdxResp index associated with the response variable
#' @param IdxShock index associated with the shock variable
#' @param Ortho Option for orthogonal outputs, for JLL models. Default is FALSE.
#'
#' @keywords internal

BuildCI_Yields <- function(NumOutBounds, NumOutPE, Lab_Int, ModelType, Economies, IdxResp, IdxShock, Ortho = FALSE) {
  # a) Confidence Bounds
  if (ModelType %in% c("JPS original", "JPS global", "GVAR single")) {
    Paras <- NumOutBounds$Yields[[Lab_Int]][[ModelType]][[Economies]]
    LL <- Paras$INF[, IdxResp, IdxShock]
    UU <- Paras$SUP[, IdxResp, IdxShock]
    MM <- Paras$MED[, IdxResp, IdxShock]
  } else if ((ModelType %in% c("JLL original", "JLL No DomUnit", "JLL joint Sigma")) & Ortho == 1) {
    Paras_Ortho <- NumOutBounds$Yields[[Lab_Int]][[ModelType]]$Ortho
    LL <- Paras_Ortho$INF[, IdxResp, IdxShock]
    UU <- Paras_Ortho$SUP[, IdxResp, IdxShock]
    MM <- Paras_Ortho$MED[, IdxResp, IdxShock]
  } else {
    Paras <- NumOutBounds$Yields[[Lab_Int]][[ModelType]]
    LL <- Paras$INF[, IdxResp, IdxShock]
    UU <- Paras$SUP[, IdxResp, IdxShock]
    MM <- Paras$MED[, IdxResp, IdxShock]
  }

  # b) Point estimate
  if (ModelType %in% c("JPS original", "JPS global", "GVAR single")) {
    PP <- NumOutPE[[Lab_Int]][[ModelType]][[Economies]]$Yields[, IdxResp, IdxShock]
  } else if (ModelType %in% c("JLL original", "JLL No DomUnit", "JLL joint Sigma")) {
    Para_PE <- NumOutPE[[Lab_Int]][[ModelType]]$Yields
    if (Ortho) {
      PP <- Para_PE$Ortho[, IdxResp, IdxShock]
    } else {
      PP <- Para_PE$NonOrtho[, IdxResp, IdxShock]
    }
  } else {
    PP <- NumOutPE[[Lab_Int]][[ModelType]]$Yields[, IdxResp, IdxShock]
  }

  return(list(LL = LL, MM = MM, PP = PP, UU = UU))
}

###############################################################################################
#' Compute the confidence bounds around the P-dynamics and bond yields for FEVD and GFEVD
#'
#' @param ModelBootstrap numerical output set from the bootstrap analysis
#' @param quants quantile of the confidence bounds
#' @param FacDim dimension of the risk factor set
#' @param YieDim dimension of the bond yield set
#' @param ModelType desired model type
#' @param Economies Economies that are part of the economic system
#' @param ndraws number of draws
#' @param Horiz horizon of numerical outputs
#'
#' @keywords internal

ComputeBounds_FEVDandGFEVD <- function(ModelBootstrap, quants, FacDim, YieDim, ModelType, Economies, ndraws, Horiz) {
  LabFEVD <- c("FEVD", "GFEVD")

  # 1) Factors
  NumOutBounds_Fac <- FactorBounds_FEVDandGFEVD(ModelBootstrap, quants, ModelType, ndraws, Horiz, FacDim, LabFEVD, Economies)

  # 2) Yields
  # NOTE: in order to avoid over-complicating the code, in the function below the arguments "YieDim" and "FacDim"
  # are swapped wrt to the original function design. This helps to accommodate the code for the different output
  # structures in both IRF and FEVD.
  NumOutBounds_Yie <- YieldBounds_FEVDandGFEVD(
    ModelBootstrap, quants, ModelType, ndraws, Horiz, YieDim,
    FacDim, LabFEVD, Economies
  )
  # Export output
  Out <- list(Factors = NumOutBounds_Fac, Yields = NumOutBounds_Yie)

  return(Out)
}

###############################################################################################
#' Compute the confidence bounds for the model bond P-dynamics-related outputs
#'
#' @param ModelBootstrap numerical output set from the bootstrap analysis
#' @param quants quantile of the confidence bounds
#' @param ModelType desired model type
#' @param ndraws number of draws
#' @param Horiz horizon of numerical outputs
#' @param FacDim dimension of the risk factor set
#' @param LabFEVD vector containing the labels "FEVD" and "GFEVD"
#' @param Economies Economies that are part of the economic system
#'
#' @keywords internal

FactorBounds_FEVDandGFEVD <- function(ModelBootstrap, quants, ModelType, ndraws, Horiz, FacDim, LabFEVD, Economies) {
  NumOutBounds_Fac <- list()

  # 1) For models estimated on a country-by-country basis
  if (ModelType %in% c("JPS original", "JPS global", "GVAR single")) {
    C <- length(Economies)

    for (nn in 1:length(LabFEVD)) {
      for (i in 1:C) {
        DrawSet <- ModelBootstrap$NumOutDraws[[LabFEVD[nn]]][[ModelType]][[Economies[i]]]
        DimLabelsFac <- dimnames(DrawSet[[1]]$Factors)

        NumOutBounds_CS <- FacQuantile_bs(DrawSet, LabFEVD[nn], ndraws, quants, Horiz, FacDim, DimLabelsFac, ModelType)
        NumOutBounds_Fac[[LabFEVD[nn]]][[ModelType]][[Economies[i]]] <- NumOutBounds_CS
      }
    }
  } else {
    # 2) For models estimated for countries jointly
    JLL_ModLabel <- c("JLL original", "JLL No DomUnit", "JLL joint Sigma")

    for (nn in 1:length(LabFEVD)) {
      DrawSet <- ModelBootstrap$NumOutDraws[[LabFEVD[nn]]][[ModelType]]

      if (ModelType %in% JLL_ModLabel) {
        DimLabelsFac <- dimnames(DrawSet[[1]]$Factors$NonOrtho)
      } else {
        DimLabelsFac <- dimnames(DrawSet[[1]]$Factors)
      }

      NumOutBounds_AllEco <- FacQuantile_bs(DrawSet, LabFEVD[nn], ndraws, quants, Horiz, FacDim, DimLabelsFac, ModelType)
      NumOutBounds_Fac[[LabFEVD[nn]]][[ModelType]] <- NumOutBounds_AllEco
    }

    # Orthogonalized version (JLL-based models)
    if (ModelType %in% JLL_ModLabel) {
      LabFEVD_Ortho <- c("FEVD_Ortho", "GFEVD_Ortho")

      for (nn in 1:length(LabFEVD_Ortho)) {
        DrawSet <- ModelBootstrap$NumOutDraws[[LabFEVD[nn]]][[ModelType]]
        DimLabelsFac <- dimnames(DrawSet[[1]]$Factors$Ortho)
        NumOutBounds_AllEco_Ortho <- FacQuantile_bs(DrawSet, LabFEVD_Ortho[nn], ndraws, quants, Horiz, FacDim, DimLabelsFac, ModelType, Ortho = TRUE)

        NumOutBounds_Fac[[LabFEVD[nn]]][[ModelType]]$Ortho <- NumOutBounds_AllEco_Ortho
      }
    }
  }
  return(NumOutBounds_Fac)
}
######################################################################################################
#' Compute the confidence bounds for the model bond yield-related outputs
#'
#' @param ModelBootstrap numerical output set from the bootstrap analysis
#' @param quants quantile of the confidence bounds
#' @param ModelType desired model type
#' @param ndraws number of draws
#' @param Horiz horizon of numerical outputs
#' @param FacDim dimension of the risk factor set
#' @param YieDim dimension of the bond yield set
#' @param LabFEVD vector containing the labels "FEVD" and "GFEVD"
#' @param Economies Economies that are part of the economic system
#'
#' @keywords internal

YieldBounds_FEVDandGFEVD <- function(ModelBootstrap, quants, ModelType, ndraws, Horiz, FacDim, YieDim, LabFEVD,
                                     Economies) {
  NumOutBounds_Yields <- list()

  # 1) For models estimated on a country-by-country basis
  if (ModelType %in% c("JPS original", "JPS global", "GVAR single")) {
    C <- length(Economies)

    for (nn in 1:length(LabFEVD)) {
      for (i in 1:C) {
        DrawSet <- ModelBootstrap$NumOutDraws[[LabFEVD[nn]]][[ModelType]][[Economies[i]]]
        DimLabelsYields <- dimnames(DrawSet[[1]]$Yields)

        NumOutBounds_CS <- YieldQuantile_bs(
          DrawSet, LabFEVD[nn], ndraws, quants, Horiz, FacDim, YieDim,
          DimLabelsYields, ModelType
        )

        NumOutBounds_Yields[[LabFEVD[nn]]][[ModelType]][[Economies[i]]] <- NumOutBounds_CS
      }
    }
  } else {
    # 2) For models estimated for countries jointly
    JLL_ModLabel <- c("JLL original", "JLL No DomUnit", "JLL joint Sigma")

    for (nn in 1:length(LabFEVD)) {
      DrawSet <- ModelBootstrap$NumOutDraws[[LabFEVD[nn]]][[ModelType]]

      if (ModelType %in% JLL_ModLabel) {
        DimLabelsYields <- dimnames(DrawSet[[1]]$Yields$NonOrtho)
      } else {
        DimLabelsYields <- dimnames(DrawSet[[1]]$Yields)
      }

      NumOutBounds_AllEco <- YieldQuantile_bs(DrawSet, LabFEVD[nn], ndraws, quants, Horiz, FacDim, YieDim, DimLabelsYields, ModelType)
      NumOutBounds_Yields[[LabFEVD[nn]]][[ModelType]] <- NumOutBounds_AllEco
    }

    # Orthogonalized version (JLL-based models)
    if (ModelType %in% JLL_ModLabel) {
      LabFEVD_Ortho <- c("FEVD_Ortho", "GFEVD_Ortho")

      for (nn in 1:length(LabFEVD_Ortho)) {
        DrawSet <- ModelBootstrap$NumOutDraws[[LabFEVD[nn]]][[ModelType]]
        DimLabelsYields <- dimnames(DrawSet[[1]]$Yields$Ortho)
        NumOutBounds_AllEco_Ortho <- YieldQuantile_bs(DrawSet, LabFEVD_Ortho[nn], ndraws, quants, Horiz, FacDim,
          YieDim, DimLabelsYields, ModelType,
          Ortho = TRUE
        )

        NumOutBounds_Yields[[LabFEVD[nn]]][[ModelType]]$Ortho <- NumOutBounds_AllEco_Ortho
      }
    }
  }
  return(NumOutBounds_Yields)
}
######################################################################################################
#' Extract graphs of interest (bootstrap version)
#'
#' @param InputsForOutputs list containing the desired inputs for the construction of IRFs, GIRFs, FEVDs, and GFEVDs
#' @param ModelType desired model type
#'
#' @keywords internal

WishGraphs_FEVDandGFEVD_Boot <- function(InputsForOutputs, ModelType) {
  # 1) Factors
  WishGraphFac <- c(
    InputsForOutputs[[ModelType]]$FEVD$WishGraphs$RiskFactorsBootstrap,
    InputsForOutputs[[ModelType]]$GFEVD$WishGraphs$RiskFactorsBootstrap
  )

  if (ModelType %in% c("JLL original", "JLL No DomUnit", "JLL joint Sigma")) {
    WishGraphFac_Ortho <- c(
      InputsForOutputs[[ModelType]]$FEVD$WishGraphsOrtho$RiskFactorsBootstrap,
      InputsForOutputs[[ModelType]]$GFEVD$WishGraphsOrtho$RiskFactorsBootstrap
    )
  } else {
    WishGraphFac_Ortho <- NULL
  }

  # 2) Yields
  WishGraphYields <- c(
    InputsForOutputs[[ModelType]]$FEVD$WishGraphs$YieldsBootstrap,
    InputsForOutputs[[ModelType]]$GFEVD$WishGraphs$YieldsBootstrap
  )

  if (any(ModelType %in% c("JLL original", "JLL No DomUnit", "JLL joint Sigma"))) {
    WishGraphiYields_Ortho <- c(
      InputsForOutputs[[ModelType]]$FEVD$WishGraphsOrtho$YieldsBootstrap,
      InputsForOutputs[[ModelType]]$GFEVD$WishGraphsOrtho$YieldsBootstrap
    )
  } else {
    WishGraphiYields_Ortho <- NULL
  }

  Out <- list(
    Fac = WishGraphFac, Fac_Ortho = WishGraphFac_Ortho, Yields = WishGraphYields,
    Yields_Ortho = WishGraphiYields_Ortho
  )

  return(Out)
}
