#' Generate the graphical outputs for the selected models (Point estimate)
#'
#' @param ModelType A character vector indicating the model type to be estimated.
#' @param ModelPara List of model parameter estimates (See the \code{\link{Optimization}} function)
#' @param NumOut list of computed outputs containing the model fit, IRFs, FEVDs, GIRFs, GFEVDs and Term Premia
#' @param InputsForOutputs list containing the desired inputs for the construction of the desired output
#' @param Economies A character vector containing the names of the economies included in the system.
#' @param FactorLabels A list of character vectors with labels for all variables in the model.
#' @param Folder2save Folder path where the outputs will be stored.
#' @param verbose Logical flag controlling function messaging.
#'
#' @importFrom ggplot2 ggplot theme_classic scale_x_date element_rect scale_fill_brewer theme_minimal scale_x_continuous scale_y_continuous ggplot_build ggplot_gtable expansion
#' @importFrom cowplot plot_grid ggdraw draw_label
#'
#' @keywords internal

GraphicalOutputs <- function(ModelType, ModelPara, NumOut, InputsForOutputs, Economies, FactorLabels,
                             Folder2save, verbose) {
  # Generate the graph paths and the graph folders
  dirs <- c("Outputs", paste0("Outputs/", ModelType), paste0("Outputs/", ModelType, "/Point Estimate"))
  lapply(file.path(Folder2save, dirs), dir.create, recursive = TRUE, showWarnings = FALSE)

  PathsGraphs <- FolderCreationPoint(ModelType, Economies, Folder2save)

  # 1) Plot the set of risk factors
  if (verbose) message("2.3) Generating the graphs of interest")
  RiskFactorsGraphs(
    ModelType, InputsForOutputs[[ModelType]]$RiskFactors$WishGraphs, ModelPara, Economies,
    FactorLabels, Folder2save, verbose
  )

  if (ModelType %in% c("JPS original", "JPS global", "GVAR single")) {
    for (i in 1:length(Economies)) {
      dir.create(paste(Folder2save, "/Outputs/", ModelType, "/Point Estimate/Model ", Economies[i], sep = ""))
    }
  }

  # 2) Model fit
  Fitgraphs(
    ModelType, InputsForOutputs[[ModelType]]$Fit$WishGraphs, ModelPara, NumOut, Economies,
    PathsGraphs, Folder2save, verbose
  )

  # 3) IRF and GIRF
  if (any(ModelType == c("JLL original", "JLL No DomUnit", "JLL joint Sigma"))) {
    OutType <- c("IRF", "GIRF", "IRF Ortho", "GIRF Ortho")
  } else {
    OutType <- c("IRF", "GIRF")
  }

  for (j in 1:length(OutType)) {
    WG <- Wished_Graphs_IRFandGIRF(InputsForOutputs, OutType[j], ModelType)
    IRFandGIRFgraphs(ModelType, NumOut, WG$RiskGraphs, WG$YieldGraphs, InputsForOutputs[[ModelType]]$IRF$horiz,
      PathsGraphs,
      OutputType = OutType[j], Economies, Folder2save, verbose
    )
  }

  # 4) FEVD and GFEVD
  if (any(ModelType == c("JLL original", "JLL No DomUnit", "JLL joint Sigma"))) {
    OutType <- c("FEVD", "GFEVD", "FEVD Ortho", "GFEVD Ortho")
  } else {
    OutType <- c("FEVD", "GFEVD")
  }

  for (j in 1:length(OutType)) {
    WG <- Wished_Graphs_FEVDandGFEVD(InputsForOutputs, OutType[j], ModelType)
    FEVDandGFEVDgraphs(
      ModelType, NumOut, WG$RiskGraphs, WG$YieldGraphs, InputsForOutputs[[ModelType]]$FEVD$horiz,
      PathsGraphs, OutType[j], Economies, Folder2save, verbose
    )
  }

  # 5) Term premia decomposition
  TPDecompGraph(
    ModelType, NumOut, ModelPara, InputsForOutputs[[ModelType]]$RiskPremia$WishGraphs,
    InputsForOutputs$UnitMatYields, Economies, PathsGraphs, Folder2save, verbose
  )

  if (verbose) {
    message(paste("Desired graphs are saved in your chosen directory. Please, check:", Folder2save, "\n"))
  }
}

#########################################################################################################
#################################### RISK FACTORS GRAPHS ################################################
#########################################################################################################
#' Spanned and unspanned factors plot
#'
#' @param ModelType character. Estimated model type. Permissible choices: "JPS original", "JPS global", "GVAR single", "JPS multi", "GVAR multi", "JLL original", "JLL No DomUnit", "JLL joint Sigma".
#' @param WishRFgraphs logical. Set TRUE to generate graphs, FALSE otherwise.
#' @param ModelOutputs list. Model parameter estimates (see \code{\link{Optimization}}).
#' @param Economies character vector. Names of the \code{C} economies included in the system.
#' @param FactorLabels list. Labels for all variables in the model.
#' @param Folder2save character. Folder path where the outputs will be stored.
#' @param verbose logical. Flag controlling function messaging.
#'
#' @examples
#' data("ParaSetEx")
#' # Adapt factor labels according to the example
#' ModelType <- "JPS original"
#' Economy <- "Brazil"
#' FacLab <- LabFac(N = 1, DomVar = "Eco_Act", GlobalVar = "Gl_Eco_Act", Economy, ModelType)
#'
#' RiskFactorsGraphs(ModelType,
#'   WishRFgraphs = FALSE, ParaSetEx, Economy, FacLab,
#'   Folder2save = NULL, verbose = FALSE
#' )
#'
#' @section Available Methods:
#' - `autoplot(object, type = "RiskFactors")`
#'
#' @export

RiskFactorsGraphs <- function(ModelType, WishRFgraphs, ModelOutputs, Economies, FactorLabels, Folder2save, verbose) {
  if (!WishRFgraphs) {
    if (verbose) message("No graphs for risk factor dynamics were generated")
    return(invisible(NULL))
  }

  if (verbose) message(" ** Risk Factor dynamics")
  # Custom ggplot theme
  theme_custom <- function() {
    theme_classic() +
      theme(
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_text(size = 10)
      )
  }

  G <- length(FactorLabels$Global)
  N <- length(FactorLabels$Spanned)
  M <- length(FactorLabels$Domestic) - N
  C <- length(Economies)

  SepQ_Lab <- c("JPS original", "JPS global", "GVAR single")

  # Extract Factors
  if (ModelType %in% SepQ_Lab) {
    FactorList <- lapply(Economies, function(e) {
      X <- ModelOutputs[[ModelType]][[e]]$Inputs$AllFactors
      if (e != Economies[1] && G > 0) X <- X[-seq_len(G), ] # Remove global factors for second economy onwards
      return(X)
    })
    Factors <- do.call(rbind, FactorList)
  } else {
    Factors <- ModelOutputs[[ModelType]]$Inputs$AllFactors
  }

  # Ensure column names match expected format
  FactorNames <- c(FactorLabels$Global, FactorLabels$Domestic)
  if (!is.null(Folder2save)) { # useful command for the implementation of the autoplot method
    Folder2save <- paste(Folder2save, "/Outputs", "/", ModelType, sep = "")
  }

  T_dim <- ncol(Factors) # Time periods

  # Convert matrix to a properly structured dataframe
  SS <- as.data.frame(t(Factors)) # Transpose correctly
  SS$TimeSpan <- seq_len(T_dim)
  SS <- SS[, c("TimeSpan", colnames(SS)[-ncol(SS)])] # Reorder so TimeSpan is first

  # Extract column names after transposition
  nmFactors <- colnames(SS)[-1]

  # Generate plots
  plot_list <- vector("list", G + N + M)
  plot_list_no_legend <- vector("list", G + N + M)

  # a) Global Factors
  TimeSpan <- Value <- Legend <- NULL # Define global variables for R CMD check

  for (i in seq_len(G)) {
    g <- ggplot(SS, aes(x = TimeSpan)) +
      geom_line(aes(y = get(nmFactors[i]))) +
      geom_hline(yintercept = 0) +
      ggtitle(FactorLabels$Global[i]) +
      theme_custom()

    plot_list[[i]] <- g
    plot_list_no_legend[[i]] <- g + theme(legend.position = "none")
  }

  # b) Domestic Factors
  IDX <- outer(1:(N + M), 1:C, function(j, i) (G + j) + (N + M) * (i - 1))

  for (j in seq_len(N + M)) {
    plot_data <- data.frame(TimeSpan = SS$TimeSpan)

    for (i in seq_len(C)) {
      plot_data[[Economies[i]]] <- SS[[nmFactors[IDX[j, i]]]]
    }

    plot_data_long <- stats::reshape(
      plot_data,
      direction = "long",
      varying = names(plot_data)[names(plot_data) != "TimeSpan"],
      v.names = "Value",
      timevar = "Legend",
      times = names(plot_data)[names(plot_data) != "TimeSpan"],
      idvar = "TimeSpan"
    )


    d <- ggplot(plot_data_long, aes(x = TimeSpan, y = Value, color = Legend)) +
      geom_line() +
      geom_hline(yintercept = 0) +
      ggtitle(FactorLabels$Domestic[j]) +
      theme_custom() +
      labs(colour = "Legend")

    plot_list[[G + j]] <- d
    plot_list_no_legend[[G + j]] <- d + theme(legend.position = "none")
  }

  # c) Get common legend
  # Get data from the first domestic factor plot
  legend_data <- plot_list[[G + 1]]$data
  # Format legend
  legend_plot <- ggplot(legend_data, aes(x = TimeSpan, y = Value, color = Legend)) +
    geom_line(linewidth = 1) + # Make lines slightly thicker for legend clarity
    labs(color = "Legend") +
    theme_minimal() +
    theme(
      legend.direction = "horizontal",
      legend.position = "bottom",
      legend.justification = "center",
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 8),
      legend.background = element_rect(color = "black", fill = "white"),
      axis.title = element_blank(),
      axis.text = element_blank(),
      panel.grid = element_blank()
    ) +
    scale_x_continuous(breaks = NULL) +
    scale_y_continuous(breaks = NULL)

  # Force build the plot to ensure legend is created
  built_plot <- ggplot_build(legend_plot)
  gtable_plot <- ggplot_gtable(built_plot)

  # Find the legend guide box
  guide_index <- which(sapply(gtable_plot$grobs, function(x) grepl("guide-box", x$name, fixed = TRUE)))
  LegendSubPlot <- gtable_plot$grobs[[guide_index[1]]]

  # Build the graph subplots:
  subplots <- plot_grid(plotlist = plot_list_no_legend, ncol = 3)

  # Create a new subplot with graphs + legend:
  FactorsGraph <- plot_grid(LegendSubPlot, subplots, ncol = 1, rel_heights = c(.1, 1))

  # Save final graph at the desried folder:
  if (!is.null(Folder2save)) {
    suppressMessages(ggplot2::ggsave(FactorsGraph, filename = paste0("RiskFactors", ".png"), path = Folder2save))
    print(FactorsGraph)
  }

  return(FactorsGraph)
}

######################################################################################################
##################################### 1) FIT ########################################################
#####################################################################################################
#' Model fit graphs for all models
#'
#' @param ModelType character. Estimated model type.
#' @param WishFitgraphs logical. Set TRUE to generate fit graphs, FALSE otherwise.
#' @param ModelPara list. Model parameter estimates (see \code{\link{Optimization}}).
#' @param NumOut list. Outputs containing model fit, IRFs, FEVDs, GIRFs, GFEVDs and Term premia.
#' @param Economies character vector. Names of the economies included in the system.
#' @param PathsGraphs character. Path of the folder in which the graphs will be saved.
#' @param Folder2save character. Desired folder path to save outputs.
#' @param verbose logical. Flag controlling function messaging.
#'
#' @importFrom ggplot2 ggplot theme_classic scale_x_date element_rect
#'
#' @examples
#' data("ParaSetEx")
#' data("NumOutEx")
#' ModelType <- "JPS original"
#' Economy <- "Brazil"
#' Fitgraphs(ModelType,
#'   WishFitgraphs = TRUE, ParaSetEx, NumOutEx, Economy, PathsGraphs = NULL,
#'   Folder2save = NULL, verbose = FALSE
#' )
#'
#' @section Available Methods:
#' - `autoplot(object, type = "Fit")`
#'
#' @export

Fitgraphs <- function(ModelType, WishFitgraphs, ModelPara, NumOut, Economies, PathsGraphs, Folder2save, verbose) {
  if (!WishFitgraphs) {
    if (verbose) message("No graphs for bond yields fit were generated")
    return(invisible(NULL))
  } else {
    if (verbose) message(" ** Fit of bond yields")

    C <- length(Economies)

    Autoplot_List <- list()

    # 1) Models estimated individually
    if (any(ModelType == c("JPS original", "JPS global", "GVAR single"))) {
      for (i in 1:C) {
        if (!is.null(Folder2save)) { # useful command for the implementation of the autoplot method
          dir.create(paste(Folder2save, "/Outputs/", ModelType, "/Point Estimate/Model ", Economies[i], "/Fit", sep = ""))
        }
        mat <- (ModelPara[[ModelType]][[Economies[i]]]$Inputs$mat) * 12
        J <- length(mat)

        YieldData <- ModelPara[[ModelType]][[Economies[i]]]$Inputs$Y
        ModelFit <- NumOut$Fit[[ModelType]][[Economies[i]]]$`Yield Fit`
        ModelImplied <- NumOut$Fit[[ModelType]][[Economies[i]]]$`Yield Model Implied`

        # Make graphs
        p <- Fit_Subplot(YieldData, ModelFit, ModelImplied, J, mat, Economies[i], ModelType, PathsGraphs)
        Autoplot_List[[Economies[i]]] <- p
      }
    } else {
      # 2) Models estimated jointly
      if (!is.null(Folder2save)) { # useful command for the implementation of the autoplot method
        dir.create(paste(Folder2save, "/Outputs/", ModelType, "/Point Estimate/Fit", sep = "")) # Folder creation
      }
      mat <- (ModelPara[[ModelType]]$Inputs$mat) * 12
      J <- length(mat)

      Idx0 <- 0
      for (i in 1:C) {
        Idx <- Idx0
        IdxGrpahs <- (Idx + 1):(Idx + J)

        YieldData <- ModelPara[[ModelType]]$Inputs$Y[IdxGrpahs, ]
        ModelFit <- NumOut$Fit$`Yield Fit`[IdxGrpahs, ]
        ModelImplied <- NumOut$Fit$`Yield Model Implied`[IdxGrpahs, ]

        # Make graphs
        p <- Fit_Subplot(YieldData, ModelFit, ModelImplied, J, mat, Economies[i], ModelType, PathsGraphs)
        Autoplot_List[[Economies[i]]] <- p

        Idx0 <- (Idx + J)
      }
    }
  }

  return(Autoplot_List)
}

######################################################################################################
##################################### 2) IRF and GIRF ################################################
#####################################################################################################
#' IRF and GIRF graphs for all models
#'
#' @param ModelType character. Estimated model type.Permissible choices: "JPS original", "JPS global", "GVAR single", "JPS multi", "GVAR multi", "JLL original", "JLL No DomUnit", "JLL joint Sigma".
#' @param NumOut list. Computed outputs containing model fit, IRFs, FEVDs, GIRFs, GFEVDs and term premia.
#' @param WishPdynamicsgraphs logical. Set TRUE to generate risk factor graphs, FALSE otherwise.
#' @param WishYieldsgraphs logical. Set TRUE to generate bond yield graphs, FALSE otherwise.
#' @param IRFhoriz integer. Desired horizon of analysis for the IRFs.
#' @param PathsGraphs character. Path of the folder in which the graphs will be saved.
#' @param OutputType character. Available options: "IRF", "GIRF", "IRF Ortho", "GIRF Ortho".
#' @param Economies character vector. Names of the \code{C} economies included in the system.
#' @param Folder2save character. Folder path where the outputs will be stored.
#' @param verbose logical. Flag controlling function messaging.
#'
#' @importFrom ggplot2 ggplot geom_line labs geom_hline theme ggtitle theme_classic
#'
#' @examples
#' data("NumOutEx")
#' ModelType <- "JPS original"
#' Economy <- "Brazil"
#' IRFhoriz <- 20
#' irf_Out <- IRFandGIRFgraphs(ModelType, NumOutEx,
#'   WishPdynamicsgraphs = FALSE, WishYieldsgraphs = TRUE, IRFhoriz,
#'   PathsGraphs = NULL, OutputType = "GIRF", Economy, Folder2save = NULL,
#'   verbose = FALSE
#' )
#'
#' @section Available Methods:
#' - `autoplot(object, type = "IRF_Factor")`, `autoplot(object, type = "IRF_Yields")`,
#'   `autoplot(object, type = "GIRF_Yields")`, `autoplot(object, type = "GIRF_Yields")`.
#'   For JLL-based models: `autoplot(object, type = "IRF_Factor-_Ortho")`,\cr
#'   `autoplot(object, type = "IRF_Yields_Ortho")`, `autoplot(object, type = "GIRF_Yields_Ortho")`,
#'   `autoplot(object, type = "GIRF_Yields_Ortho")`.
#'
#' @export

IRFandGIRFgraphs <- function(ModelType, NumOut, WishPdynamicsgraphs, WishYieldsgraphs, IRFhoriz,
                             PathsGraphs, OutputType, Economies, Folder2save, verbose) {
  if (!WishPdynamicsgraphs & !WishYieldsgraphs) {
    if (verbose) message(paste("No graphs for", OutputType, "were generated"))

    if (any(ModelType == c("JLL original", "JLL No DomUnit", "JLL joint Sigma")) &&
      any(ModelType == c("IRF Ortho", "GIRF Ortho"))) {
      if (verbose) message(paste("No graphs for", OutputType, "were generated (orthogonolized version) \n"))
    }
    return(invisible(NULL))
  } else {
    if (OutputType == "IRF") {
      if (verbose) message(" ** IRFs")
      Lab_Fac <- "IRFFactors_shock_to_"
      Lab_Yield <- "IRFYields_shock_to_"
    } else if (OutputType == "GIRF") {
      if (verbose) message(" ** GIRFs")
      Lab_Fac <- "GIRFFactors_shock_to_"
      Lab_Yield <- "GIRFYields_shock_to_"
    } else if (OutputType == "IRF Ortho") {
      if (verbose) message(" ** IRFs-Ortho")
      Lab_Fac <- "IRFFactors_shock_to_"
      Lab_Yield <- "IRFYields_shock_to_"
    } else {
      if (verbose) message(" ** GIRFs-Ortho")
      Lab_Fac <- "GIRFFactors_shock_to_"
      Lab_Yield <- "GIRFYields_shock_to_"
    }

    OutList <- list()

    ################ 1) Estimation done for countries individually ################
    if (any(ModelType == c("JPS original", "JPS global", "GVAR single"))) {
      K <- dim(NumOut$IRF[[ModelType]][[Economies[1]]]$Factors)[2] # Total number of risk factors
      C <- length(Economies)

      for (i in 1:C) {
        # a) Folder Creation
        if (!is.null(Folder2save)) {
          FolderPrep_IRFs(OutputType, WishPdynamicsgraphs, WishYieldsgraphs, Economies[i], ModelType, Folder2save)
        }
        # b) Recast IRFs or GIRFs
        IRFset <- BuildIRFlist(NumOut, Economies[i], ModelType, IRFhoriz, K, OutputType)
        nmFactors <- names(IRFset$IRFFactors[[1]])[-1] # Factor names
        nmYields <- names(IRFset$IRFYields[[1]])[-1] # Factor names

        # c) Graph of Factors
        if (WishPdynamicsgraphs == 1) {
          if (!is.null(Folder2save)) {
            PathAdj <- AdjustPathIRFs(OutputType, "Factors", PathsGraphs, Economies[i], ModelType)
          }

          for (h in 1:K) {
            p_list <- IRFandGIRFs_Format_Fac(IRFset$IRFFactors[[h]])
            subplots <- plot_grid(plotlist = p_list, ncol = 3) # Gather the graphs in sub-plots
            if (!is.null(Folder2save)) {
              suppressMessages(ggplot2::ggsave(subplots,
                file = paste0(Lab_Fac, nmFactors[[h]], "_Model_", Economies[i], ".png"),
                path = PathAdj
              )) # Save sub-plots
              print(subplots)
            }

            OutList[[Economies[i]]][[nmFactors[h]]] <- subplots
          }
        }

        # d) Graph of yields
        if (WishYieldsgraphs == 1) {
          if (!is.null(Folder2save)) {
            PathAdj <- AdjustPathIRFs(OutputType, "Yields", PathsGraphs, Economies[i], ModelType)
          }

          for (l in 1:K) {
            p_list <- IRFandGIRFs_Format_Yields(IRFset$IRFYields[[l]])
            subplots <- plot_grid(plotlist = p_list, ncol = 3) # Gather the graphs in sub-plots
            if (!is.null(Folder2save)) {
              suppressMessages(ggplot2::ggsave(subplots,
                file = paste0(Lab_Yield, nmFactors[[l]], "_Model_", Economies[i], ".png"),
                path = PathAdj
              )) # Save sub-plots
              print(subplots)
            }

            OutList[[Economies[i]]][[nmFactors[l]]] <- subplots
          }
        }
      }
    } else {
      ################ 2) Estimation done for countries jointly ###############################
      if (any(OutputType == c("IRF", "GIRF"))) {
        # Total number of risk factors
        if (any(ModelType == c("JLL original", "JLL No DomUnit", "JLL joint Sigma"))) {
          K <- dim(NumOut$IRF[[ModelType]]$Factors$Ortho)[2]
        } else {
          K <- dim(NumOut$IRF[[ModelType]]$Factors)[2]
        }

        # a) Folder Creation
        if (!is.null(Folder2save)) {
          FolderPrep_IRFs(OutputType, WishPdynamicsgraphs, WishYieldsgraphs, Economies, ModelType, Folder2save)
        }

        # b) Recast IRFs or GIRFs
        IRFset <- BuildIRFlist(NumOut, Economies, ModelType, IRFhoriz, K, OutputType)

        nmFactors <- names(IRFset$IRFFactors[[1]])[-1] # Factor names
        nmYields <- names(IRFset$IRFYields[[1]])[-1] # Factor names

        # c) Graph of Factors
        if (WishPdynamicsgraphs == 1) {
          if (!is.null(Folder2save)) {
            PathAdj <- AdjustPathIRFs(OutputType, "Factors", PathsGraphs, Economies, ModelType)
          }
          for (h in 1:K) {
            p_list <- IRFandGIRFs_Format_Fac(IRFset$IRFFactors[[h]])
            subplots <- plot_grid(plotlist = p_list, ncol = 3) # Gather the graphs in sub-plots
            if (!is.null(Folder2save)) {
              suppressMessages(ggplot2::ggsave(subplots,
                file = paste0(Lab_Fac, nmFactors[[h]], ".png"),
                path = PathAdj
              )) # Save sub-plots
              print(subplots)
            }

            OutList[[nmFactors[h]]] <- subplots
          }
        }

        # d) Graph of yields
        if (WishYieldsgraphs == 1) {
          if (!is.null(Folder2save)) {
            PathAdj <- AdjustPathIRFs(OutputType, "Yields", PathsGraphs, Economies, ModelType)
          }
          for (l in 1:K) {
            p_list <- IRFandGIRFs_Format_Yields(IRFset$IRFYields[[l]])
            subplots <- plot_grid(plotlist = p_list, ncol = 3) # Gather the graphs in sub-plots
            if (!is.null(Folder2save)) {
              suppressMessages(ggplot2::ggsave(subplots, file = paste0(Lab_Yield, nmFactors[[l]], ".png"), path = PathAdj)) # Save sub-plots
              print(subplots)
            }

            OutList[[nmFactors[l]]] <- subplots
          }
        }
      }

      ###################### 3) JLL orthoonalized ############################################################
      if (any(OutputType == c("IRF Ortho", "GIRF Ortho"))) {
        K <- dim(NumOut$IRF[[ModelType]]$Factors$Ortho)[2]

        # a) Folder Creation
        if (!is.null(Folder2save)) {
          FolderPrep_IRFs(OutputType, WishPdynamicsgraphs, WishYieldsgraphs, Economies, ModelType, Folder2save)
        }

        # b) Recast IRFs or GIRFs
        IRFset <- BuildIRFlist(NumOut, Economies, ModelType, IRFhoriz, K, OutputType)

        nmFactors <- names(IRFset$IRFFactors[[1]])[-1] # Factor names
        nmYields <- names(IRFset$IRFYields[[1]])[-1] # Factor names

        # a.2) Graph of Factors
        if (WishPdynamicsgraphs == 1) {
          if (!is.null(Folder2save)) {
            PathAdj <- AdjustPathIRFs(OutputType, "Factors", PathsGraphs, Economies, ModelType)
          }
          for (h in 1:K) {
            p_list <- IRFandGIRFs_Format_Fac(IRFset$IRFFactors[[h]])
            subplots <- plot_grid(plotlist = p_list, ncol = 3) # sub-plots
            if (!is.null(Folder2save)) {
              suppressMessages(ggplot2::ggsave(subplots,
                file = paste0(Lab_Fac, nmFactors[[h]], "ORTHO", ".png"),
                path = PathAdj
              ))
              print(subplots)
            }
            OutList[[nmFactors[h]]] <- subplots
          }
        }

        # d) Graph of yields
        if (WishYieldsgraphs == 1) {
          if (!is.null(Folder2save)) {
            PathAdj <- AdjustPathIRFs(OutputType, "Yields", PathsGraphs, Economies, ModelType)
          }
          for (l in 1:K) {
            p_list <- IRFandGIRFs_Format_Yields(IRFset$IRFYields[[l]])
            subplots <- plot_grid(plotlist = p_list, ncol = 3) # Gather the graphs in sub-plots
            if (!is.null(Folder2save)) {
              suppressMessages(ggplot2::ggsave(subplots,
                file = paste0(Lab_Yield, nmFactors[[l]], "ORTHO", ".png"),
                path = PathAdj
              )) # Save sub-plots
              print(subplots)
            }

            OutList[[nmFactors[l]]] <- subplots
          }
        }
      }
    }
  }

  return(OutList)
}

######################################################################################################
##################################### 3) FEVD and GFEVD ##############################################
#####################################################################################################
#' FEVD and GFEVD graphs for all models
#'
#' @param ModelType character. Estimated model type. Permissible choices: "JPS original", "JPS global", "GVAR single", "JPS multi", "GVAR multi", "JLL original", "JLL No DomUnit", "JLL joint Sigma".
#' @param NumOut list. Computed outputs containing model fit, IRFs, FEVDs, GIRFs, GFEVDs and Term premia.
#' @param WishPdynamicsgraphs logical. Set TRUE to generate risk factor graphs, FALSE otherwise.
#' @param WishYieldsgraphs logical. Set TRUE to generate bond yield graphs, FALSE otherwise.
#' @param FEVDhoriz integer. Desired horizon of analysis for the FEVDs.
#' @param PathsGraphs character. Path of the folder in which the graphs will be saved.
#' @param OutputType character. Available options: "FEVD", "GFEVD", "FEVD Ortho", "GFEVD Ortho".
#' @param Economies character vector. Names of the \code{C} economies included in the system.
#' @param Folder2save character. Folder path where the outputs will be stored.
#' @param verbose logical. Flag controlling function messaging.
#'
#' @importFrom ggplot2 ggplot theme_classic aes element_text theme labs ggtitle element_blank geom_col geom_hline
#'
#' @examples
#' data("NumOutEx")
#' ModelType <- "JPS original"
#' Economy <- "Brazil"
#' FEVDhoriz <- 20
#' FEVDandGFEVDgraphs(ModelType, NumOutEx,
#'   WishPdynamicsgraphs = FALSE, WishYieldsgraphs = TRUE, FEVDhoriz,
#'   PathsGraphs = NULL, OutputType = "FEVD", Economy,
#'   Folder2save = NULL, verbose = FALSE
#' )
#'
#' @section Available Methods:
#' - `autoplot(object, type = "FEVD_Factor")`, `autoplot(object, type = "FEVD_Yields")`,
#'   `autoplot(object, type = "GFEVD_Yields")`, `autoplot(object, type = "GFEVD_Yields")`.
#'   For JLL-based models: `autoplot(object, type = "FEVD_Factor-_Ortho")`, \cr
#'   `autoplot(object, type = "FEVD_Yields_Ortho")`, `autoplot(object, type = "GFEVD_Yields_Ortho")`,
#'   `autoplot(object, type = "GFEVD_Yields_Ortho")`.
#'
#' @export

FEVDandGFEVDgraphs <- function(ModelType, NumOut, WishPdynamicsgraphs, WishYieldsgraphs, FEVDhoriz,
                               PathsGraphs, OutputType, Economies, Folder2save, verbose) {
  if (!WishPdynamicsgraphs & !WishYieldsgraphs) {
    if (verbose) message(paste("No graphs for", OutputType, "were generated"))

    if (any(ModelType == c("JLL original", "JLL No DomUnit", "JLL joint Sigma")) &&
      any(ModelType == c("FEVD Ortho", "GFEVD Ortho"))) {
      if (verbose) message(paste("No graphs for", OutputType, "were generated (orthogonolized version)"))
    }
    return(invisible(NULL))
  } else {
    if (OutputType == "FEVD") {
      if (verbose) message(" ** FEVDs")
      Lab_Fac <- "FEVDFactors_"
      Lab_Yield <- "FEVDYields_"
    } else if (OutputType == "GFEVD") {
      if (verbose) message(" ** GFEVDs")
      Lab_Fac <- "GFEVDFactors_"
      Lab_Yield <- "GFEVDYields_"
    } else if (OutputType == "FEVD Ortho") {
      if (verbose) message(" ** FEVDs-Ortho")
      Lab_Fac <- "FEVDFactors_ORTHO_"
      Lab_Yield <- "FEVDYields_ORTHO_"
    } else {
      if (verbose) message(" ** GFEVDs-Ortho")
      Lab_Fac <- "GFEVDFactors_ORTHO_"
      Lab_Yield <- "GFEVDYields_ORTHO_"
    }

    Autoplot_List <- list()

    ################ 1) Estimation done for countries individually ################
    if (any(ModelType == c("JPS original", "JPS global", "GVAR single"))) {
      K <- dim(NumOut$FEVD[[ModelType]][[Economies[1]]]$Factors)[[2]] # Total number of risk factors
      J <- dim(NumOut$FEVD[[ModelType]][[Economies[1]]]$Yields)[3]

      C <- length(Economies)

      for (i in 1:C) {
        # a) Folder Creation
        if (!is.null(Folder2save)) {
          FolderPrep_FEVDs(OutputType, WishPdynamicsgraphs, WishYieldsgraphs, Economies[i], ModelType, Folder2save)
        }
        # b) Recast FEVDs or GFEVDs
        FEDset <- BuildFEVDlist(NumOut, Economies[i], ModelType, FEVDhoriz, K, J, OutputType)

        # c) Graph of Factors
        if (WishPdynamicsgraphs == 1) {
          PathAdj <- AdjustPathFEVDs(OutputType, "Factors", PathsGraphs, Economies[i], ModelType)
          nmFactors <- colnames(NumOut$FEVD[[ModelType]][[Economies[i]]]$Factors) # Factor names

          for (h in 1:K) {
            pp <- FEVDandGFEVDs_Graphs(
              OutputType, FEDset$FEVDFactors[[h]], nmFactors[h], Lab_Fac,
              PathAdj
            )

            Autoplot_List[[Economies[i]]][[nmFactors[h]]] <- pp
          }
        }

        # d) Graph of Yields
        if (WishYieldsgraphs == 1) {
          PathAdj <- AdjustPathFEVDs(OutputType, "Yields", PathsGraphs, Economies[i], ModelType)
          nmYields <- dimnames(NumOut$FEVD[[ModelType]][[Economies[i]]]$Yields)[[3]] # Yield labels

          for (h in 1:J) {
            pp <- FEVDandGFEVDs_Graphs(
              OutputType, FEDset$FEVDYields[[h]], nmYields[h], Lab_Yield,
              PathAdj
            )

            Autoplot_List[[Economies[i]]][[nmYields[h]]] <- pp
          }
        }
      }
    } else {
      ################ 2) Estimation done for countries jointly ###############################
      if (any(OutputType == c("FEVD", "GFEVD"))) {
        # Total number of risk factors
        if (any(ModelType == c("JLL original", "JLL No DomUnit", "JLL joint Sigma"))) {
          K <- dim(NumOut$FEVD[[ModelType]]$Factors$NonOrtho)[2]
          CJ <- dim(NumOut$GFEVD[[ModelType]]$Yields$NonOrtho)[3]
          nmFactors <- colnames(NumOut$FEVD[[ModelType]]$Factors$NonOrtho)
          nmYields <- dimnames(NumOut$FEVD[[ModelType]]$Yields$NonOrtho)[[3]]
        } else {
          K <- dim(NumOut$FEVD[[ModelType]]$Factors)[[2]]
          CJ <- dim(NumOut$FEVD[[ModelType]]$Yields)[[3]]
          nmFactors <- colnames(NumOut$FEVD[[ModelType]]$Factors) # Factor names
          nmYields <- dimnames(NumOut$FEVD[[ModelType]]$Yields)[[3]] # Yield labels
        }

        # a) Folder Creation
        if (!is.null(Folder2save)) {
          FolderPrep_FEVDs(OutputType, WishPdynamicsgraphs, WishYieldsgraphs, Economies, ModelType, Folder2save)
        }

        # b) Recast FEVDs or GFEVDs
        FEDset <- BuildFEVDlist(NumOut, Economies, ModelType, FEVDhoriz, K, CJ, OutputType)

        # c) Graph of Factors
        if (WishPdynamicsgraphs == 1) {
          PathAdj <- AdjustPathFEVDs(OutputType, "Factors", PathsGraphs, Economies, ModelType)

          for (h in 1:K) {
            pp <- FEVDandGFEVDs_Graphs(
              OutputType, FEDset$FEVDFactors[[h]], nmFactors[h], Lab_Fac,
              PathAdj
            )
            Autoplot_List[[nmFactors[h]]] <- pp
          }
        }

        # d) Graph of Yields
        if (WishYieldsgraphs == 1) {
          PathAdj <- AdjustPathFEVDs(OutputType, "Yields", PathsGraphs, Economies, ModelType)
          for (h in 1:CJ) {
            pp <- FEVDandGFEVDs_Graphs(
              OutputType, FEDset$FEVDYields[[h]], nmYields[h], Lab_Yield,
              PathAdj
            )
            Autoplot_List[[nmYields[h]]] <- pp
          }
        }
      }

      ###################### 3) JLL orthoonalized ############################################################
      if (any(OutputType == c("FEVD Ortho", "GFEVD Ortho"))) {
        # Total number of risk factors
        K <- dim(NumOut$FEVD[[ModelType]]$Factors$Ortho)[2]
        CJ <- dim(NumOut$GFEVD[[ModelType]]$Yields$Ortho)[3]
        nmFactors <- colnames(NumOut$FEVD[[ModelType]]$Factors$Ortho)
        nmYields <- dimnames(NumOut$FEVD[[ModelType]]$Yields$Ortho)[[3]]

        # a) Folder Creation
        if (!is.null(Folder2save)) {
          FolderPrep_FEVDs(OutputType, WishPdynamicsgraphs, WishYieldsgraphs, Economies, ModelType, Folder2save)
        }

        # b) Recast FEVDs or GFEVDs
        FEDset <- BuildFEVDlist(NumOut, Economies, ModelType, FEVDhoriz, K, CJ, OutputType)

        # c) Graph of Factors
        if (WishPdynamicsgraphs == 1) {
          PathAdj <- AdjustPathFEVDs(OutputType, "Factors", PathsGraphs, Economies, ModelType)
          for (h in 1:K) {
            pp <- FEVDandGFEVDs_Graphs(
              OutputType, FEDset$FEVDFactors[[h]], nmFactors[h], Lab_Fac,
              PathAdj
            )
            Autoplot_List[[nmFactors[h]]] <- pp
          }
        }

        # d) Graph of Yields
        if (WishYieldsgraphs == 1) {
          PathAdj <- AdjustPathFEVDs(OutputType, "Yields", PathsGraphs, Economies, ModelType)
          for (h in 1:CJ) {
            pp <- FEVDandGFEVDs_Graphs(
              OutputType, FEDset$FEVDYields[[h]], nmYields[h], Lab_Yield,
              PathAdj
            )
            Autoplot_List[[nmYields[h]]] <- pp
          }
        }
      }
    }
  }

  return(Autoplot_List)
}

#####################################################################################################################
######################################## 4) TERM PREMIA DECOMPOSITION ##################################################
#####################################################################################################################
#' Term Premia decomposition graphs for all models
#'
#' @param ModelType character. Estimated model type. Permissible choices: "JPS original", "JPS global", "GVAR single", "JPS multi", "GVAR multi", "JLL original", "JLL No DomUnit", "JLL joint Sigma".
#' @param NumOut list. Computed outputs containing model fit, IRFs, FEVDs, GIRFs, GFEVDs and risk premia.
#' @param ModelPara list. Model parameter estimates (see \code{\link{Optimization}}).
#' @param WishRPgraphs logical. Set TRUE to generate term premia graphs, FALSE otherwise.
#' @param UnitYields character. "Month" if yields are in months, "Year" if in years.
#' @param Economies character vector. Names of the \code{C} economies included in the system.
#' @param PathsGraphs character. Path of the folder in which the graphs will be saved.
#' @param Folder2Save character. Folder path where the outputs will be stored.
#' @param verbose logical. Flag controlling function messaging.
#'
#' @examples
#' data("ParaSetEx")
#' data("NumOutEx")
#' ModelType <- "JPS original"
#' Economy <- "Brazil"
#' UnitYields <- "Month"
#' TPDecompGraph(ModelType, NumOutEx, ParaSetEx,
#'   WishRPgraphs = FALSE, UnitYields, Economy,
#'   PathsGraphs = NULL, Folder2Save = NULL, verbose = FALSE
#' )
#'
#' @section Available Methods:
#' - `autoplot(object, type = "TermPremia")`
#'
#' @export

TPDecompGraph <- function(ModelType, NumOut, ModelPara, WishRPgraphs, UnitYields, Economies, PathsGraphs,
                          Folder2Save, verbose) {
  if (!WishRPgraphs) {
    if (verbose) message("No graphs for term-premia were generated")
    return(invisible(NULL))
  }

  if (verbose) message(" ** Term premia")

  # Function to create a single plot
  create_plot <- function(Data, ExpCom, TP, TimeSpan, title) {
    DataGraph <- data.frame(Data, ExpCom, TP, TimeSpan)
    ggplot(DataGraph) +
      geom_line(aes(x = TimeSpan, y = Data, color = GraphLegend[1])) +
      geom_line(aes(x = TimeSpan, y = ExpCom, color = GraphLegend[2])) +
      geom_line(aes(x = TimeSpan, y = TP, color = GraphLegend[3])) +
      labs(color = "Legend") +
      scale_color_manual(values = c("black", "#0072B2", "#D55E00")) +
      ggtitle(title) +
      theme_classic() +
      theme(
        plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 6)
      ) +
      geom_hline(yintercept = 0)
  }

  # Preliminary work
  SepQ_Lab <- c("JPS original", "JPS global", "GVAR single")
  isSepQ <- ModelType %in% SepQ_Lab

  # Folder creation
  if (!is.null(Folder2Save)) {
    if (!isSepQ) {
      dir.create(paste(Folder2Save, "/Outputs/", ModelType, "/Point Estimate/TermPremia", sep = ""), showWarnings = FALSE)
    }
  }

  # Extract common parameters
  Para_Set <- ModelPara[[ModelType]]

  if (isSepQ) {
    dt <- Para_Set[[Economies[1]]]$Inputs$dt
    mat <- Para_Set[[Economies[1]]]$Inputs$mat
    T_dim <- ncol(Para_Set[[Economies[1]]]$Inputs$Y)
  } else {
    dt <- Para_Set$Inputs$dt
    mat <- Para_Set$Inputs$mat
    Yields <- Para_Set$Inputs$Y
    T_dim <- ncol(Yields)
  }

  J <- length(mat)
  C <- length(Economies)
  TPdecomp <- NumOut$TermPremiaDecomp

  # Adjust maturity units
  k <- ifelse(UnitYields == "Month", 12, 1)
  YLab <- ifelse(UnitYields == "Month", "Months", "Years")
  matAdjUnit <- mat * k

  # Graph title, legends, and axes
  Dates <- 1:T_dim
  matRPLab <- paste(matAdjUnit, YLab)
  GraphLegend <- c("Data", "Expected Component", "Term Premium")

  # Main loop for generating plots
  AutoplotExport <- list()

  for (i in 1:C) {
    if (isSepQ) {
      if (!is.null(Folder2Save)) {
        dir.create(paste(Folder2Save, "/Outputs/", ModelType, "/Point Estimate/Model ", Economies[i], "/TermPremia", sep = ""), showWarnings = FALSE)
      }
      GraphPath <- PathsGraphs[[ModelType]]$TermPremia[[Economies[i]]]
      Yields <- Para_Set[[Economies[i]]]$Inputs$Y
      YieldData <- t(Yields) * 100
    } else {
      GraphPath <- PathsGraphs[[ModelType]]$TermPremia
      IdxRP <- 1:J + J * (i - 1)
      YieldData <- t(Yields[IdxRP, ] * 100)
    }

    plot_list_RP <- list()
    plot_list_no_legend_RP <- list()

    for (h in 1:length(matRPLab)) {
      Data <- YieldData[, h]
      ExpCom <- TPdecomp$RiskPremia[["Expected Component"]][[Economies[i]]][, h]
      TP <- TPdecomp$RiskPremia[["Term Premia"]][[Economies[i]]][, h]

      p <- create_plot(Data, ExpCom, TP, Dates, matRPLab[h])
      plot_list_RP[[h]] <- p
      plot_list_no_legend_RP[[h]] <- p + theme(legend.position = "none")
    }

    # Common title for each model
    title <- ggdraw() + draw_label(Economies[i], fontface = "bold", size = 14)

    # Generate legend
    built_plot <- ggplot_build(
      plot_list_RP[[1]] + # Start with your original plot
        theme(
          legend.direction = "horizontal", legend.position = "bottom", legend.justification = "center",
          legend.title = element_blank(),
          legend.text = element_text(size = 8)
        )
    )

    gtable_plot <- ggplot_gtable(built_plot)

    # Find the legend guide box
    guide_index <- which(sapply(gtable_plot$grobs, function(x) grepl("guide-box", x$name, fixed = TRUE)))
    LegendSubPlot_RP <- gtable_plot$grobs[[guide_index[1]]]

    # Build the graph subplots
    subplots_CS <- plot_grid(plotlist = plot_list_no_legend_RP, ncol = 3)
    subplots2save <- plot_grid(LegendSubPlot_RP, subplots_CS, ncol = 1, rel_heights = c(0.1, 1))

    if (!is.null(Folder2Save)) {
      suppressMessages(ggplot2::ggsave(subplots2save,
        file = paste0("TermPremia_", Economies[i], "_", ModelType, ".png"),
        path = GraphPath
      ))
      print(subplots2save)
    }
    AutoplotExport[[Economies[i]]] <- subplots2save
  }

  return(AutoplotExport)
}

######################################################################################################
##################################### AUXILIARY FUNCTIONS #############################################
######################################################################################################
#' Build subplot for fitted yields
#'
#' @param YieldData Time series of bond yields
#' @param ModelFit Time series of fitted bond yields
#' @param ModelImplied Time series of model-implied bond yields
#' @param MatLength number of country-specific maturities
#' @param mat vector of maturities
#' @param Economies Economies of the economic system
#' @param ModelType Desired estimated model
#' @param PathsGraphs Path to save the graphs
#'
#' @importFrom ggplot2 scale_color_manual
#'
#' @keywords internal

Fit_Subplot <- function(YieldData, ModelFit, ModelImplied, MatLength, mat, Economies, ModelType, PathsGraphs) {
  if (any(ModelType == c("JPS original", "JPS global", "GVAR single"))) {
    AdjPath <- PathsGraphs[[ModelType]]$Fit[[Economies]]
  } else {
    AdjPath <- PathsGraphs[[ModelType]]$Fit
  }

  # Graph titles
  matmonths <- paste(mat, "months - ")
  GraphTitles <- paste(matmonths, Economies)

  GraphLegend <- c("Data", "Model Fit", "Model-Implied")

  T_dim <- ncol(YieldData)
  # Prepare plots
  p <- c()
  plot_list_FIT <- list()
  plot_list_no_legend_FIT <- list()

  for (j in 1:MatLength) {
    Data <- YieldData[j, -1]
    Modfit <- ModelFit[j, -1]
    ModImp <- ModelImplied[j, -1]
    TimeSpan <- 1:(T_dim - 1)

    YieldCompar <- data.frame(Data, Modfit, ModImp, TimeSpan)

    # Graph: Fit x Data
    p <- ggplot(data = YieldCompar, aes(x = TimeSpan)) +
      geom_line(aes(y = Data, color = GraphLegend[1]), linewidth = 0.5) +
      geom_line(aes(y = Modfit, color = GraphLegend[2]), linewidth = 0.5) +
      geom_line(aes(y = ModImp, color = GraphLegend[3]), linetype = "dashed") +
      labs(color = "Legend") +
      scale_color_manual(values = c("black", "#0072B2", "#D55E00")) +
      ggtitle(GraphTitles[j]) +
      theme_classic() +
      theme(
        plot.title = element_text(size = 10, face = "bold", hjust = 0.5), axis.title.x = element_blank(),
        axis.title.y = element_blank(), axis.text.x = element_text(size = 8)
      ) +
      geom_hline(yintercept = 0)

    plot_list_FIT[[j]] <- p # All plots with legends

    plot_list_no_legend_FIT[[j]] <- p + theme(legend.position = "none") # Hide legends from all plots
  }

  # Generate legend
  # Legend settings:
  built_plot <- ggplot_build(
    plot_list_FIT[[1]] + # Start with the original plot
      theme(
        legend.direction = "horizontal",
        legend.position = "bottom",
        legend.justification = "center",
        legend.title = element_blank(),
        legend.text = element_text(size = 8)
      )
  )

  gtable_plot <- ggplot_gtable(built_plot)

  # Find the legend guide box
  guide_index <- which(sapply(gtable_plot$grobs, function(x) grepl("guide-box", x$name, fixed = TRUE)))
  LegendSubPlot <- gtable_plot$grobs[[guide_index[1]]]

  # Build the graph subplots:
  subplots_FIT <- plot_grid(plotlist = plot_list_no_legend_FIT, ncol = 3)
  subplots2save_FIT <- plot_grid(LegendSubPlot, subplots_FIT, ncol = 1, rel_heights = c(.1, 1))
  if (!is.null(AdjPath)) {
    suppressMessages(ggplot2::ggsave(subplots2save_FIT,
      file = paste0("Yields_Fit_", Economies, ".png"),
      path = AdjPath
    ))
    print(subplots2save_FIT)
  }

  return(subplots2save_FIT)
}
#########################################################################################################
#' Create folders for storing IRFs and GIRFs
#'
#' @param OutputType available options are "IRF", "GIRF", "IRF Ortho" and "GIRF Ortho"
#' @param WishPdynamicsgraphs binary variable specifing whether the user whishes IRFs and/or GIRFs of risk factors
#' @param WishYieldsgraphs binary variable specifing whether the user whishes IRFs and/or GIRFs of bond yields
#' @param Economies Set of economies that are part of the economic system
#' @param ModelType Desired modem type
#' @param Folder2Save Folder path where the outputs will be stored.
#'
#' @keywords internal

FolderPrep_IRFs <- function(OutputType, WishPdynamicsgraphs, WishYieldsgraphs, Economies, ModelType, Folder2Save) {
  ################# 1) SINGLE-COUNTRY CONTRY MODELS #################
  if (any(ModelType == c("JPS original", "JPS global", "GVAR single"))) {
    # A) General Folder Creation
    if (OutputType == "IRF") {
      dir.create(paste(Folder2Save, "/Outputs/", ModelType, "/Point Estimate/Model ", Economies, "/IRF", sep = ""))
    } else {
      dir.create(paste(Folder2Save, "/Outputs/", ModelType, "/Point Estimate/Model ", Economies, "/GIRF", sep = ""))
    }

    # a.1) Folders for graph of Factors
    if (WishPdynamicsgraphs == 1) {
      if (OutputType == "IRF") {
        dir.create(paste(Folder2Save, "/Outputs/", ModelType, "/Point Estimate/Model ", Economies,
          "/IRF/Factors",
          sep = ""
        ))
      } else {
        dir.create(paste(Folder2Save, "/Outputs/", ModelType, "/Point Estimate/Model ", Economies,
          "/GIRF/Factors",
          sep = ""
        ))
      }
    }

    # a.2) Folders for graph of yields
    if (WishYieldsgraphs == 1) {
      if (OutputType == "IRF") {
        dir.create(paste(Folder2Save, "/Outputs/", ModelType, "/Point Estimate/Model ", Economies,
          "/IRF/Yields",
          sep = ""
        ))
      } else {
        dir.create(paste(Folder2Save, "/Outputs/", ModelType, "/Point Estimate/Model ", Economies,
          "/GIRF/Yields",
          sep = ""
        ))
      }
    }
  } else if (any(ModelType == c("JPS multi", "GVAR multi"))) {
    ################# 2) JOINT COUNTRY CONTRY MODELS #################
    # A) General Folder Creation
    if (OutputType == "IRF") {
      dir.create(paste(Folder2Save, "/Outputs/", ModelType, "/Point Estimate/IRF", sep = ""))
    } else {
      dir.create(paste(Folder2Save, "/Outputs/", ModelType, "/Point Estimate/GIRF", sep = ""))
    }

    # a.1) Folders for graph of Factors
    if (WishPdynamicsgraphs == 1) {
      if (OutputType == "IRF") {
        dir.create(paste(Folder2Save, "/Outputs/", ModelType, "/Point Estimate/IRF/Factors", sep = ""))
      } else {
        dir.create(paste(Folder2Save, "/Outputs/", ModelType, "/Point Estimate/GIRF/Factors", sep = ""))
      }
    }

    # a.2) Folders for graph of yields
    if (WishYieldsgraphs == 1) {
      if (OutputType == "IRF") {
        dir.create(paste(Folder2Save, "/Outputs/", ModelType, "/Point Estimate/IRF/Yields", sep = ""))
      } else {
        dir.create(paste(Folder2Save, "/Outputs/", ModelType, "/Point Estimate/GIRF/Yields", sep = ""))
      }
    }

    ################# 3) JLL SPECIFIC MODELS #################
  } else {
    # A) General Folder Creation
    if (OutputType == "IRF") {
      dir.create(paste(Folder2Save, "/Outputs/", ModelType, "/Point Estimate/IRF", sep = ""))
    } else if (OutputType == "GIRF") {
      dir.create(paste(Folder2Save, "/Outputs/", ModelType, "/Point Estimate/GIRF", sep = ""))
    }

    ########### NON-ORTHO #####################################
    if (any(OutputType == c("IRF", "GIRF"))) {
      # a.1) Folders for graph of Factors
      if (WishPdynamicsgraphs == 1) {
        if (OutputType == "IRF") {
          dir.create(paste(Folder2Save, "/Outputs/", ModelType, "/Point Estimate/IRF/Factors", sep = ""))
        } else {
          dir.create(paste(Folder2Save, "/Outputs/", ModelType, "/Point Estimate/GIRF/Factors", sep = ""))
        }
      }

      # a.2) Folders for graph of yields
      if (WishYieldsgraphs == 1) {
        if (OutputType == "IRF") {
          dir.create(paste(Folder2Save, "/Outputs/", ModelType, "/Point Estimate/IRF/Yields", sep = ""))
        } else {
          dir.create(paste(Folder2Save, "/Outputs/", ModelType, "/Point Estimate/GIRF/Yields", sep = ""))
        }
      }
    } else {
      ########### ORTHO #####################################
      if (WishPdynamicsgraphs == 1) {
        if (OutputType == "IRF Ortho") {
          dir.create(paste(Folder2Save, "/Outputs/", ModelType, "/Point Estimate/IRF/Factors/Ortho", sep = ""))
        } else {
          dir.create(paste(Folder2Save, "/Outputs/", ModelType, "/Point Estimate/GIRF/Factors/Ortho", sep = ""))
        }
      }

      # a.2) Folders for graph of yields
      if (WishYieldsgraphs == 1) {
        if (OutputType == "IRF Ortho") {
          dir.create(paste(Folder2Save, "/Outputs/", ModelType, "/Point Estimate/IRF/Yields/Ortho", sep = ""))
        } else {
          dir.create(paste(Folder2Save, "/Outputs/", ModelType, "/Point Estimate/GIRF/Yields/Ortho", sep = ""))
        }
      }
    }
  }
}
#########################################################################################################
#' Generate paths to save IRFs/GIRFs graphs
#'
#' @param OutputType available options are "IRF" and "GIRF"
#' @param ElementType available options are "Factors" and "Yields"
#' @param PathsGraphs desired path to save the graphs
#' @param Economies Economies of the economic system
#' @param ModelType Desired estimated model type
#'
#' @keywords internal

AdjustPathIRFs <- function(OutputType, ElementType, PathsGraphs, Economies, ModelType) {
  # 1) Models estimated separately
  if (any(ModelType == c("JPS original", "JPS global", "GVAR single"))) {
    if (OutputType == "IRF") {
      if (ElementType == "Factors") {
        PathAdj <- PathsGraphs[[ModelType]]$IRF[[Economies]]$Factors
      } else {
        PathAdj <- PathsGraphs[[ModelType]]$IRF[[Economies]]$Yields
      } # Yields
    } else if (OutputType == "GIRF") {
      if (ElementType == "Factors") {
        PathAdj <- PathsGraphs[[ModelType]]$GIRF[[Economies]]$Factors
      } else {
        PathAdj <- PathsGraphs[[ModelType]]$GIRF[[Economies]]$Yields
      }
    }
  } else {
    # 2) Models estimated jointly
    if (OutputType == "IRF") { # IRF
      if (ElementType == "Factors") {
        PathAdj <- PathsGraphs[[ModelType]]$IRF$Factors
      } else {
        PathAdj <- PathsGraphs[[ModelType]]$IRF$Yields
      } # Yields
    } else if (OutputType == "GIRF") { # GIRF
      if (ElementType == "Factors") {
        PathAdj <- PathsGraphs[[ModelType]]$GIRF$Factors
      } else {
        PathAdj <- PathsGraphs[[ModelType]]$GIRF$Yields
      }

      # 3) Exclusively for JLL models
      # IRF ortho
    } else if (OutputType == "IRF Ortho") {
      if (ElementType == "Factors") {
        PathAdj <- PathsGraphs[[ModelType]]$IRF[["Factors Ortho"]]
      } else {
        PathAdj <- PathsGraphs[[ModelType]]$IRF[["Yields Ortho"]]
      } # Yields

      # GIRF ortho
    } else {
      if (ElementType == "Factors") {
        PathAdj <- PathsGraphs[[ModelType]]$GIRF[["Factors Ortho"]]
      } else {
        PathAdj <- PathsGraphs[[ModelType]]$GIRF[["Yields Ortho"]]
      } # Yields
    }
  }

  return(PathAdj)
}
##########################################################################################################
#' Build the list of IRF and GIRF for both factors and bond yields
#'
#' @param NumOut list of computed outputs containing the model fit, IRFs, FEVDs, GIRFs, GFEVDs and Term Premia
#' @param Economies Economies of the economic system
#' @param ModelType Desired model type
#' @param IRFhoriz time-horizon of the IRF and GIRF
#' @param FacDim dimension of the risk factor vector
#' @param OutputType available option are 'IRF' and 'GIRF'
#'
#' @keywords internal

BuildIRFlist <- function(NumOut, Economies, ModelType, IRFhoriz, FacDim, OutputType) {
  Horiz <- 0:(IRFhoriz - 1)
  IRFFactors <- list()
  IRFYields <- list()

  # 1) Extract IRFs
  if (OutputType == "IRF") {
    # a) Models estimated individually
    if (any(ModelType == c("JPS original", "JPS global", "GVAR single"))) {
      for (g in 1:FacDim) { # Recast the IRFs into a data-frame
        IRFFactors[[g]] <- data.frame(cbind(Horiz, NumOut$IRF[[ModelType]][[Economies]]$Factors[, , g]))
        IRFYields[[g]] <- data.frame(cbind(Horiz, NumOut$IRF[[ModelType]][[Economies]]$Yields[, , g]))
      }
      # b) Models estimated jointly
    } else if (any(ModelType == c("JPS multi", "GVAR multi"))) {
      for (g in 1:FacDim) { # Outputs in data-frame format
        IRFFactors[[g]] <- data.frame(cbind(Horiz, NumOut$IRF[[ModelType]]$Factors[, , g]))
        IRFYields[[g]] <- data.frame(cbind(Horiz, NumOut$IRF[[ModelType]]$Yields[, , g]))
      }
    } else {
      # c) JLL-based setups (Non-Orthogonalized version)
      for (g in 1:FacDim) {
        IRFFactors[[g]] <- data.frame(cbind(Horiz, NumOut$IRF[[ModelType]]$Factors$NonOrtho[, , g]))
        IRFYields[[g]] <- data.frame(cbind(Horiz, NumOut$IRF[[ModelType]]$Yields$NonOrtho[, , g]))
      }
    }

    # 2) Extract GIRFs
  } else if (OutputType == "GIRF") {
    # a) Models estimated individually
    if (any(ModelType == c("JPS original", "JPS global", "GVAR single"))) {
      for (g in 1:FacDim) { # Recast the IRFs into a data-frame
        IRFFactors[[g]] <- data.frame(cbind(Horiz, NumOut$GIRF[[ModelType]][[Economies]]$Factors[, , g]))
        IRFYields[[g]] <- data.frame(cbind(Horiz, NumOut$GIRF[[ModelType]][[Economies]]$Yields[, , g]))
      }
      # b) Models estimated jointly
    } else if (any(ModelType == c("JPS multi", "GVAR multi"))) {
      for (g in 1:FacDim) { # Outputs in data-frame format
        IRFFactors[[g]] <- data.frame(cbind(Horiz, NumOut$GIRF[[ModelType]]$Factors[, , g]))
        IRFYields[[g]] <- data.frame(cbind(Horiz, NumOut$GIRF[[ModelType]]$Yields[, , g]))
      }
      # c) JLL-based setups (non-Orthogonalized version)
    } else {
      for (g in 1:FacDim) {
        IRFFactors[[g]] <- data.frame(cbind(Horiz, NumOut$GIRF[[ModelType]]$Factors$NonOrtho[, , g]))
        IRFYields[[g]] <- data.frame(cbind(Horiz, NumOut$GIRF[[ModelType]]$Yields$NonOrtho[, , g]))
      }
    }

    # JLL-based setups (orthogonalized version)
  } else if (OutputType == "IRF Ortho") {
    for (g in 1:FacDim) { # Recast outputs in a data-frame
      IRFFactors[[g]] <- data.frame(cbind(Horiz, NumOut$IRF[[ModelType]]$Factors$Ortho[, , g]))
      IRFYields[[g]] <- data.frame(cbind(Horiz, NumOut$IRF[[ModelType]]$Yields$Ortho[, , g]))
    }
  } else { # GIRF Ortho
    for (g in 1:FacDim) { # Recast outputs in a data-frame
      IRFFactors[[g]] <- data.frame(cbind(Horiz, NumOut$GIRF[[ModelType]]$Factors$Ortho[, , g]))
      IRFYields[[g]] <- data.frame(cbind(Horiz, NumOut$GIRF[[ModelType]]$Yields$Ortho[, , g]))
    }
  }

  nmFactors <- names(IRFFactors[[1]]) # Factor names
  nmYields <- names(IRFYields[[1]]) # Yield names

  return(list(IRFFactors = IRFFactors, IRFYields = IRFYields))
}
##########################################################################################################
#' Gather data for IRFs and GIRFs grahs (version "Factors")
#'
#' @param IRFFac Data-frame with basic features of a single IRF for risk factors
#'
#' @keywords internal

IRFandGIRFs_Format_Fac <- function(IRFFac) {
  nmFactors <- names(IRFFac) # Factor names
  K <- length(nmFactors) - 1
  plot_list <- list()


  for (x in 1:K) {
    p <- c()
    mapping <- stats::setNames(list(IRFFac[[nmFactors[1]]], IRFFac[[nmFactors[x + 1]]]),
                              c("x", "y"))
    p <- ggplot(IRFFac, do.call(aes, mapping)) +
     geom_line() +
      ggtitle(nmFactors[x + 1])
    p <- p + theme_classic() + theme(
      plot.title = element_text(size = 6, face = "bold", hjust = 0.5),
      axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size = 8),
      axis.text.y = element_text(size = 4)
    )
    plot_list[[x]] <- p
  }

  return(plot_list)
}
#####################################################################################################
#' Gather data for IRFs and GIRFs grahs (version "Yields")
#'
#' @param IRFYields Data-frame with basic features of a single IRF for yields
#'
#' @keywords internal

IRFandGIRFs_Format_Yields <- function(IRFYields) {
  nmYields <- names(IRFYields) # Yields names
  J <- length(nmYields) - 1

  plot_list <- list()

  for (x in 1:J) { # Generate graph-by-graph
    p <- c()

    mapping <- stats::setNames(list(IRFYields[[nmYields[1]]], IRFYields[[nmYields[x + 1]]]),
                               c("x", "y"))

    p <- ggplot(IRFYields, do.call(aes, mapping)) +
      geom_line() +
      labs(title = nmYields[x + 1])
    p <- p + theme_classic() + theme(
      plot.title = element_text(size = 6, face = "bold", hjust = 0.5),
      axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size = 8),
      axis.text.y = element_text(size = 4)
    )
    plot_list[[x]] <- p
  }

  return(plot_list)
}
##########################################################################################################
#' Extract list of desired graph features (IRFs anc GIRFs)
#'
#' @param InputsForOutputs List of inputs for outputs
#' @param OutType Output types "IRF", "GIRF" and "IRF Ortho"
#' @param ModelType desired model type
#'
#' @keywords internal

Wished_Graphs_IRFandGIRF <- function(InputsForOutputs, OutType, ModelType) {
  # 1) JLL models
  if (any(ModelType == c("JLL original", "JLL No DomUnit", "JLL joint Sigma"))) {
    if (OutType == "IRF") {
      RiskGraphs <- InputsForOutputs[[ModelType]]$IRF$WishGraphs$RiskFactors
      YieldGraphs <- InputsForOutputs[[ModelType]]$IRF$WishGraphs$Yields
    } else if (OutType == "GIRF") {
      RiskGraphs <- InputsForOutputs[[ModelType]]$GIRF$WishGraphs$RiskFactors
      YieldGraphs <- InputsForOutputs[[ModelType]]$GIRF$WishGraphs$Yields
    } else if (OutType == "IRF Ortho") {
      RiskGraphs <- InputsForOutputs[[ModelType]]$IRF$WishGraphsOrtho$RiskFactors
      YieldGraphs <- InputsForOutputs[[ModelType]]$IRF$WishGraphsOrtho$Yields
    } else { # GIRF ortho
      RiskGraphs <- InputsForOutputs[[ModelType]]$GIRF$WishGraphsOrtho$RiskFactors
      YieldGraphs <- InputsForOutputs[[ModelType]]$GIRF$WishGraphsOrtho$Yields
    }
  } else {
    # 2) All other models
    if (OutType == "IRF") {
      RiskGraphs <- InputsForOutputs[[ModelType]]$IRF$WishGraphs$RiskFactors
      YieldGraphs <- InputsForOutputs[[ModelType]]$IRF$WishGraphs$Yields
    } else if (OutType == "GIRF") {
      RiskGraphs <- InputsForOutputs[[ModelType]]$GIRF$WishGraphs$RiskFactors
      YieldGraphs <- InputsForOutputs[[ModelType]]$GIRF$WishGraphs$Yields
    }
  }
  return(list(RiskGraphs = RiskGraphs, YieldGraphs = YieldGraphs))
}
#################################################################################################
#' Create folders for storing IRFs and GIRFs
#'
#' @param OutputType available options are "IRF", "GIRF", "IRF Ortho" and "GIRF Ortho"
#' @param WishPdynamicsgraphs binary variable specifing whether the user whishes IRFs and/or GIRFs of risk factors
#' @param WishYieldsgraphs binary variable specifing whether the user whishes IRFs and/or GIRFs of bond yields
#' @param Economies Set of economies that are part of the economic system
#' @param ModelType Desired modem type
#' @param Folder2Save Folder path where the outputs will be stored.
#'
#' @keywords internal

FolderPrep_FEVDs <- function(OutputType, WishPdynamicsgraphs, WishYieldsgraphs, Economies, ModelType, Folder2Save) {
  ################# 1) SINGLE-COUNTRY CONTRY MODELS #################
  if (any(ModelType == c("JPS original", "JPS global", "GVAR single"))) {
    # A) General Folder Creation
    if (OutputType == "FEVD") {
      dir.create(paste(Folder2Save, "/Outputs/", ModelType, "/Point Estimate/Model ", Economies, "/FEVD", sep = ""))
    } else {
      dir.create(paste(Folder2Save, "/Outputs/", ModelType, "/Point Estimate/Model ", Economies, "/GFEVD", sep = ""))
    }

    # a.1) Folders for graph of Factors
    if (WishPdynamicsgraphs == 1) {
      if (OutputType == "FEVD") {
        dir.create(paste(Folder2Save, "/Outputs/", ModelType, "/Point Estimate/Model ", Economies,
          "/FEVD/Factors",
          sep = ""
        ))
      } else {
        dir.create(paste(Folder2Save, "/Outputs/", ModelType, "/Point Estimate/Model ", Economies,
          "/GFEVD/Factors",
          sep = ""
        ))
      }
    }

    # a.2) Folders for graph of yields
    if (WishYieldsgraphs == 1) {
      if (OutputType == "FEVD") {
        dir.create(paste(Folder2Save, "/Outputs/", ModelType, "/Point Estimate/Model ", Economies,
          "/FEVD/Yields",
          sep = ""
        ))
      } else {
        dir.create(paste(Folder2Save, "/Outputs/", ModelType, "/Point Estimate/Model ", Economies,
          "/GFEVD/Yields",
          sep = ""
        ))
      }
    }
  } else if (any(ModelType == c("JPS multi", "GVAR multi"))) {
    ################# 2) JOINT COUNTRY CONTRY MODELS #################
    # A) General Folder Creation
    if (OutputType == "FEVD") {
      dir.create(paste(Folder2Save, "/Outputs/", ModelType, "/Point Estimate/FEVD", sep = ""))
    } else {
      dir.create(paste(Folder2Save, "/Outputs/", ModelType, "/Point Estimate/GFEVD", sep = ""))
    }

    # a.1) Folders for graph of Factors
    if (WishPdynamicsgraphs == 1) {
      if (OutputType == "FEVD") {
        dir.create(paste(Folder2Save, "/Outputs/", ModelType, "/Point Estimate/FEVD/Factors", sep = ""))
      } else {
        dir.create(paste(Folder2Save, "/Outputs/", ModelType, "/Point Estimate/GFEVD/Factors", sep = ""))
      }
    }

    # a.2) Folders for graph of yields
    if (WishYieldsgraphs == 1) {
      if (OutputType == "FEVD") {
        dir.create(paste(Folder2Save, "/Outputs/", ModelType, "/Point Estimate/FEVD/Yields", sep = ""))
      } else {
        dir.create(paste(Folder2Save, "/Outputs/", ModelType, "/Point Estimate/GFEVD/Yields", sep = ""))
      }
    }

    ################# 3) JLL SPECIFIC MODELS #################
  } else {
    # A) General Folder Creation
    if (OutputType == "FEVD") {
      dir.create(paste(Folder2Save, "/Outputs/", ModelType, "/Point Estimate/FEVD", sep = ""))
    } else if (OutputType == "GFEVD") {
      dir.create(paste(Folder2Save, "/Outputs/", ModelType, "/Point Estimate/GFEVD", sep = ""))
    }

    ########### NON-ORTHO #####################################
    if (any(OutputType == c("FEVD", "GFEVD"))) {
      # a.1) Folders for graph of Factors
      if (WishPdynamicsgraphs == 1) {
        if (OutputType == "FEVD") {
          dir.create(paste(Folder2Save, "/Outputs/", ModelType, "/Point Estimate/FEVD/Factors", sep = ""))
        } else {
          dir.create(paste(Folder2Save, "/Outputs/", ModelType, "/Point Estimate/GFEVD/Factors", sep = ""))
        }
      }

      # a.2) Folders for graph of yields
      if (WishYieldsgraphs == 1) {
        if (OutputType == "FEVD") {
          dir.create(paste(Folder2Save, "/Outputs/", ModelType, "/Point Estimate/FEVD/Yields", sep = ""))
        } else {
          dir.create(paste(Folder2Save, "/Outputs/", ModelType, "/Point Estimate/GFEVD/Yields", sep = ""))
        }
      }
    } else {
      ########### ORTHO #####################################
      if (WishPdynamicsgraphs == 1) {
        if (OutputType == "FEVD Ortho") {
          dir.create(paste(Folder2Save, "/Outputs/", ModelType, "/Point Estimate/FEVD/Factors/Ortho", sep = ""))
        } else {
          dir.create(paste(Folder2Save, "/Outputs/", ModelType, "/Point Estimate/GFEVD/Factors/Ortho", sep = ""))
        }
      }

      # a.2) Folders for graph of yields
      if (WishYieldsgraphs == 1) {
        if (OutputType == "FEVD Ortho") {
          dir.create(paste(Folder2Save, "/Outputs/", ModelType, "/Point Estimate/FEVD/Yields/Ortho", sep = ""))
        } else {
          dir.create(paste(Folder2Save, "/Outputs/", ModelType, "/Point Estimate/GFEVD/Yields/Ortho", sep = ""))
        }
      }
    }
  }
}
##########################################################################################################
#' Build the list of IRF and GIRF for both factors and bond yields
#'
#' @param NumOut list of computed outputs containing the model fit, IRFs, FEVDs, GIRFs, GFEVDs and Term Premia
#' @param Economies Economies of the economic system
#' @param ModelType Desired model type
#' @param FEVDhoriz time-horizon of the FEVD and GFEVD
#' @param FacDim dimension of the risk factor vector
#' @param YieldsDim dimension of the model set of yields
#' @param OutputType available option are 'FEVD' and 'GFEVD'
#'
#' @keywords internal

BuildFEVDlist <- function(NumOut, Economies, ModelType, FEVDhoriz, FacDim, YieldsDim, OutputType) {
  Horiz <- 0:(FEVDhoriz - 1)
  FEVDFactors <- list()
  FEVDYields <- list()
  index <- NULL

  # 1) Extract IRFs
  if (OutputType == "FEVD") {
    # a) Models estimated individually
    if (any(ModelType == c("JPS original", "JPS global", "GVAR single"))) {
      # Factors
      for (g in 1:FacDim) {
        FEVDFactors[[g]] <- data.frame(NumOut$FEVD[[ModelType]][[Economies]]$Factors[, , g])
        FEVDFactors[[g]]$index <- rownames(FEVDFactors[[g]])
        stacked_data <- utils::stack(FEVDFactors[[g]], select = -index)
        FEVDFactors[[g]] <- data.frame(
          index = FEVDFactors[[g]]$index,
          variable = stacked_data$ind,
          value = stacked_data$values
        )
        FEVDFactors[[g]]$index <- factor(FEVDFactors[[g]]$index, levels = 1:(FEVDhoriz))
      }

      # Yields
      for (g in 1:YieldsDim) {
        FEVDYields[[g]] <- data.frame(NumOut$FEVD[[ModelType]][[Economies]]$Yields[, , g])
        FEVDYields[[g]]$index <- rownames(FEVDYields[[g]])
        stacked_data <- utils::stack(FEVDYields[[g]], select = -index)
        FEVDYields[[g]] <- data.frame(
          index = FEVDYields[[g]]$index,
          variable = stacked_data$ind,
          value = stacked_data$values
        )
        FEVDYields[[g]]$index <- factor(FEVDYields[[g]]$index, levels = 1:(FEVDhoriz))
      }

      # b) Models estimated jointly
    } else if (any(ModelType == c("JPS multi", "GVAR multi"))) {
      # Factors
      for (g in 1:FacDim) {
        FEVDFactors[[g]] <- data.frame(NumOut$FEVD[[ModelType]]$Factors[, , g])
        FEVDFactors[[g]]$index <- rownames(FEVDFactors[[g]])
        stacked_data <- utils::stack(FEVDFactors[[g]], select = -index)
        FEVDFactors[[g]] <- data.frame(
          index = FEVDFactors[[g]]$index,
          variable = stacked_data$ind,
          value = stacked_data$values
        )
        FEVDFactors[[g]]$index <- factor(FEVDFactors[[g]]$index, levels = 1:(FEVDhoriz))
      }

      # Yields
      for (g in 1:YieldsDim) {
        FEVDYields[[g]] <- data.frame(NumOut$FEVD[[ModelType]]$Yields[, , g])
        FEVDYields[[g]]$index <- rownames(FEVDYields[[g]])
        stacked_data <- utils::stack(FEVDYields[[g]], select = -index)
        FEVDYields[[g]] <- data.frame(
          index = FEVDYields[[g]]$index,
          variable = stacked_data$ind,
          value = stacked_data$values
        )
        FEVDYields[[g]]$index <- factor(FEVDYields[[g]]$index, levels = 1:(FEVDhoriz))
      }
    } else {
      # c) JLL-based setups (Non-Orthogonalized version)
      # Factors
      for (g in 1:FacDim) {
        FEVDFactors[[g]] <- data.frame(NumOut$FEVD[[ModelType]]$Factors$NonOrtho[, , g])
        FEVDFactors[[g]]$index <- rownames(FEVDFactors[[g]])
        stacked_data <- utils::stack(FEVDFactors[[g]], select = -index)
        FEVDFactors[[g]] <- data.frame(
          index = FEVDFactors[[g]]$index,
          variable = stacked_data$ind,
          value = stacked_data$values
        )
        FEVDFactors[[g]]$index <- factor(FEVDFactors[[g]]$index, levels = 1:(FEVDhoriz))
      }

      # Yields
      for (g in 1:YieldsDim) {
        FEVDYields[[g]] <- data.frame(NumOut$FEVD[[ModelType]]$Yields$NonOrtho[, , g])
        FEVDYields[[g]]$index <- rownames(FEVDYields[[g]])
        stacked_data <- utils::stack(FEVDYields[[g]], select = -index)
        FEVDYields[[g]] <- data.frame(
          index = FEVDYields[[g]]$index,
          variable = stacked_data$ind,
          value = stacked_data$values
        )
        FEVDYields[[g]]$index <- factor(FEVDYields[[g]]$index, levels = 1:(FEVDhoriz))
      }
    }

    # 2) Extract GFEVD
  } else if (OutputType == "GFEVD") {
    # a) Models estimated individually
    if (any(ModelType == c("JPS original", "JPS global", "GVAR single"))) {
      for (g in 1:FacDim) {
        FEVDFactors[[g]] <- data.frame(NumOut$GFEVD[[ModelType]][[Economies]]$Factors[, , g])
        FEVDFactors[[g]]$index <- rownames(FEVDFactors[[g]])
        stacked_data <- utils::stack(FEVDFactors[[g]], select = -index)
        FEVDFactors[[g]] <- data.frame(
          index = FEVDFactors[[g]]$index,
          variable = stacked_data$ind,
          value = stacked_data$values
        )
        FEVDFactors[[g]]$index <- factor(FEVDFactors[[g]]$index, levels = 1:(FEVDhoriz))
      }

      # Yields
      for (g in 1:YieldsDim) {
        FEVDYields[[g]] <- data.frame(NumOut$GFEVD[[ModelType]][[Economies]]$Yields[, , g])
        FEVDYields[[g]]$index <- rownames(FEVDYields[[g]])
        stacked_data <- utils::stack(FEVDYields[[g]], select = -index)
        FEVDYields[[g]] <- data.frame(
          index = FEVDYields[[g]]$index,
          variable = stacked_data$ind,
          value = stacked_data$values
        )
        FEVDYields[[g]]$index <- factor(FEVDYields[[g]]$index, levels = 1:(FEVDhoriz))
      }

      # b) Models estimated jointly
    } else if (any(ModelType == c("JPS multi", "GVAR multi"))) {
      for (g in 1:FacDim) { # Outputs in data-frame format
        FEVDFactors[[g]] <- data.frame(NumOut$GFEVD[[ModelType]]$Factors[, , g])
        FEVDFactors[[g]]$index <- rownames(FEVDFactors[[g]])
        stacked_data <- utils::stack(FEVDFactors[[g]], select = -index)
        FEVDFactors[[g]] <- data.frame(
          index = FEVDFactors[[g]]$index,
          variable = stacked_data$ind,
          value = stacked_data$values
        )
        FEVDFactors[[g]]$index <- factor(FEVDFactors[[g]]$index, levels = 1:(FEVDhoriz))
      }

      # Yields
      for (g in 1:YieldsDim) {
        FEVDYields[[g]] <- data.frame(NumOut$GFEVD[[ModelType]]$Yields[, , g])
        FEVDYields[[g]]$index <- rownames(FEVDYields[[g]])
        stacked_data <- utils::stack(FEVDYields[[g]], select = -index)
        FEVDYields[[g]] <- data.frame(
          index = FEVDYields[[g]]$index,
          variable = stacked_data$ind,
          value = stacked_data$values
        )
        FEVDYields[[g]]$index <- factor(FEVDYields[[g]]$index, levels = 1:(FEVDhoriz))
      }

      # c) JLL-based setups (non-Orthogonalized version)
    } else {
      # Factors
      for (g in 1:FacDim) {
        FEVDFactors[[g]] <- data.frame(NumOut$GFEVD[[ModelType]]$Factors$NonOrtho[, , g])
        FEVDFactors[[g]]$index <- rownames(FEVDFactors[[g]])
        stacked_data <- utils::stack(FEVDFactors[[g]], select = -index)
        FEVDFactors[[g]] <- data.frame(
          index = FEVDFactors[[g]]$index,
          variable = stacked_data$ind,
          value = stacked_data$values
        )
        FEVDFactors[[g]]$index <- factor(FEVDFactors[[g]]$index, levels = 1:(FEVDhoriz))
      }

      # Yields
      for (g in 1:YieldsDim) {
        FEVDYields[[g]] <- data.frame(NumOut$GFEVD[[ModelType]]$Yields$NonOrtho[, , g])
        FEVDYields[[g]]$index <- rownames(FEVDYields[[g]])
        stacked_data <- utils::stack(FEVDYields[[g]], select = -index)
        FEVDYields[[g]] <- data.frame(
          index = FEVDYields[[g]]$index,
          variable = stacked_data$ind,
          value = stacked_data$values
        )
        FEVDYields[[g]]$index <- factor(FEVDYields[[g]]$index, levels = 1:(FEVDhoriz))
      }
    }

    # JLL-based setups (orthogonalized version)
  } else if (OutputType == "FEVD Ortho") {
    # Factors
    for (g in 1:FacDim) { # Recast outputs in a data-frame
      FEVDFactors[[g]] <- data.frame(NumOut$FEVD[[ModelType]]$Factors$Ortho[, , g])
      FEVDFactors[[g]]$index <- rownames(FEVDFactors[[g]])
      stacked_data <- utils::stack(FEVDFactors[[g]], select = -index)
      FEVDFactors[[g]] <- data.frame(
        index = FEVDFactors[[g]]$index,
        variable = stacked_data$ind,
        value = stacked_data$values
      )
      FEVDFactors[[g]]$index <- factor(FEVDFactors[[g]]$index, levels = 1:(FEVDhoriz))
    }

    # Yields
    for (g in 1:YieldsDim) {
      FEVDYields[[g]] <- data.frame(NumOut$FEVD[[ModelType]]$Yields$Ortho[, , g])
      FEVDYields[[g]]$index <- rownames(FEVDYields[[g]])
      stacked_data <- utils::stack(FEVDYields[[g]], select = -index)
      FEVDYields[[g]] <- data.frame(
        index = FEVDYields[[g]]$index,
        variable = stacked_data$ind,
        value = stacked_data$values
      )
      FEVDYields[[g]]$index <- factor(FEVDYields[[g]]$index, levels = 1:(FEVDhoriz))
    }
  } else { # GFEVD Ortho
    # Factors
    for (g in 1:FacDim) { # Recast outputs in a data-frame
      FEVDFactors[[g]] <- data.frame(NumOut$GFEVD[[ModelType]]$Factors$Ortho[, , g])
      FEVDFactors[[g]]$index <- rownames(FEVDFactors[[g]])
      stacked_data <- utils::stack(FEVDFactors[[g]], select = -index)
      FEVDFactors[[g]] <- data.frame(
        index = FEVDFactors[[g]]$index,
        variable = stacked_data$ind,
        value = stacked_data$values
      )
      FEVDFactors[[g]]$index <- factor(FEVDFactors[[g]]$index, levels = 1:(FEVDhoriz))
    }

    # Yields
    for (g in 1:YieldsDim) {
      FEVDYields[[g]] <- data.frame(NumOut$GFEVD[[ModelType]]$Yields$Ortho[, , g])
      FEVDYields[[g]]$index <- rownames(FEVDYields[[g]])
      stacked_data <- utils::stack(FEVDYields[[g]], select = -index)
      FEVDYields[[g]] <- data.frame(
        index = FEVDYields[[g]]$index,
        variable = stacked_data$ind,
        value = stacked_data$values
      )
      FEVDYields[[g]]$index <- factor(FEVDYields[[g]]$index, levels = 1:(FEVDhoriz))
    }
  }

  return(list(FEVDFactors = FEVDFactors, FEVDYields = FEVDYields))
}

##########################################################################################################
#' Generate paths to save IRFs/GIRFs graphs
#'
#' @param OutputType available options are "IRF" and "GIRF"
#' @param ElementType available options are "Factors" and "Yields"
#' @param PathsGraphs desired path to save the graphs
#' @param Economies Economies of the economic system
#' @param ModelType Desired estimated model type
#'
#' @keywords internal

AdjustPathFEVDs <- function(OutputType, ElementType, PathsGraphs, Economies, ModelType) {
  # 1) Models estimated seperatly
  if (any(ModelType == c("JPS original", "JPS global", "GVAR single"))) {
    if (OutputType == "FEVD") {
      if (ElementType == "Factors") {
        PathAdj <- PathsGraphs[[ModelType]]$FEVD[[Economies]]$Factors
      } else {
        PathAdj <- PathsGraphs[[ModelType]]$FEVD[[Economies]]$Yields
      } # Yields
    } else {
      if (ElementType == "Factors") {
        PathAdj <- PathsGraphs[[ModelType]]$GFEVD[[Economies]]$Factors
      } else {
        PathAdj <- PathsGraphs[[ModelType]]$GFEVD[[Economies]]$Yields
      }
    }

    # 2) Models estimated jointly
  } else {
    if (OutputType == "FEVD") { # FEVD
      if (ElementType == "Factors") {
        PathAdj <- PathsGraphs[[ModelType]]$FEVD$Factors
      } else {
        PathAdj <- PathsGraphs[[ModelType]]$FEVD$Yields
      } # Yields
    } else if (OutputType == "GFEVD") { # GFEVD
      if (ElementType == "Factors") {
        PathAdj <- PathsGraphs[[ModelType]]$GFEVD$Factors
      } else {
        PathAdj <- PathsGraphs[[ModelType]]$GFEVD$Yields
      }

      # 3) Exclusively for JLL models
      # FEVD ortho
    } else if (OutputType == "FEVD Ortho") {
      if (ElementType == "Factors") {
        PathAdj <- PathsGraphs[[ModelType]]$FEVD[["Factors Ortho"]]
      } else {
        PathAdj <- PathsGraphs[[ModelType]]$FEVD[["Yields Ortho"]]
      } # Yields

      # GFEVD ortho
    } else {
      if (ElementType == "Factors") {
        PathAdj <- PathsGraphs[[ModelType]]$GFEVD[["Factors Ortho"]]
      } else {
        PathAdj <- PathsGraphs[[ModelType]]$GFEVD[["Yields Ortho"]]
      } # Yields
    }
  }

  return(PathAdj)
}
##################################################################################################
#' Generates graphs for FEVDs and GFEVDs
#'
#' @param OutputType available options are "FEVD", "GFEVD", "FEVD Ortho" and "GFEVD Ortho"
#' @param FEVDlist list of FEVD and GFEVD outputs
#' @param nmVarInt name of variable of interest. Options: "Factors" and "Yields"
#' @param Lab_Fac label of the model factors
#' @param PathsGraphs Path to save graphs
#'
#' @keywords internal

FEVDandGFEVDs_Graphs <- function(OutputType, FEVDlist, nmVarInt, Lab_Fac, PathsGraphs) {
  # Merge factors with a small contributions
  FEVDlist <- MergeFEVD_graphs(FEVDlist, Threshold = 0.05)
  # Make variable labels more readable
  FEVDlist$variable <- factor(FEVDlist$variable,
                              labels = clean_labels(levels(factor(FEVDlist$variable))))

  index <- FEVDlist$index
  value <- FEVDlist$value
  variable <- FEVDlist$variable

  p <- ggplot(FEVDlist, aes(x = index, y = value, fill = variable)) +
    geom_col(position = "dodge") +
    scale_fill_brewer(palette = "Set2") + # color-friendly pallet
    scale_y_continuous(limits = c(0, NA), expand = expansion()) +
    labs(title = paste0(OutputType, " - ", nmVarInt)) +
    theme_classic() +
    theme(
      axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size = 8),
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5)
    )
  if (!is.null(PathsGraphs)) {
    suppressMessages(ggplot2::ggsave(p,
      file = paste0(Lab_Fac, nmVarInt, ".png"),
      path = PathsGraphs, width = 11, height = 7.98
    ))
    print(p)
  }
  return(p)
}
#############################################################################################
#' Extract list of desired graph features (IRFs anc GIRFs)
#'
#' @param InputsForOutputs List of inputs for outputs
#' @param OutType Output types "FEVD", "GFEVD" and "FEVD Ortho"
#' @param ModelType desired model type
#'
#' @keywords internal

Wished_Graphs_FEVDandGFEVD <- function(InputsForOutputs, OutType, ModelType) {
  # 1) JLL models
  if (ModelType %in% c("JLL original", "JLL No DomUnit", "JLL joint Sigma")) {
    if (OutType == "FEVD") {
      RiskGraphs <- InputsForOutputs[[ModelType]]$FEVD$WishGraphs$RiskFactors
      YieldGraphs <- InputsForOutputs[[ModelType]]$FEVD$WishGraphs$Yields
    } else if (OutType == "GFEVD") {
      RiskGraphs <- InputsForOutputs[[ModelType]]$GFEVD$WishGraphs$RiskFactors
      YieldGraphs <- InputsForOutputs[[ModelType]]$GFEVD$WishGraphs$Yields
    } else if (OutType == "FEVD Ortho") {
      RiskGraphs <- InputsForOutputs[[ModelType]]$FEVD$WishGraphsOrtho$RiskFactors
      YieldGraphs <- InputsForOutputs[[ModelType]]$FEVD$WishGraphsOrtho$Yields
    } else { # GIRF ortho
      RiskGraphs <- InputsForOutputs[[ModelType]]$GFEVD$WishGraphsOrtho$RiskFactors
      YieldGraphs <- InputsForOutputs[[ModelType]]$GFEVD$WishGraphsOrtho$Yields
    }
  } else {
    # 2) All other models
    if (OutType == "FEVD") {
      RiskGraphs <- InputsForOutputs[[ModelType]]$FEVD$WishGraphs$RiskFactors
      YieldGraphs <- InputsForOutputs[[ModelType]]$FEVD$WishGraphs$Yields
    } else if (OutType == "GFEVD") {
      RiskGraphs <- InputsForOutputs[[ModelType]]$GFEVD$WishGraphs$RiskFactors
      YieldGraphs <- InputsForOutputs[[ModelType]]$GFEVD$WishGraphs$Yields
    }
  }
  return(list(RiskGraphs = RiskGraphs, YieldGraphs = YieldGraphs))
}

###########################################################################################################
#' Limit the number of categories in FEVDs and GFEVDs graphs by merging small groups into other
#'
#' @param FEVDlist list of FEVD and GFEVD outputs
#' @param Threshold Threshold for merging factors. Default is 0.05.
#'
#' @keywords internal

MergeFEVD_graphs <- function(FEVDlist, Threshold = 0.05) {
  variable <- FEVDlist$variable

  mean_contrib <- tapply(FEVDlist$value, FEVDlist$variable, mean)
  small_vars <- names(mean_contrib[mean_contrib < Threshold])
  FEVD_cleaned <- transform(FEVDlist, variable = ifelse(variable %in% small_vars, "Others", as.character(variable)))
  FEVD_cleaned <- stats::aggregate(value ~ index + variable, data = FEVD_cleaned, sum)

  return(FEVD_cleaned)
}
###########################################################################################################
#' Modify variable labels to make legends more readable
#'
#' @param x variable labels prior to the modification
#' @keywords internal

clean_labels <- function(x) {
  # Replace underscores and dots with spaces or dashes
  x <- gsub("_", " ", x)
  x <- gsub("\\.", " - ", x)

  # Capitalize first letter of each word
  x <- sapply(strsplit(x, " "), function(words) {
    paste(toupper(substring(words, 1, 1)), substring(words, 2), sep = "", collapse = " ")
  })

  return(x)
}
