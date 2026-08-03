#' Collects the inputs that are used to construct the numerical and graphical outputs
#'
#' @param ModelType character. Model type to be estimated. Permissible choices: "JPS original", "JPS global", "GVAR single", "JPS multi", "GVAR multi", "JLL original", "JLL No DomUnit", "JLL joint Sigma".
#' @param Horiz numeric scalar. Desired analysis horizon for the outputs.
#' @param ListOutputWished character vector. Desired graphical outputs. Available options: "RiskFactors", "Fit", "IRF", "FEVD", "GIRF", "GFEVD", "TermPremia", "ForwardPremia".
#' @param OutputLabel character. Name of the output label to be stored.
#' @param WishStationarityQ logical. Whether to impose that the largest eigenvalue under Q is strictly smaller than 1. TRUE to impose.
#' @param DataFrequency character. Data frequency. Permissible choices: "Daily All Days", "Daily Business Days", "Weekly", "Monthly", "Quarterly", "Annually".
#' @param WishGraphYields logical. Whether to generate graphs for yields. Default is FALSE.
#' @param WishGraphRiskFactors logical. Whether to generate graphs for risk factors. Default is FALSE.
#' @param WishOrthoJLLgraphs logical. Whether to generate orthogonalized JLL-based graphs. Default is FALSE.
#' @param WishForwardPremia logical. Whether to generate forward premia graphs. Default is FALSE.
#' @param LimFP numeric vector. Maturities associated with the start and end dates of the loan.
#' @param WishBootstrap logical. Whether to perform bootstrap-based estimation. Default is FALSE.
#' @param ListBoot list. Contains bootstrap settings: methodBS ("bs", "wild", "block"), BlockLength (numeric), ndraws (numeric), pctg (numeric).
#' @param WishForecast logical. Whether to generate forecasts. Default is FALSE.
#' @param ListForecast list. Contains forecast settings: ForHoriz (numeric), t0Sample (numeric), t0Forecast (numeric), ForType ("Rolling", "Expanding").
#' @param UnitYields character. Maturity unit of yields. Options: "Month" or "Year". Default is "Month".
#'
#' @return List of necessary inputs to generate the graphs and outputs of the desired model.
#'
#' @examples
#'
#' ModelType <- "JPS original"
#' Horiz <- 100
#' DesiredOutputGraphs <- c("Fit", "GIRF", "GFEVD")
#' OutputLabel <- "Test"
#' WishStationarityQ <- TRUE
#' WishGraphRiskFac <- FALSE
#' WishGraphYields <- TRUE
#'
#' InputsList <- InputsForOutputs(
#'   ModelType, Horiz, DesiredOutputGraphs, OutputLabel,
#'   WishStationarityQ, WishGraphYields, WishGraphRiskFac
#' )
#'
#' @export

InputsForOutputs <- function(ModelType, Horiz, ListOutputWished, OutputLabel, WishStationarityQ, DataFrequency,
                             WishGraphYields = FALSE, WishGraphRiskFactors = FALSE, WishOrthoJLLgraphs = FALSE,
                             WishForwardPremia = FALSE, LimFP = NULL, WishBootstrap = FALSE, ListBoot = NULL,
                             WishForecast = FALSE, ListForecast = NULL, UnitYields = "Month") {
  OutputTypeSet <- c("RiskFactors", "Fit", "IRF", "FEVD", "GIRF", "GFEVD", "TermPremia", "ForwardPremia")
  IdxWishOut <- match(ListOutputWished, OutputTypeSet)

  if ("RiskFactors" %in% ListOutputWished) {
    WishRFGraph <- TRUE
  }
  if ("Fit" %in% ListOutputWished) {
    WishFitGraph <- TRUE
  }
  if ("TermPremia" %in% ListOutputWished) {
    WishTPGraph <- TRUE
  }

  InputsForOutputs <- list(
    "Label Outputs" = OutputLabel,
    StationaryQ = WishStationarityQ,
    UnitMatYields = UnitYields,
    DataFreq = DataFrequency,
    ForwardPremia = WishForwardPremia
  )

  # Initialization
  for (h in 1:length(OutputTypeSet)) {
    if (h == 1) {
      InputsForOutputs[[ModelType]]$RiskFactors$WishGraphs <- FALSE
    }
    if (h == 2) {
      InputsForOutputs[[ModelType]]$Fit$WishGraphs <- FALSE
    }
    if (h >= 3 & h <= 6) {
      InputsForOutputs[[ModelType]][[OutputTypeSet[h]]]$horiz <- Horiz

      InputsForOutputs[[ModelType]][[OutputTypeSet[h]]]$WishGraphs$RiskFactors <- FALSE
      InputsForOutputs[[ModelType]][[OutputTypeSet[h]]]$WishGraphs$Yields <- FALSE
      InputsForOutputs[[ModelType]][[OutputTypeSet[h]]]$WishGraphs$RiskFactorsBootstrap <- FALSE
      InputsForOutputs[[ModelType]][[OutputTypeSet[h]]]$WishGraphs$YieldsBootstrap <- FALSE

      if (any(ModelType == c("JLL original", "JLL No DomUnit", "JLL joint Sigma"))) {
        InputsForOutputs[[ModelType]][[OutputTypeSet[h]]]$WishGraphsOrtho$RiskFactors <- FALSE
        InputsForOutputs[[ModelType]][[OutputTypeSet[h]]]$WishGraphsOrtho$Yields <- FALSE
        InputsForOutputs[[ModelType]][[OutputTypeSet[h]]]$WishGraphsOrtho$RiskFactorsBootstrap <- FALSE
        InputsForOutputs[[ModelType]][[OutputTypeSet[h]]]$WishGraphsOrtho$YieldsBootstrap <- FALSE
      }
    }
    if (h > 6) {
      InputsForOutputs[[ModelType]]$RiskPremia$WishGraphs <- FALSE
      InputsForOutputs[[ModelType]]$ForwardPremia$Limits <- LimFP
    }
  }


  # Filling up the desired outputs
  for (h in IdxWishOut) {
    if (h == 1) {
      InputsForOutputs[[ModelType]]$RiskFactors$WishGraphs <- WishRFGraph
    }
    if (h == 2) {
      InputsForOutputs[[ModelType]]$Fit$WishGraphs <- WishFitGraph
    }

    if (h >= 3 & h <= 6) {
      InputsForOutputs[[ModelType]][[OutputTypeSet[h]]]$horiz <- Horiz

      InputsForOutputs[[ModelType]][[OutputTypeSet[h]]]$WishGraphs$RiskFactors <- WishGraphRiskFactors
      InputsForOutputs[[ModelType]][[OutputTypeSet[h]]]$WishGraphs$Yields <- WishGraphYields
      if (WishBootstrap == 1) {
        InputsForOutputs[[ModelType]][[OutputTypeSet[h]]]$WishGraphs$RiskFactorsBootstrap <- WishGraphRiskFactors
        InputsForOutputs[[ModelType]][[OutputTypeSet[h]]]$WishGraphs$YieldsBootstrap <- WishGraphYields
      }

      if (any(ModelType == c("JLL original", "JLL No DomUnit", "JLL joint Sigma"))) {
        if (WishOrthoJLLgraphs) {
          InputsForOutputs[[ModelType]][[OutputTypeSet[h]]]$WishGraphsOrtho$RiskFactors <- WishGraphRiskFactors
          InputsForOutputs[[ModelType]][[OutputTypeSet[h]]]$WishGraphsOrtho$Yields <- WishGraphYields
          if (WishBootstrap) {
            InputsForOutputs[[ModelType]][[OutputTypeSet[h]]]$WishGraphsOrtho$RiskFactorsBootstrap <- WishGraphRiskFactors
            InputsForOutputs[[ModelType]][[OutputTypeSet[h]]]$WishGraphsOrtho$YieldsBootstrap <- WishGraphYields
          }
        }
      }
    }

    if (h == 7) {
      InputsForOutputs[[ModelType]]$RiskPremia$WishGraphs <- WishTPGraph
    }
  }

  # 2) Bootstrap list
  if (!WishBootstrap) {
    InputsForOutputs[[ModelType]]$Bootstrap$WishBoot <- FALSE
  } else {
    InputsForOutputs[[ModelType]]$Bootstrap <- WishBootstrap
    names(InputsForOutputs[[ModelType]]$Bootstrap) <- "WishBoot"
    InputsForOutputs[[ModelType]]$Bootstrap <- append(InputsForOutputs[[ModelType]]$Bootstrap, ListBoot)
  }


  # 3) Forecasting list
  if (!WishForecast) {
    InputsForOutputs[[ModelType]]$Forecasting$WishForecast <- FALSE
  } else {
    InputsForOutputs[[ModelType]]$Forecasting <- WishForecast
    names(InputsForOutputs[[ModelType]]$Forecasting) <- "WishForecast"
    InputsForOutputs[[ModelType]]$Forecasting <- append(InputsForOutputs[[ModelType]]$Forecasting, ListForecast)
  }


  return(InputsForOutputs)
}
