#' Generates the labels for risk factors used in the model
#'
#' @param N positive integer. Number of country-specific spanned factors. Must be between 1 and 8.
#' @param DomVar character vector. Names of the domestic variables.
#' @param GlobalVar character vector. Names of the global variables.
#' @param Economies character vector. Names of the economies included in the system.
#' @param ModelType character. Model type to be estimated. Permissible choices: "JPS original", "JPS global", "GVAR single", "JPS multi", "GVAR multi", "JLL original", "JLL No DomUnit", "JLL joint Sigma".
#'
#' @return List containing the risk factor labels for spanned, domestic, star, and global variables, as well as tables for each country and all countries.
#'
#' @examples
#' N <- 2
#' DomVar <- c("inflation", "Output gap")
#' GlobalVar <- "Commodity Prices"
#' Economies <- c("U.S.", "Canada", "Germany", "Japan")
#' ModelType <- "JPS original"
#' VarLabels <- LabFac(N, DomVar, GlobalVar, Economies, ModelType)
#' @export

LabFac <- function(N, DomVar, GlobalVar, Economies, ModelType) {
  M <- length(DomVar) # Number of domestic unspanned factors
  C <- length(Economies) # Number of countries of the economic system

  FactorLabels <- list()
  FactorLabels$Spanned <- LabelsSpanned(N)
  FactorLabels$Domestic <- c(DomVar, FactorLabels$Spanned)
  FactorLabels$Star <- LabelsStar(FactorLabels$Domestic)
  FactorLabels$Global <- GlobalVar

  labelsDomVar <- c()
  labelsDomVarJLL <- c()
  for (i in seq_len(C)) {
    for (j in seq_len(N + M)) {
      labelsDomVar <- c(labelsDomVar, paste(FactorLabels$Domestic[j], Economies[i]))
      if (any(ModelType == c("JLL original", "JLL No DomUnit", "JLL joint Sigma"))) {
        labelsDomVarJLL <- c(labelsDomVarJLL, paste(FactorLabels$Domestic[j], Economies[i], "JLL"))
      }
    }
    FactorLabels$Tables[[Economies[i]]] <- labelsDomVar[(length(labelsDomVar) - (N + M) + 1):length(labelsDomVar)]
    FactorLabels$TablesJLL[[Economies[i]]] <- labelsDomVarJLL[(length(labelsDomVarJLL) - (N + M) + 1):length(labelsDomVarJLL)]
  }

  idx0 <- 0
  for (a in seq_len(C)) {
    idx1 <- idx0 + N + M
    FactorLabels$Tables$AllCountries[(idx0 + 1):idx1] <- do.call(cbind, lapply(FactorLabels$Tables[[Economies[a]]], matrix, nrow = 1))
    if (any(ModelType == c("JLL original", "JLL No DomUnit", "JLL joint Sigma"))) {
      FactorLabels$Tables$AllCountriesJLL[(idx0 + 1):idx1] <- do.call(cbind, lapply(FactorLabels$TablesJLL[[Economies[a]]], matrix, nrow = 1))
    }
    idx0 <- idx1
  }

  FactorLabels$TablesJLL <- NULL # Remove the country-specific JLL factors list

  return(FactorLabels)
}

########################################################################################################
#' Generate the labels of the spanned factors
#'
#' @param N number of spanned factors
#'
#' @keywords internal

LabelsSpanned <- function(N) {
  if (!(N %in% 1:8)) {
    stop("N, the number of country-specific spanned factors, must be an integer between 1 and 8.")
  }

  LabelsSpannedALL <- c("Level", "Slope", "Curvature", "Fourth PC", "Fifth PC", "Sixth PC", "Seventh PC", "Eighth PC")
  LabelsSpan <- LabelsSpannedALL[1:N]

  return(LabelsSpan)
}

###################################################################################################
#' Generate the labels of the star variables
#'
#' @param FactorLabels Factor labels
#'
#' @keywords internal

LabelsStar <- function(FactorLabels) {
  Label <- paste(FactorLabels, ".Star", sep = "")
  return(Label)
}
