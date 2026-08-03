#' Computes the country-specific spanned factors
#'
#' @param Yields matrix (\code{J x Td}). Bond yields for all countries.
#' @param Economies character vector. Names of the \code{C} economies included in the system.
#' @param N integer. Desired number of country-specific spanned factors (maximum allowed is \code{N = J}).
#'
#' @return matrix. Contains the \code{N} spanned factors for all countries in the system (\code{CJ x Td}).
#'
#' @section General Notation:
#' \itemize{
#'   \item \code{Td}: model time series dimension
#'   \item \code{C}: number of countries in the system
#'   \item \code{N}: number of country-specific spanned factors
#'   \item \code{J}: number of bond yields per country used in estimation
#' }
#'
#' @examples
#' data(Yields)
#' Economies <- c("China", "Brazil", "Mexico", "Uruguay")
#' N <- 3
#' SpaFact_TS <- Spanned_Factors(Yields, Economies, N)
#'
#' @export


Spanned_Factors <- function(Yields, Economies, N) {
  # Validate inputs
  if (!is.matrix(Yields)) {
    stop("Yields must be a matrix.")
  }
  if (!is.character(Economies) || length(Economies) == 0) {
    stop("Economies must be a non-empty character vector.")
  }
  if (!is.numeric(N) || N <= 0 || N != round(N)) {
    stop("N must be a positive integer.")
  }

  C <- length(Economies) # number of the economies of the system

  # Label of the spanned factors
  LabelsSpannedALL <- c("Level", "Slope", "Curvature", "Fourth PC", "Fifth PC", "Sixth PC", "Seventh PC", "Eighth PC")
  LabelsSpanned <- LabelsSpannedALL[1:N]

  # Spanned factors:
  P <- c()

  for (i in 1:C) {
    Idx <- grepl(Economies[i], rownames(Yields)) # Extract the yields of a single country
    Y_CS <- Yields[Idx, ] # Country-specific yields
    W <- pca_weights_one_country(Y_CS, Economies[i]) # Weight matrix
    if (N > nrow(Y_CS)) {
      stop("The Number of country-specific spanned factors cannot exceed the number of country-specific bond yields.")
    }
    W <- W[1:N, ] * 100 # Select only the set of weights that will be used to compute the first N PCs
    if (i == 1) {
      P <- W %*% Y_CS # Country-specific spanned factor
      row.names(P) <- paste(LabelsSpanned, Economies[i])
    } else {
      Ptemp <- W %*% Y_CS # Country-specific spanned factor
      row.names(Ptemp) <- paste(LabelsSpanned, Economies[i])
      P <- rbind(P, Ptemp)
    }
  }

  return(P)
}
