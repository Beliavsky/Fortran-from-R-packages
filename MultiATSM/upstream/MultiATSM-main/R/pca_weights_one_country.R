#' Computes the PCA weights for a single country
#'
#' @param Yields matrix (\code{J x Td}). Bond yields for a single country.
#' @param Economy character. Name of the economy.
#'
#' @return matrix (\code{J x J}). Eigenvectors of the variance-covariance matrix of yields.
#'
#' @section General Notation:
#' \itemize{
#'   \item \code{Td}: model time series dimension
#'   \item \code{J}: number of bond yields per country used in estimation
#' }
#'
#' @examples
#' data(Yields)
#' Economy <- "Mexico"
#' pca_weights <- pca_weights_one_country(Yields, Economy)
#'
#' @export


pca_weights_one_country <- function(Yields, Economy) {
  # Extract the yields of a single country
  Idx <- grepl(Economy, rownames(Yields))
  Y_CS <- Yields[Idx, ] # Country-specific yields

  # Store the eigenvectors of the matrix of yields
  H <- eigen(stats::cov(t(Y_CS)))
  W <- t(H$vectors)

  J <- nrow(Y_CS) # Total number of maturities
  M <- round(stats::median(1:J)) # Median of the maturity range


  # Adjust the weights to ease the interpretability of the spanned factors
  # a) LEVEL: if the weights of the level are negative, then invert the sign
  # (ideally, high values of level should be interpreted as a high values for all yields)
  if (all(W[1, ] < 0)) {
    W[1, ] <- W[1, ] * (-1)
  }
  # b) SLOPE: if the slope is downward sloping, then invert the sign
  # (ideally, high values of slope should be interpreted as a steep high curve)
  if (W[2, 1] > W[2, J]) {
    W[2, ] <- W[2, ] * (-1)
  }
  # c) Curvature: if the curvature is U-shapped, then invert the sign
  # (ideally, the curvature factor should capture an effect that is linked to mid-term maturities)
  if (W[3, 1] > W[3, M] & W[3, J] > W[3, M]) {
    W[3, ] <- W[3, ] * (-1)
  }


  return(W)
}
