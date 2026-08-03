#' Compute the cross-section loadings of yields of a canonical A0_N model
#'
#' @param mat       vector of maturities (J x 1).
#' @param K1XQ      risk-neutral feedback matrix (N x N)
#' @param ModelType string-vector containing the label of the model to be estimated
#' @param r0        the long run risk neutral (scalar)
#' @param SSX       covariance matrix of the latent states (N x N)
#' @param Economies string-vector containing the names of the economies which are part of the economic system
#'
#'
#' @references
#'  \itemize{
#'  \item Dai, Q. and Singleton, K. (2000). "Specification Analysis of Affine Term Structure Models". The Journal of Finance.
#'  \item Joslin, S., Singleton, K. and Zhu, H. (2011). "A new perspective on Gaussian dynamic term structure models".
#'  The Review of Financial Studies.
#'  }
#'
#' @keywords internal

Get__BnXAnX <- function(mat, K1XQ, ModelType, r0 = NULL, SSX = NULL, Economies) {
  # Preliminary checks - Input validations
  if (!is.matrix(K1XQ) && !is.numeric(K1XQ)) {
    stop("'K1XQ' must be either a numeric or a matrix.")
  }

  if (!is.null(SSX) && (!is.matrix(SSX) || nrow(SSX) != ncol(SSX))) {
    stop("'SSX' must be a square matrix if provided.")
  }

  if (!is.null(r0) && !is.numeric(r0)) {
    stop("'r0' must be numeric or NULL.")
  }

  # 0) Prelimiary work
  Lab_SingleQ <- c("JPS original", "JPS global", "GVAR single")

  K <- max(mat) # longest maturity of interest
  J <- length(mat)
  if (any(ModelType == Lab_SingleQ)) {
    N <- if (is.matrix(K1XQ)) nrow(K1XQ) else 1 # Number of latent factors
  } else {
    C <- length(Economies)
    NC <- dim(K1XQ)[1] # Total number of latent factors of the system
    N <- NC / C # number of spanned factors per country
  }

  if (N %% 1 != 0) stop("Number of factors per country N is not an integer.")

  # 1) Compute dX
  # a) For models estimate separetely
  if (any(ModelType == Lab_SingleQ)) {
    dX <- rep(1, times = N)
    dX <- c(dX)

    # b) For models estimated jointly
  } else {
    dxCS <- rep(1, N)
    dX <- matrix(0, nrow = C, ncol = NC)
    idx0 <- 0

    for (j in 1:C) {
      idx1 <- idx0 + N
      dX[j, (idx0 + 1):idx1] <- dxCS
      idx0 <- idx1
    }
  }

  # 2) Compute BnX and AnX
  LoadX <- Compute_BnX_AnX(mat, N, K1XQ, r0, dX, SSX, Economies, ModelType, Lab_SingleQ)

  # 3) Compute the adjstment term B_adj, since: AnX = r0 + B_adj
  # For models estimate separetely
  if (any(ModelType == Lab_SingleQ)) {
    if (is.null(r0)) {
      B_adj <- c()
    } else {
      B_adj <- LoadX$AnX - rep(r0, times = J)
    }
    # For models estimated jointly
  } else {
    idx0 <- 0
    if (is.null(r0)) {
      B_adj <- as.numeric(vector())
    } else {
      B_adj <- c()
      for (h in 1:C) {
        idx1 <- idx0 + J
        B_adj[(idx0 + 1):idx1] <- LoadX$AnX[(idx0 + 1):idx1] - t(t(rep(r0[h], times = J)))
        idx0 <- idx1
      }
    }
  }

  return(list(BnX = LoadX$BnX, AnX = LoadX$AnX, B_adj = B_adj))
}


##########################################################################################################
#' Compute the latent loading AnX and BnX
#'
#' @param mat vector of maturities (J x 1). Maturities are in multiples of the discrete interval used in the model
#' @param N  number of country-specific spanned factors
#' @param K1XQ risk neutral feedback matrix (N x N)
#' @param r0 the long run risk neutral mean of the short rate (scalar)
#' @param dX state loadings for the one-period rate (1xN). Default is a vector of ones
#' @param SSX the covariance matrix of the errors (N x N)
#' @param Economies string-vector containing the names of the economies which are part of the economic system
#' @param ModelType string-vector containing the label of the model to be estimated
#' @param Lab_SingleQ string-vector containing the labels of the models estimated on a country-by-country basis
#'
#' @keywords internal

Compute_BnX_AnX <- function(mat, N, K1XQ, r0, dX, SSX, Economies, ModelType, Lab_SingleQ) {
  K <- max(mat)
  # For models estimate separately
  if (any(ModelType == Lab_SingleQ)) {
    BnX <- matrix(0, nrow = K + 1, ncol = N)
    AnX <- matrix(0, nrow = K + 1, ncol = 1)

    # Generate the loadings from a canonical A_0(N) model:
    for (k in seq_len(K)) {
      BnX[k + 1, ] <- dX + BnX[k, ] %*% K1XQ # Note: BnX depends only on K1XQ and dx NOT r0 or SSX.
      if (is.null(r0)) {
        AnX <- as.numeric(vector())
      } else {
        AnX[k + 1] <- r0 + AnX[k] - 0.5 * BnX[k, ] %*% SSX %*% t(t(BnX[k, ]))
      }
    }

    # Adjust the loading for the maturity
    BnX <- BnX[mat + 1, ] / mat
    AnX <- AnX[mat + 1] / mat

    # For models estimated jointly
  } else {
    C <- length(Economies)
    idx0 <- 0
    for (i in 1:C) {
      AnXCS <- matrix(0, nrow = K + 1, ncol = 1) # Country-specific A
      BnXCS <- matrix(0, nrow = K + 1, ncol = N * C) # Country-specific B
      for (k in 1:K) {
        idx1 <- idx0 + N
        BnXCS[k + 1, ] <- dX[i, ] + BnXCS[k, ] %*% K1XQ # Note: BnX depends only on K1XQ and dx NOT r0 or SSX.
        if (is.null(r0)) {
          AnXCS <- as.numeric(vector())
        } else {
          AnXCS[k + 1] <- r0[i] + AnXCS[k] - 0.5 * BnXCS[k, ] %*% SSX %*% t(t(BnXCS[k, ]))
        }
      }
      if (i == 1) {
        BnX <- BnXCS[mat + 1, (idx0 + 1):idx1] / mat # Adjust the loading for the maturity
        if (N == 1) {
          BnX <- matrix(BnX)
        }

        AnXCS <- AnXCS[mat + 1] / mat
        AnX <- t(t(AnXCS))

        idx0 <- idx1
      } else {
        BnXCS <- BnXCS[mat + 1, (idx0 + 1):idx1] / mat

        if (N == 1) {
          BnXCS <- matrix(BnXCS)
        }
        BnX <- adiag(BnX, BnXCS)

        AnXCS <- t(t(AnXCS[mat + 1] / mat))
        AnX <- rbind(AnX, AnXCS)

        idx0 <- idx1
      }
    }
  }


  return(list(AnX = AnX, BnX = BnX))
}
