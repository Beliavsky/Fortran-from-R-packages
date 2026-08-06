brick.wall.2d <- function (x, method = "modwt") {
  wf <- attributes(x)$wavelet
  m <- wave.filter(wf)$length
  for (i in names(x)) {
    j <- as.numeric(substr(i, 3, 3))
    if (method == "dwt") {
      n <- ceiling((m - 2) * (1 - 1/2^j))
    } else {
      n <- (2^j - 1) * (m - 1)
    }
    n <- min(n, nrow(x[[i]]))
    x[[i]][1:n, ] <- NA
    x[[i]][, 1:n] <- NA
  }
  
  return(x)
}

#' Wavelet Analysis of Images
#' 
#' Produces an estimate of the multiscale variance with approximate 
#' confidence intervals using the 2D MODWT.
#' 
#' The wavelet variance is basically the average of the squared wavelet 
#' coefficients across each scale and direction of an image. As shown 
#' in Mondal and Percival (2012), the wavelet variance is a 
#' scale-by-scale decomposition of the variance for a stationary spatial
#' process, and certain non-stationary spatial processes.
#' 
#' @param x image
#' @param p (one minus the) two-sided p-value for the confidence interval
#' @return Data frame with 3J+1 rows.
#' @author B. Whitcher
#' @references Mondal, D. and D. B. Percival (2012). Wavelet variance 
#' analysis for random fields on a regular lattice. \emph{IEEE 
#' Transactions on Image Processing} \bold{21}, 537–549.
#'
wave.variance.2d <- function(x, p = 0.025) {
  
  # The unbiased estimator ignores those coefficients affected by the boundary
  x_bw <- brick.wall.2d(x)
  x_ss <- unlist(lapply(x_bw, FUN = function(v) sum(v * v, na.rm = TRUE)))
  x_length <- unlist(lapply(x_bw, FUN = function(v) sum(! is.na(v))))
  wave_var <- x_ss / x_length
  
  edof <- rep(NA, length(x))
  names(edof) <- names(x)
  wf_length <- wave.filter(attributes(x)$wavelet)$length
  # from Section 3.3 Confidence intervals in Geilhufe et al. (2013)
  for (i in names(x)) {
    j <- as.integer(substr(i, 3, 3))
    Lj <- (2^j - 1) * (wf_length - 1) + 1
    Nj <- nrow(x[[i]]) - Lj + 1
    Mj <- ncol(x[[i]]) - Lj + 1
    pad_with_zeros <- matrix(0, nrow = 2^(trunc(log2(Nj)) + 2), ncol = 2^(trunc(log2(Mj)) + 2))
    pad_with_zeros[1:Nj, 1:Mj] <- x[[i]][Lj:nrow(x[[i]]), Lj:ncol(x[[i]])]
    sW <- fft(fft(pad_with_zeros) * Conj(fft(pad_with_zeros)), inverse = TRUE) / prod(dim(pad_with_zeros)) / Nj / Mj
    sigma_W <- sum(sW^2)
    if (Nj * Mj > 128) {
      edof[i] <- 2 * Nj * Mj * wave_var[i]^2 / Re(sigma_W)
    } else {
      edof[i] <- max((Nj * Mj) / (2^j * 2^j), 1)
    }
  }
  
  data.frame(
    value = wave_var,
    level = as.integer(substr(names(x), 3, 3)),
    direction = factor(
      toupper(substr(names(x), 1, 2)),
      levels = c("LH", "HL", "HH", "LL"),
      labels = c("Horizontal", "Vertical", "Diagonal", "Approximation")
    ),
    ci_lower = unlist(edof) * wave_var / qchisq(1 - p, unlist(edof)),
    ci_upper = unlist(edof) * wave_var / qchisq(p, unlist(edof))
  )
}
