#' @title Joint Loss Function
#' @description Computes the joint (VaR, ES) loss
#' @param r Vector of returns
#' @param q Vector of quantiles
#' @param e Vector of expected shortfalls
#' @param alpha Probability level
#' @param g1 1, 2, see \link{G1_fun}
#' @param g2 1, 2, 3, 4, 5, see \link{G2_curly_fun}, \link{G2_fun}
#' @param return_mean If TRUE returns the average tick loss, else the individual values
#' @references Fissler and Ziegel (2016)
#' @export
esr_loss <- function(r, q, e, alpha, g1 = 2L, g2 = 1L, return_mean = TRUE) {
  G1_q <- G_vec(z = q, g = "G1", type = g1)
  G1_r <- G_vec(z = r, g = "G1", type = g1)
  G2_curly_e <- G_vec(z = e, g = "G2_curly", type = g2)
  G2_e <- G_vec(z = e, g = "G2", type = g2)

  loss <- ((r <= q) - alpha) * G1_q - (r <= q) * G1_r +
    G2_e * (e - q + (q - r) * (r <= q)/alpha) - G2_curly_e

  if (return_mean) {
    mean(loss)
  } else {
    loss
  }
}
