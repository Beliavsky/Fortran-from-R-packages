.G_function_names <- function(g1, g2) {
  if (g1 == 1) {
    g1_name <- "G1(z) = z"
  } else if (g1 == 2) {
    g1_name <- "G1(z) = 0"
  }
  if (g2 == 1) {
    g2_name <- "G2(z) = -1/z, z < 0"
  } else if (g2 == 2) {
    g2_name <- "G2(z) = 0.5 / sqrt(z), z < 0"
  } else if (g2 == 3) {
    g2_name <- "G2(z) = 1 / z^2, z < 0"
  } else if (g2 == 4) {
    g2_name <- "G2(z) = 1 / (1 + exp(-z))"
  } else if (g2 == 5) {
    g2_name <- "G2(z) = exp(z)"
  }
  c(g1_name, g2_name)
}
