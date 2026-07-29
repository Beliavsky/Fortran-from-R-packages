tsmooth.Alg <- function(y, p, mu, bStart, alg, method) {
  results <- tsmoothCalc(y, p, mu, lookup$lookup_alg["Mcf", alg],
                     lookup$lookup_alg["InfR", alg],
                     bStart, lookup$lookup_alg["bvc", alg], 1, 0.05, method)
  results
}
