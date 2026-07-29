msmoothCalc <- function (y, p = c(1, 3), mu = c(0, 1, 2, 3), bStart = 0.15,
                     alg = c("A", "B", "N", "NA", "NAM", "NM", "O", "OA",
                             "OAM", "OM"),
                     method = c("lpr", "kr")) {

  result <- tsmooth.Alg(y, p, mu, bStart, alg, method)

  result
}
