#'
#'@importFrom grDevices col2rgb colors
#'

color.name <- function(col) {
  col.rgb <- c(col2rgb(col))
  colors.all.rgb <- col2rgb(colors())
  eucl.d <- sqrt((col.rgb[[1]] - colors.all.rgb[1, ])^2 +
    (col.rgb[[2]] - colors.all.rgb[2, ])^2 +
    (col.rgb[[3]] - colors.all.rgb[3, ])^2)
  col.ind <- which(eucl.d == min(eucl.d))[[1]]
  color.out <- gsub("[0-9]+", "", colors()[[col.ind]])
  color.out
}
