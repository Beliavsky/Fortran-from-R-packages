#'
#'@importFrom grDevices col2rgb rgb
#'

col.alpha <- function(col, alpha = 1) {
  col.rgb <- col2rgb(col) / 255
  red <- col.rgb[1]
  green <- col.rgb[2]
  blue <- col.rgb[3]
  outp <- rgb(red = red, green = green, blue = blue, alpha = alpha)
  outp
}
