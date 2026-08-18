histDens <- function( x, breaks = "Sturges", ... ) {
  h <- hist( x, breaks = breaks, plot = FALSE )
  d <- density( x )
  hist( x, breaks = breaks, prob = TRUE,
    ylim = c( 0, max( c( h$density, d$y ) ) ), ... )
  lines( d, lwd = 2 )
  invisible( NULL )
}
