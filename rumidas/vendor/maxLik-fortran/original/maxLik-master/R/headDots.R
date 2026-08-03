### paste head of vector, and if some of it is left out, add '...' to it.

headDots <- function(x, max.cols) {
   s <- paste(head(x, max.cols), collapse=", ")
   if(length(x) > max.cols) {
      s <- paste(s, "...")
   }
   s
}
