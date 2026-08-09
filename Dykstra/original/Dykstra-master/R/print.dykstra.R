print.dykstra <-
  function(x, n = 5L, ...){
    
    cat("converged: ", x$converged, " (", x$iter, " iterations)", sep = "")
    if(length(x$solution) <= n) {
      cat("\n solution:", round(x$solution, getOption("digits"))) 
    } else {
      cat("\n solution:", "vector of length", length(x$solution))
    }
    cat("\n    value:", round(x$value, getOption("digits")),"\n")
  }