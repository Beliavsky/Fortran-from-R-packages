maximMessage <- function(code) {
   message <- switch(code,
                     "1" = "gradient close to zero (gradtol)",
                     "2" = "successive function values within tolerance limit (tol)",
                     "3" = paste("Last step could not find a value above the",
                     "current.\nBoundary of parameter space?",
                     " \nConsider switching to a more robust optimisation method temporarily."),
                     "4" = "Iteration limit exceeded (iterlim)",
                     "5" = "Infinite value",                        
                     "6" = "Infinite gradient",                        
                     "7" = "Infinite Hessian",
                     "8" = "successive function values within relative tolerance limit (reltol)",
                     "9" = paste("Gradient did not change,",
                                 "cannot improve BFGS approximation for the Hessian.\n",
                                 "Use different optimizer and/or analytic gradient."),
                     "10" = "Lost patience (SG_patience)",
                     "100" = "Initial value out of range.",
                     paste("Code", code))
   return(message)
}
