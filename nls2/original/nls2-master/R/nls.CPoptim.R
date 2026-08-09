nls.CPoptim <- function(formula, data = parent.frame(), start, 
	control = nls.control(),
	trace = FALSE, weights, subset, ...) { 

   qdata <- substitute(data)

   stopifnot(NROW(start) == 2L)
   attr(data, "reference") <- NULL
   f <- function(st) {
	   names(st) <- names(start)
	   mod <- nlsModel(formula, data, st)
	   mod$deviance()
    }

    if (is.null(control$FEmax)) control$FEmax <- 5000 * NCOL(start)
    if (is.null(control$sampleSize)) control$sampleSize <- 1000
    
    res <- CPoptim::CPoptim(f, unlist(start[1, ]), unlist(start[2, ]), control$FEmax,
      control$sampleSize)
    co <- res$sample$x[which.min(res$sample$y), ]
    names(co) <- names(start)
    m <- nlsModel(formula, data, co)
    call <- match.call()
    structure(list(m = m, call = call, data = qdata), class = "nls")

}

# library(CPoptim)
# out <- nls.CPoptim(demand ~ a + b * Time, data = BOD, start = 
#  data.frame(a = as.numeric(1:2), b = as.numeric(1:2)))
