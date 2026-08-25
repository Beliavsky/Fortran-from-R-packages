# Miscellaneous functions of R package mnormt
# Author: Adelchi Azzalini
# 
sample_Mardia_measures <- function(data, correct=FALSE)
{# Computes measures of Mardia (1970, Biometrika, vol.57, pp.519-530) and
 # (1974, Sankhya ser.B, vol.36, pp.115-128 - includes equations quoted below).
 # R function written in 2020.
  y <- data.matrix(data)
  y <- y[complete.cases(y), , drop=FALSE]
  n <- nrow(y)
  if(n < 4) stop("condition n>3 is required")
  d <- ncol(y)
  m <- if(d>1) apply(y, 2, mean) else mean(y)
  y0 <- y - outer(rep(1, n), m)
  t.y0 <- t(y0)
  f <- if(correct) 1 else (n-1)/n 
  S <-  var(y0)*f  # f=(n-1)/n gives results as on p.127 of Mardia (1974)
  Sinv <- pd.solve(S)
  b1 <- mean((y0 %*% Sinv %*% t.y0)^3)                       # eq.(2.2)
  b2 <- mean(colSums((Sinv %*%  t.y0) * t.y0)^2)             # eq.(2.4)
  g1 <- n*(n-1)*b1/(n-2)^2                                   # eq.(2.10)
  g2 <- (n-1)*((n+1)*b2 - (n-1)*d*(d+2))/((n-2)*(n-3))       # eq.(2.11)
  k <- (d+1)*(n+1)*(n+3)/(n*((n+1)*(d+1)-6))                 # eq.(5.5)
  p.b1 <- 1 - pchisq(b1*n*k/6, df=d*(d+1)*(d+2)/6)
  p.b2 <- if((n-d-1) > 0) {                                  # next 3 lines \
    r <- sqrt((n+3)*(n+5)/(n-d-1))                           # use eq.(5.6)  
    std.b2 <- ((n+1)*b2-d*(d+2)*(n-1))*r/((8*d*(d+2)*(n-3)*(n-d-1)))
    p.b2 <- 2*pnorm(-abs(std.b2))
    } else  NA
  return(c(b1=b1, b2=b2, g1=g1, g2=g2, p.b1=p.b1, p.b2=p.b2, n=n))
}
#
plot_fxy <- function(f, xlim, ylim, ..., npt=51, grf, grpar) {
  # 2022-06-01 
  dots <- list(...)
  dnm <- names(dots)    # 'dots' names
  funct <- if(is.character(f)) get(f, inherits=TRUE) else f
  ff <- formals(funct) 
  fnm <- names(ff[-1])  # 'formal' names
  nmm <- match(dnm, fnm)
  if(any(is.na(nmm))) stop("some arguments do not match those of 'f'")
  if(length(npt) < 2) npt <- rep(npt, 2)
  n1 <- npt[1]
  n2 <- npt[2]
  x <- x1 <- if(length(xlim) > 2) xlim else { if(length(xlim) ==2)  
             seq(min(xlim), max(xlim), length=n1) else stop("invalid xlim") }
  y <- x2 <- if(length(ylim) > 2) ylim else { if(length(ylim) ==2)  
             seq(min(ylim), max(ylim), length=n2) else stop("invalid ylim") }
  x1.x2 <- cbind(rep(x1, n2), as.vector(matrix(x2, n1, n2, byrow=TRUE)))
  X <- matrix(x1.x2, n1 * n2, 2, byrow = FALSE)
  f.val <- matrix(funct(X, ...), n1, n2)
  if(missing(grf)) grf <- "contour"  
  if(!(grf %in% c("contour", "filled.contour", "persp", "image")))
     stop(gettextf("'%s' is not a valid 'grf' function", grf))
  cmd.char <- paste(paste("graphics::", grf, "(x, y, f.val", sep=""),
                    if(missing(grpar)) ")" else paste(",", grpar, ")"))
  eval(parse(text=cmd.char))                  
  invisible(list(x=x, y=y, pts=X, f.values=f.val))
}
