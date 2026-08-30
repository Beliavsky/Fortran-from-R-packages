which.tile <- function(x,y,tl){
  u  <- c(x,y)
  nt <- length(tl)
  d2 <- numeric(nt)
  for(i in 1:nt) {
    d2[i] <- sum((u-tl[[i]]$pt)^2)
  }
  which.min(d2)
}
