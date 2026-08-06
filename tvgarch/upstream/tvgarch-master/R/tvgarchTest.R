tvgarchTest <- function (y, xtv = NULL, alpha = 0.05)
{
  y <- as.matrix(y)
  garch11Est <- tvgarch(y = y, order.g = NULL)
  theta_hat <- coef.tvgarch(garch11Est)
  h <- as.matrix(fitted.tvgarch(garch11Est))
  n <- length(y)
  y2 <- y^2
  if (!is.null(xtv)) {
    if (length(xtv) != n) {
      stop("Length of xtv unequal to the number of observations.")
    }
    xtv <- matrix(xtv, n, 1)
    if (is.null(colnames(xtv))) colnames(xtv) <- "xtv"
  }
  if (is.null(xtv))  {
    xtv <- matrix((1:n)/n, n, 1)
    colnames(xtv) <- "time"
  }
  z <- as.numeric(y2/h-1)
  v <- cbind(1, y2, h)
  v <- rbind(c(1, mean(y2), mean(y2)), v)
  v <- v[1:n,]
  one <- rep(1, n) 
  RSS0 <- t(z) %*% z
  inv_h <- as.numeric(1/h)
  dhdw <- matrix(0, n, 3)
  dhdw[1,] <- v[1,]
  for(i in 2:n){
    dhdw[i,] <- v[i,] + theta_hat[3]*dhdw[(i-1),]
  }
  z0 <- one
  z1 <- xtv
  z2 <- cbind(xtv, xtv^2)
  z3 <- cbind(xtv, xtv^2, xtv^3)
  x0 <- cbind(inv_h * dhdw, z0)
  x1 <- cbind(inv_h * dhdw, z0, z1)
  x2 <- cbind(inv_h * dhdw, z0, z2)
  x3 <- cbind(inv_h * dhdw, z0, z3)
  ' 
    Non-Robust LM Test
  '
  RSS10 <- sum(lm(z ~ x0 -1)$residuals^2)
  RSS11 <- sum(lm(z ~ x1 -1)$residuals^2)
  RSS12 <- sum(lm(z ~ x2 -1)$residuals^2)
  RSS13 <- sum(lm(z ~ x3 -1)$residuals^2)
  LM1 <- n * (RSS0 - RSS11) / RSS0 
  LM3 <- n * (RSS0 - RSS13) / RSS0 
  df1 <- NCOL(z1)
  df3 <- NCOL(z3)
  pval1 <- pchisq(LM1, df1, lower.tail = FALSE)
  pval3 <- pchisq(LM3, df3, lower.tail = FALSE)
  LMB3 <- n * (RSS12 - RSS13) / RSS12
  dfB3 <- 1
  pvalB3 <- pchisq(LMB3, dfB3, lower.tail = FALSE)
  LMB2 <- n * (RSS11 - RSS12) / RSS11
  dfB2 <- 1
  pvalB2 <- pchisq(LMB2, dfB2, lower.tail = FALSE)
  LMB1 <- n * (RSS0 - RSS11) / RSS0
  dfB1 <- 1
  pvalB1 <- pchisq(LMB1, dfB1, lower.tail = FALSE)
  mat = round(rbind(c(LM3,pval3), c(LMB3,pvalB3), c(LMB2,pvalB2), 
                    c(LMB1,pvalB1)), 4)
  colnames(mat) = c("NonRobTR2", "p-value")
  rownames(mat) = c("H0:B3=B2=B1=0", "H03:B3=0", "H02:B2=0|B3=0", 
                    "H01:B1=0|B3=B2=0")
  ' 
    Robust LM Test
  '
  Y33 <- cbind(z0, z3) 
  X33 <- inv_h * dhdw
  C33 <- lm(Y33 ~ X33 - 1)
  ResRob <- z * C33$residuals
  one <- rep(1, n)
  RSS0Rob <- sum(lm(one ~ ResRob - 1)$residuals^2)
  LMRob <- n - RSS0Rob
  pvalRob <- pchisq(LMRob, df3, lower.tail = FALSE)
  Y3 <- xtv^3
  X33 <- cbind(inv_h * dhdw, z0, z2)
  C3 <- lm(Y3 ~ X33 -1)
  ResRob3 <- z * C3$residuals
  RSS3Rob <- sum(lm(one ~ ResRob3 -1)$residuals^2)
  LMRob3 <- n - RSS3Rob
  pvalRob3 <- pchisq(LMRob3, dfB3, lower.tail = FALSE)
  Y2 <- xtv^2
  X22 <- cbind(inv_h * dhdw, z0, z1)
  C2 <- lm(Y2 ~ X22 -1)
  ResRob2 <- z * C2$residuals
  RSS2Rob <- sum(lm(one ~ ResRob2 -1)$residuals^2)
  LMRob2 <- n - RSS2Rob
  pvalRob2 <- pchisq(LMRob2, dfB2, lower.tail = FALSE)
  Y1 <- cbind(z0, z1)
  X11 <- inv_h * dhdw
  C1 <- lm(Y1 ~ X11 -1)
  ResRob1 <- z * C1$residuals
  RSS1Rob <- sum(lm(one ~ ResRob1 -1)$residuals^2)
  LMRob1 <- n - RSS1Rob
  pvalRob1 <- pchisq(LMRob1, dfB1, lower.tail = FALSE)
  matRob = round(rbind(c(LMRob,pvalRob), c(LMRob3,pvalRob3), c(LMRob2,pvalRob2), 
                       c(LMRob1,pvalRob1)),4)
  colnames(matRob) = c("RobTR2", "p-value")
  rownames(matRob) = c("H0:B3=B2=B1=0", "H03:B3=0", "H02:B2=0|B3=0", 
                       "H01:B1=0|B3=B2=0")
  order.g <- as.matrix(NA)
  colnames(order.g) <- "Single"
  rownames(order.g) <- paste("No. of locations (alpha = ", alpha, ")", sep = "")
  if (pvalRob < alpha) {
    if (pvalRob2 < pvalRob1 && pvalRob2 < pvalRob3 && pvalRob2 < alpha) {
      order.g[1,1] <- 2
    }
    if (pvalRob1 < pvalRob3 && pvalRob1 < alpha) order.g[1,1] <- 1
    if (pvalRob3 < pvalRob1 && pvalRob3 < alpha) order.g[1,1] <- 3
  }
  if (pvalRob >= alpha) {
    order.g[1,1] <- 0
    warning("The unconditional variance is constant.")
  }
  results <- list(mat = mat, matRob = matRob, xtv = xtv, order.g = order.g,
                  garch11 = garch11Est, date = date())
  class(results) <- "tvgarchTest"
  return(results)
}