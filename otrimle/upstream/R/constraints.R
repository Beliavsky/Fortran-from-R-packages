.GssERC <- function(values, erc, Lmin, Lmax, sumtau, P, G) {
  tol <- sqrt(.Machine$double.eps)
  iter.max <- 99
  iter <- 0L
  grc <- {
    sqrt(5) - 1
  }/2
  a <- log(max(Lmin, .Machine$double.eps))
  b <- log(min(Lmax, .Machine$double.xmax))
  tL <- values
  l <- a
  x1 <- a + {
    1 - grc
  } * {
    b - a
  }
  x2 <- a + grc * {
    b - a
  }
  iter <- iter + 1L
  l <- exp(x1)
  tL <- values
  tL[values < l] <- l
  tL[values > erc * l] <- l * erc
  f1 <- sum(sumtau * .colSums({
    log(tL) + values/tL
  }, m = P, n = G, na.rm = FALSE))
  iter <- iter + 1L
  l <- exp(x2)
  tL <- values
  tL[values < l] <- l
  tL[values > erc * l] <- l * erc
  f2 <- sum(sumtau * .colSums({
    log(tL) + values/tL
  }, m = P, n = G, na.rm = FALSE))
  while ({
    abs(b - a) > tol
  } & {
    iter < iter.max
  }) {
    if (f1 > f2) {
      a <- x1
      x1 <- x2
      x2 <- a + grc * {
        b - a
      }
      f1 <- f2
      iter <- iter + 1L
      l <- exp(x2)
      tL <- values
      tL[values < l] <- l
      tL[values > erc * l] <- l * erc
      f2 <- sum(sumtau * .colSums({
        log(tL) + values/tL
      }, m = P, n = G, na.rm = FALSE))
    }
    else {
      b <- x2
      x2 <- x1
      x1 <- a + {
        1 - grc
      } * {
        b - a
      }
      f2 <- f1
      l <- exp(x1)
      iter <- iter + 1L
      tL <- values
      tL[values < l] <- l
      tL[values > erc * l] <- l * erc
      f1 <- sum(sumtau * .colSums({
        log(tL) + values/tL
      }, m = P, n = G, na.rm = FALSE))
    }
  }
  l <- exp({
    a + b
  }/2)
  tL <- values
  tL[values < l] <- l
  tL[values > erc * l] <- l * erc
  return(tL)
}
.EqNPC <- function(x, icd, sumtau.old, pr, psi, N, npr.max) {
  if (icd * x == 0) {
    return(-N * npr.max)
  }
  else {
    return(sum(icd * x/{
      icd * x + {
        {
          1 - x
        }/{
          N - sumtau.old[1]
        }
      } * {
        {
          psi - pr[1] * icd
        } * N
      }
    }) - N * npr.max)
  }
}
.EqERC <- function(x, data, icd, M, L, V, pr, erc, npr.max, N, G, P) {
  gausscost <- {
    2 * pi
  }^{
    -0.5 * P
  }
  L[L <= x] <- x
  L[L >= x * erc] <- x * erc
  qtau <- matrix(0, nrow = N, ncol = {
    G + 1
  }, byrow = TRUE, dimnames = NULL)
  qtau[, 1] <- pr[1] * icd
  Xpp <- matrix(0, nrow = P, ncol = P, byrow = TRUE, dimnames = NULL)
  Delta <- matrix(0, nrow = N, ncol = P, byrow = TRUE, dimnames = NULL)
  smd <- rep(0, times = N)
  for (j in 1:G) {
    Xpp <- V[, , j] %*% {
      {
        1/L[, j]
      } * t(V[, , j])
    }
    Delta <- {
      data - matrix(M[, j], nrow = N, ncol = P, byrow = TRUE, dimnames = NULL)
    }
    smd <- .rowSums({
      Delta %*% Xpp
    } * Delta, m = N, n = P, na.rm = TRUE)
    qtau[, 1 + j] <- pr[1 + j] * gausscost * {
      {
        prod(L[, j])
      }^{
        -0.5
      }
    } * {
      exp(-0.5 * smd)
    }
  }
  ans <- sum(qtau[, 1]/.rowSums(qtau, m = N, n = {
    G + 1
  }, na.rm = TRUE)) - N * npr.max
  return(ans)
}
