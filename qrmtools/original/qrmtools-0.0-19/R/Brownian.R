### Dependent Brownian and related motions #####################################

##' @title Simulating Paths from Brownian Motions, Geometric Brownian Motions
##'        and Brownian Bridges
##' @param N number of paths
##' @param t (n+1)-vector (t_0,...,t_n) with 0 = t_0 < ... < t_n containing
##'        the time points at which to simulate
##' @param d dimension of the stochastic process
##' @param U (N * n, d)-matrix of copula realizations used as increments of the
##'        stochastic process.
##' @param drift numeric(1) or d-vector of drifts mu (risk-neutral drifts:
##'        r - volas^2/2 for the risk-free interest rate r)
##' @param vola numeric(1) or d-vector of volatilities sigma
##' @param type character vector indicating the type of stochastic process
##'        sampled (Brownian motion, geometric Brownian motion or Brownian
##'        bridge)
##' @param init numeric(1) or d-vector of initial values used if type = "GBM"
##' @return (N, n+1, d)-array
##' @author Marius Hofert
rBrownian <- function(N, t, d = 1, U = matrix(runif(N * n * d), ncol = d),
                      drift = 0, vola = 1, type = c("BM", "GBM", "BB"), init = 1)
{
    ## Checks and setup
    if(length(drift) == 1) drift <- rep(drift, d)
    if(length(vola) == 1) vola <- rep(vola, d)
    if(length(init) == 1) init <- rep(init, d)
    stopifnot(N >= 1, (n <- length(t)-1) >= 0, t[1] == 0, t >= 0, (delta.t <- diff(t)) > 0,
              d >= 1, (dim.U <- dim(U)) == c(N * n, d), length(drift) == d,
              length(vola) == d, vola > 0, length(init) == d, init > 0)
    ## Main
    type <- match.arg(type)
    switch(type,
           "BM" = {
               ## X_{t_k,j} = \mu_j t_k + \sigma_j * W_{t_k,j}
               ##           = \mu_j t_k + \sigma_j * \sum_{i=1}^{k} \sqrt{t_{i}-t_{i-1}} * \Phi^{-1}(U_{i,j})
               ## Preliminaries
               dmnms <- list("path" = NULL, "time" = NULL, "component" = NULL) # dimension names for returned arrays
               sqrt.delta.t <- matrix(sqrt(delta.t), nrow = n, ncol = d) # (n, d)-matrix of repeated n-vectors delta.t
               Y <- array(, dim = c(N, n, d), dimnames = dmnms) # for W_{t_k,j} summands (added up to obtain W_{t_k,j})
               for(i in 1:N) # iterate over paths and convert U to (N, n, d)-array of increments
                   Y[i,,] <- sqrt.delta.t * qnorm(U[n * (i-1) + 1:n, ]) # for fixed i: (n, d)-matrix
               ## Build the process
               mu <- rep(drift, each = N) # treated as an (N, d)-matrix
               sigma <- rep(vola, each = N) # treated as an (N, d)-matrix
               W.t <- matrix(0, nrow = N, ncol = d) # (N, d)-matrix to store W_{t_k} for a fixed k in the next for loop
               X <- array(, dim = c(N, n+1, d), dimnames = dmnms)
               X[,1,] <- 0 # deal with t_0 = t[1]; this is why we need t[1] == 0 and not just t[1] >= 0 (value of X!)
               for(k in 1:n) { # iterate over time points
                   W.t <- W.t + Y[,k,] # update Y cumulative sums (to get W_{t_k,j} for all j = 1,...,d and all paths)
                   X[,k+1,] <- mu * t[k+1] + sigma * W.t # for fixed k: (N, d)-matrix
               }
               X
           },
           "GBM" = {
               ## S_{t_k,j} = S_{0,j} * e^{X_{t_k,j}}
               X <- rBrownian(N, t = t, d = d, U = U, drift = drift, vola = vola, type = "BM") # (N, n+1, d)-array
               S.0 <- array(rep(init, each = N * (n+1)), dim = c(N, n+1, d)) # (N, n+1, d)-array
               S.0 * exp(X) # note: X[,1,] is 0, so S[,1,] is S_0 (as it should)
           },
           "BB" = {
               ## B_{t_k,j} = X_{t_k,j} - (t/T) X_{t_n,j}
               dmnms <- list("path" = NULL, "time" = NULL, "component" = NULL) # dimension names for returned arrays
               X <- rBrownian(N, t = t, d = d, U = U, drift = drift, vola = vola, type = "BM") # (N, n+1, d)-array
               t.T <- t / t[n+1] # (n+1)-vector (0, t_1/T, ..., t_n/T = 1)
               X.T <- X[,n+1,] # (N, d)-matrix
               B <- array(, dim = c(N, n+1, d), dimnames = dmnms)
               for(k in 1:(n+1)) { # by far the easiest version; both sapply() and array() require reorderings
                   B[,k,] <- X[,k,] - t.T[k] * X.T
               }
               B
           })
}


##' @title Extract the Increments from Given Paths of (Geometric) Brownian Motions
##' @param x (n+1)-vector containing an (n+1)-path of a (geometric) Brownian motion
##'        or (n+1, d)-matrix containing an (n+1)-path of d (geometric) Brownian motions
##'        or (N, n+1, d)-array containing N (n+1)-paths of d (geometric) Brownian motions
##' @param t (n+1)-vector (t_0,...,t_n) with 0 = t_0 < ... < t_n containing the time
##'        points at which to filter
##' @param drift d-vector of drifts mu (risk-neutral drifts: r - volas^2/2 for
##'        the risk-free interest rate r)
##' @param vola d-vector of volatilities sigma
##' @param type character vector indicating the type of stochastic process
##'        filtered (Brownian motion or geometric Brownian motion)
##' @return (N, n, d)-array
##' @author Marius Hofert
##' @note x[.,,.] and t have to be of the same length n + 1 (e.g., as returned by rBrownian())
##'       => could have used 'n' instead of 'n+1', but then not consistent with rBrownian()
deBrowning <- function(x, t, drift = 0, vola = 1, type = c("BM", "GBM"))
{
    ## Checks
    stopifnot((n <- length(t)-1) >= 0, t[1] == 0, t >= 0, (delta.t <- diff(t)) > 0, vola > 0)
    if(n < 1) stop("Need at least two time points.") # to compute one difference; note: n = length(t) - 1 => length(t) >= 2
    dlen <- length(dim(x))
    if(dlen == 0) {
        x <- cbind(x) # convert x to an (n+1, d)-matrix for d = 1
        dlen <- 2
    }
    if(dlen == 2) {
        x <- array(x, dim = c(1, n+1, ncol(x))) # convert x to a (N, n+1, d)-array for N = 1
        dlen <- 3
    }
    if(dlen != 3)
        stop("'x' must be a length(t)-vector, (length(t), d)-matrix or (N, length(t), d)-array.")

    ## Preliminaries
    dm <- dim(x)
    N <- dm[1]
    d <- dm[3]
    if(length(drift) == 1) drift <- rep(drift, d)
    if(length(vola) == 1) vola <- rep(vola, d)
    stopifnot(length(drift) == d, length(vola) == d)
    if(dm[2] != n + 1) stop("Need dim(x)[2] == length(t).")

    ## Deal with GBM case
    type <- match.arg(type)
    if(type == "GBM") {
        ## x is S here and S_{t,j} = S_{0,j} * \exp(X_{t,j})
        ## Recover the X's from the given S (= x) via X_{t,j} = log(S_{t,j} / S_{0,j})
        X <- log(sweep(x, MARGIN = c(1,3), STATS = x[,1,], FUN = "/")) # (N, n+1, d); X_{t_0,j} = 0,..., X_{t_n,j}
        ## Recover the Z's from the X's
        return(deBrowning(X, t = t, drift = drift, vola = vola, type = "BM"))
    }

    ## Now we are in the BM case
    ## Get rid of drift and vola via W_{t_k,j} = \sum_{i=1}^k \sqrt{t_{i}-t_{i-1}} Z_{i,j} = (X_{t_k,j} - \mu_j t_k) / \sigma_j
    drift.t. <- outer(t, drift) # (n+1, d)-matrix
    drift.t <- array(rep(drift.t., each = N),         dim = c(N, n+1, d)) # (N, n+1, d)-array
    vola    <- array(rep(vola,     each = N * (n+1)), dim = c(N, n+1, d)) # (N, n+1, d)-array
    W <- (x - drift.t) / vola # (N, n+1, d)-array

    ## Recover the Z's from the W's via W_{t_k,j} = \sum_{i=1}^k \sqrt{t_{i}-t_{i-1}} Z_{i,j}
    delta.W <- aperm(apply(W, c(1, 3), diff), perm = c(2, 1, 3)) # (N, n, d)-array
    stopifnot(dim(delta.W) == c(N, n, d))
    sweep(delta.W, MARGIN = 2, STATS = sqrt(delta.t), FUN = "/") # use that Z_{k+1} = (W_{t_{k+1}} - W_{t_k}) / \sqrt{t_{k+1}-t_{k}}
}
