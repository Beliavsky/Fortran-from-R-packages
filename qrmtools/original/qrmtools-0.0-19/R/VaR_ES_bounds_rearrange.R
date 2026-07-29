### Tools for computing the VaR_alpha and ES_alpha bounds ######################

### 1 Worst/best VaR, best ES in the inhomogeneous case ########################

### 1.1 Auxiliary functions and workhorses #####################################

##' @title Determine the indices which order any increasing (!) vector y
##'        oppositely to x
##' @param x vector
##' @return order(order(x, decreasing = TRUE)) (= N+1-rank(x))
##' @author Marius Hofert
##' @note For convergence of rearrange() it is crucial to have a stable sorting
##'       procedure underlying (as then no swaps on ties back and forth until
##'       eternity take place which decreases the probability of non-convergence).
##'       The various methods like qsort() in C or rsort_with_index() are *not*
##'       stable. In the way we need it here, rank(, ties.method = "last") would
##'       be as well, but internally uses order() and thus is not faster.
##'       However, we can make order() faster by using R_orderVector1() instead
##'       of R_orderVector() (available from 3.2.3 onwards)
##'       => For d = 1000 and N = 16384, this only brought an improvement of 1.3%,
##'          though.
##'       If that turns out to be faster (in the future), use:
##'       indices_opp_ordered_to <- if(getRversion() >= "3.2.3")
##'       {
##'          function(x) .Call(C_indices_opp_ordered_to, x)
##'       } else {
##'          function(x) order(order(x, decreasing = TRUE))
##'       }
indices_opp_ordered_to <- function(x) order(order(x, decreasing = TRUE))

##' @title Compute the number of columns oppositely ordered to the sum of all others
##' @param x (N, d)-matrix
##' @return number of columns oppositely ordered to the sum of all others
##' @author Marius Hofert
##' @note Same time-saving tricks as behind rearrange(), RA() and ARA() (work
##'       with list of columns of x)
num_of_opp_ordered_cols <- function(x) {
    x.rs <- .rowSums(x, nrow(x), ncol(x)) # faster than rowSums()
    x.lst <- .Call(col_split, x) # to avoid indexing the jth column, we work with a list!
    x.lst.sorted <- lapply(x.lst, sort.int) # sorting is only necessary once!
    sum(vapply(seq_len(ncol(x)),
               function(j) {
                   xj <- x.lst[[j]]
                   all(x.lst.sorted[[j]][indices_opp_ordered_to(x.rs - xj)]
                       == xj)
               }, NA))
}

##' @title Basic rearrangement function for (A)RA
##' @param X (N, d)-matrix \underline{X}^\alpha or \overline{X}^\alpha
##' @param tol tolerance to determine (the individual) convergence;
##'        if NULL, column rearrangements are done until the matrix doesn't
##'        change anymore n.lookback consecutive times
##' @param tol.type character string indicating the tolerance function used
##'        ("relative" or "absolute")
##' @param n.lookback number of column rearrangements to look back for deciding
##'        about 'convergence' (must be a number in {1,.., max.ra-1}, typically
##'        around d).
##' @param max.ra maximal number of column rearrangements considered
##' @param method character indicating whether worst VaR, best VaR or best ES
##'        is approximated (determines optimizing function)
##' @param sample logical indicating whether each column of the working
##'        matrix is randomly permuted before the rearrangements begin
##' @param is.sorted logical indicating whether X is columnwise sorted in
##'        increasing order
##' @param trace logical indicating whether the underlying matrix is
##'        printed after each rearrangement step
##' @param ... additional arguments passed to the internal optimization
##'        function optim.fun()
##' @return list containing the
##'         1) Computed (lower or upper [depending on X]) bound for (worst or
##'            best [depending on method]) VaR or best ES
##'         2) (Individual) tolerance reached (i.e., the absolute/relative change
##'            of the row sum when looking back n.lookback column rearrangements)
##'         3) Logical indicating whether the algorithm has converged
##'         4) Vectors of optimal (minimal for worst VaR; maximal for best VaR;
##'            ES_alpha for best ES) row sums after each considered column
##'            rearrangement
##'         5) The (optimally) rearranged (N, d)-matrix
##'         6) The (averaged) rows of the (optimally) rearranged matrix leading
##'            to the final optimal row sum
##' @author Marius Hofert and Kurt Hornik
##' @note - We use "<= tol" to determine convergence instead of "< tol" as
##'         this then also nicely works with "= 0" (if tol = 0) which stops in
##'         case the matrices are identical (no change at all).
##'       - We conduct checks of convergence after rearranging each column after the
##'         dth (not only after rearranging all d columns)
##'       - The columns of X have to be given in increasing order if is.sorted = TRUE
##'       - No checking here due to speed! Note that max.ra must be > ncol(X)
##'       - Idea: The column rearrangements should provide the maximal minimal
##'               row sum (for worst VaR), the minimal maximal row sum
##'               (for best VaR) or the minimal ES (for best ES).
rearrange <- function(X, tol = 0, tol.type = c("relative", "absolute"),
                      n.lookback = ncol(X), max.ra = Inf,
                      method = c("worst.VaR", "best.VaR", "best.ES"),
                      sample = TRUE, is.sorted = FALSE, trace = FALSE, ...)
{
    ## Setup
    N <- nrow(X)
    d <- ncol(X)
    tol.type <- match.arg(tol.type)
    method <- match.arg(method)

    ## Define helper functions
    optim.fun <- switch(method,
    "worst.VaR" = {
        min
    },
    "best.VaR" = {
        max
    },
    "best.ES" = {
        level <- NULL # make CRAN check happy
        stopifnot(hasArg(level)) # check if 'level' has been provided (via '...')
        function(x) ES_np(x, level = list(...)$level)
    },
    stop("Wrong 'method'"))
    tol.fun <- if(tol.type == "absolute") {
        function(x, y) abs(x-y)
    } else {
        function(x, y) abs((x-y)/y)
    }

    ## Tracing
    if(trace) {
        B <- X
        colnames(B) <- rep("", d)
        print(B)
    }

    ## Determine X as a list (X.lst), its row sums (X.rs) and sorted version (X.lst.sorted)
    X.lst <- .Call(col_split, X)
    X.lst.sorted <- if(is.sorted) {
        X.lst # split the given X (since it's already sorted) into columns
    } else {
        lapply(X.lst, sort) # list of (resampled) columns of X
    }
    if(sample) X.lst <- lapply(X.lst, sample) # list of (resampled) columns of X
    X.rs <- .rowSums(do.call(cbind, X.lst), N, d) # row sums of X

    ## Setup for iteration
    num.cols.no.change <- 0 # number of consecutively rearranged columns with no change
    len.opt.row.sums <- 64 # length of vector of optimal (minimal/maximal/ES) row sums after each rearranged column
    opt.row.sums <- numeric(len.opt.row.sums) # vector (will be doubled in size if necessary; faster than c()ing to it all the time)
    is.null.tol <- is.null(tol)

    ## Rearrange columns one at a time
    iter <- 0 # current iteration number
    j <- 0 # current column number
    while (TRUE) {

        ## Update the running indices
        iter <- iter + 1 # current iteration number (in IN)
        j <- if(j >= d) 1 else j+1 # current column

        ## Update the working 'matrix'
        Y.lst <- X.lst # define 'matrix' Y (former 'matrix' X) to work with
        Y.rs <- X.rs # row sum of Y (= row sum of X)

        ## Oppositely order the jth column to the sum of all others
        yj <- Y.lst[[j]] # pick out jth column (could contain -/+Inf)
        rs.mj <- Y.rs - yj # sum over all other columns (but the jth)
        ## Note: Consider a specific row and suppose Y.rs = -/+Inf.
        ##       If yj = -/+Inf, rs.mj contains NaN (as Inf - Inf = NaN).
        ##       Now it depends whether yj is the only -/+Inf in that row
        ##       (=> rs.mj should be finite, but we can't recover it from Y.rs)
        ##       or not (=> rs.mj should be -/+Inf).
        ##       => rearrange() should only be called with finite values;
        ##          the speed-up idea (and working on columns only) does not
        ##          easily seem to allow for an extension to -/+Inf values.
        yj. <- X.lst.sorted[[j]][indices_opp_ordered_to(rs.mj)] # oppositely reorder Y_j
        ## Note: The elements of X.lst.sorted are sorted in increasing order
        ##       which is required for oppositely reordering them
        ## Update the working 'matrix' and vector of row sums
        Y.lst[[j]] <- yj. # update with rearranged jth column
        Y.rs <- rs.mj + yj. # update row sum of Y

        ## Tracing
        if(trace) {
            B <- do.call(cbind, Y.lst)
            colnames(B) <- rep("", d)
            no.change <- identical(yj, yj.)
            colnames(B)[j] <- if(no.change) "=" else "|"
            B <- cbind(B, rs.mj, sum = .rowSums(B, m = N, n = d))
            colnames(B)[d+1] <- "-col" # paste0("-",j)
            print(B)
        }

        ## Update the vector of computed optimal row sums
        opt.rs.cur.col <- optim.fun(Y.rs) # compute new optimal row sum
        if(iter > len.opt.row.sums) { # if iter in {128 + 1, 256 + 1, ...} => double the size
            opt.row.sums <- c(opt.row.sums, numeric(len.opt.row.sums))
            len.opt.row.sums <- 2 * len.opt.row.sums
        }
        opt.row.sums[iter] <- opt.rs.cur.col # append it

        ## Check convergence
        ## Idea: After a column has been rearranged, compute the tol (and thus
        ##       determine convergence) between the optimal row sum
        ##       after that rearrangement and from n.lookback steps before (for the default:
        ##       when that column was rearranged the last time). The earliest we check for
        ##       convergence is when iter > d.
        ## Note: - This is a bit more elegant than the original RA which checked only
        ##         on j = d, not after rearranging *each* column.
        ##       - Checking only *two* consecutive columns led to a bad behavior for ARA()
        ##         in some cases (e.g., real OpRisk data): Both the individual and the joint
        ##         relative tolerances were satisfied but far off (with reltol[1] = 0.001).
        ##         Of course one could check n.lookback consecutive columns for *all of them
        ##         together* to fulfill the 'convergence' criterion, but then what's the reached
        ##         tolerance tol if more than two columns are involved? Maybe the maximum
        ##         tolerance computed over all previous n.lookback-many rearranged columns?
        ##         There's probably no gain in doing that.
        if(is.null.tol) { # tol = NULL
            num.cols.no.change <- if(identical(yj, yj.)) num.cols.no.change + 1 else 0
            if(num.cols.no.change == n.lookback) { # => matrix has not changed in n.lookback consecutive col rearrangements
                tol. <- 0 # as there was no change
                tol.reached <- TRUE # as we reached 'no change' in n.lookback consecutive steps (we don't care whether max.ra has been reached)
                break
            } else { # check whether we have to stop due to max.ra
                if(iter == max.ra) { # need max.ra > n.lookback (as otherwise (*) is wrong)
                    ## Note: iter = number of columns we have already rearranged
                    opt.rs.n.lookback.col.ago <- opt.row.sums[iter-n.lookback] # (*)
                    tol. <- tol.fun(opt.rs.cur.col, opt.rs.n.lookback.col.ago) # compute the attained tolerance (in comparison to n.lookback column rearrangements ago)
                    tol.reached <- FALSE # as num.cols.no.change < n.lookback
                    break
                }
            }
        } else { # tol numeric >= 0
            if(iter > n.lookback) { # ... now we can look back n.lookback rearrangements
                opt.rs.n.lookback.col.ago <- opt.row.sums[iter-n.lookback]
                tol. <- tol.fun(opt.rs.cur.col, opt.rs.n.lookback.col.ago) # compute the attained tolerance (in comparison to n.lookback rearrangements ago)
                tol.reached <- tol. <= tol
                if(iter == max.ra || tol.reached) break # also here we need max.ra > n.lookback; see (*)
            }
        }

        ## Updates for the next column rearrangement
        X.lst <- Y.lst # update the working 'matrix'
        X.rs <- Y.rs # update the row sums

    } # while()

    ## Determine (column means of) the row(s) of the final rearranged matrix
    ## which correspond to the optimal row sum
    X.rearranged <- do.call(cbind, Y.lst) # rearranged X
    opt.rs.ind <- which(Y.rs == opt.rs.cur.col) # indices of rows which lead to optimal row sums
    X.rearranged.opt.row <- colMeans(X.rearranged[opt.rs.ind,, drop = FALSE]) # average over multiple rows if optimum is not unique

    ## Return
    list(bound = opt.rs.cur.col, # computed bound (\underline{s}_N or \overline{s}_N)
         tol = tol., # tolerance for the computed bound
         converged = tol.reached, # indicating whether converged
         opt.row.sums = opt.row.sums[1:iter], # the computed optimal row sums after each rearrangement
         X.rearranged = X.rearranged, # the (single) rearranged matrix X
         X.rearranged.opt.row = X.rearranged.opt.row) # (averaged) rows of X.rearranged leading to the final optimal row sum
}

##' @title Basic rearrangement function for (A)BRA
##' @param X see rearrange()
##' @param tol see rearrange() but can't be NULL (has to be >= 0) as we would
##'        need to check after each block rearrangement whether all columns are
##'        oppositely ordered to the sum of all others (or whether all possible
##'        blocks are oppositely ordered to all others) -- too time-consuming.
##' @param tol.type character string indicating the tolerance function used
##'        ("absolute" or "relative"); here the default is "absolute" as we
##'        minimize the variance of the row sums and this variance can be 0
##'        (we would thus divide by 0 for the relative variance)
##' @param n.lookback number of block rearrangements to look back for deciding
##'        about 'convergence' (via considering the absolute/relative change
##'        in the row sum variance; must be in {1,.., max.ra-1})
##' @param max.ra maximal number of block rearrangements considered
##' @param method character indicating whether worst VaR, best VaR or best ES
##'        is approximated (determines optimizing function)
##' @param sample see rearrange()
##' @param ... additional arguments passed to the internal optimization
##'        function optim.fun()
##' @return list containing the
##'         1) Computed (lower or upper [depending on X]) bound for (worst or
##'            best [depending on method]) VaR or best ES
##'         2) (Individual) tolerance reached (i.e., the absolute/relative change
##'            of the row sum variance when looking back n.lookback block
##'            rearrangements
##'         3) Logical indicating whether the algorithm has converged
##'         4) Vectors of optimal (minimal for worst VaR; maximal for best VaR;
##'            ES_alpha for best ES) row variances after each considered block
##'            rearrangement
##'         5) The (optimally) rearranged (N, d)-matrix
##'         6) The (averaged) rows of the (optimally) rearranged matrix leading
##'            to the final optimal row sum
##' @author Marius Hofert and Martin Stefanik
##' @note - This is a modified version of Bernard and McLeish (2014)
##'       - Idea: The variance of the row sums is minimized through block rearrangements,
##'               which should provide the maximal minimal row sum (for worst VaR), the
##'               minimal maximal row sum (for best VaR) or the minimal ES (for best ES).
block_rearrange <- function(X, tol = 0, tol.type = c("absolute", "relative"),
                            n.lookback = ncol(X), max.ra = Inf,
                            method = c("worst.VaR", "best.VaR", "best.ES"),
                            sample = TRUE, trace = FALSE, ...)
{
    ## Setup
    N <- nrow(X)
    d <- ncol(X)
    tol.type <- match.arg(tol.type)
    method <- match.arg(method)

    ## Define helper functions
    optim.fun <- switch(method,
    "worst.VaR" = {
        min
    },
    "best.VaR" = {
        max
    },
    "best.ES" = {
        level <- NULL # make CRAN check happy
        stopifnot(hasArg(level)) # check if 'level' has been provided (via '...')
        function(x) ES_np(x, level = list(...)$level)
    },
    stop("Wrong 'method'"))
    tol.fun <- if(tol.type == "absolute") {
        function(x, y) abs(x-y)
    } else {
        function(x, y) abs((x-y)/y)
    }

    ## Tracing
    if(trace) {
        B <- X
        colnames(B) <- rep("", d)
        print(B)
    }

    ## Sample, row sums, row sum variance
    if(sample) X <- apply(X, 2, sample)
    X.rs <- .rowSums(X, m = N, n = d)

    ## Setup for iteration
    len.opt.row.sums <- 64
    opt.row.sums <- numeric(len.opt.row.sums) # vector (will be doubled in size if necessary; faster than c()ing to it all the time)
    rs.var <- numeric(n.lookback) # variances over the last n.lookback column rearrangements

    ## Rearrange blocks one at a time
    iter <- 0 # current iteration number
    opt.row.vars <- numeric(n.lookback)
    while (TRUE) {

        ## Update the running indices
        iter <- iter + 1

        ## Sample a random two-set partition
        bsize <- sample(seq_len(d-1), size = 1) # draw a number from {1,..,d-1} (block size)
        if(bsize > d/2) bsize <- d - bsize # speed-up; rearrange smaller block
        block <- sample(seq_len(d), size = bsize) # draw bsize numbers from {1,..,d} (actual block of size bsize)
        rs.block <- .rowSums(X[, block], m = N, n = bsize) # row sum of the block
        ## Note: We can't avoid the call of .rowSums() here if blocks of size > 1 are
        ##       rearranged => performance drawback in comparison to rearrange().
        rs.complement.block <- X.rs - rs.block # row sum of the complement block (= set of remaining columns)

        ## Oppositely order all columns belonging to block 'block' to the row
        ## sums of the complement block
        if(trace) B <- X # grab out matrix before rearranging
        ord <- order(rs.block)
        X.block.ordered <- X[ord, block, drop = FALSE] # (*)
        oord <- indices_opp_ordered_to(rs.complement.block)
        X[, block] <- X.block.ordered[oord,]
        ## Concerning (*): - In rearrange() we could work with X.lst.sorted and thus use
        ##                   presorted columns. Avoiding X[order(rs.block), block]
        ##                   doesn't seem possible here as rs.block can change all the time.
        ##                   This is one reason why block_rearrange() is slower than rearrange().
        ##                 - Another reason is that we work with matrices here, not lists
        ##                   (the effect of this is probably rather minor given the sorting)
        ## Update the vector of row sums
        X.rs <- rs.complement.block + rs.block[ord][oord] # formerly: rs.complement.block + .rowSums(X[, block], m = N, n = bsize)
        ## Note: rs.block[ord][oord] = .rowSums(X[, block], m = N, n = bsize) (but avoiding to
        ##       call *expensive* .rowSums() -- even if .rowSums() is already C code, see ./src/main/array.c -> do_colsum)

        ## Tracing
        if(trace) {
            colnames(B) <- rep("", d)
            colnames(B)[block] <- vapply(block, function(j)
                if(identical(X[,j], B[,j])) "=" else "|", character(0))
            B <- cbind(B, rs.complement.block, sum = .rowSums(B, m = N, n = d))
            colnames(B)[d+1] <- "-block"
            print(B)
        }

        ## Update the vector of computed optimal row sums
        opt.rs.cur.col <- optim.fun(X.rs) # compute new optimal row sum
        if(iter > len.opt.row.sums) { # if iter in {128 + 1, 256 + 1, ...} => double the size
            opt.row.sums <- c(opt.row.sums, numeric(len.opt.row.sums))
            len.opt.row.sums <- 2 * len.opt.row.sums
        }
        opt.row.sums[iter] <- opt.rs.cur.col # append it

        ## Assess convergence
        var.X.rs <- var(X.rs)
        if(iter == max.ra) {
            tol. <- tol.fun(var.X.rs, opt.row.vars.old)
            tol.reached <- FALSE
            break
        } else if(iter > n.lookback) {
            ii <- (iter %% n.lookback) + 1
            opt.row.vars.cur <- var.X.rs
            opt.row.vars.old <- opt.row.vars[ii]
            opt.row.vars[ii] <- opt.row.vars.cur
            tol. <- tol.fun(opt.row.vars.cur, opt.row.vars.old)
            if(tol. <= tol) {
                tol.reached <- TRUE
                break
            }
        } else {
            ii <- (iter %% n.lookback) + 1
            opt.row.vars[ii] <- var.X.rs
        }
    }

    ## Determine (column means of) the row(s) of the final rearranged matrix
    ## which correspond to the optimal row sum
    X.rearranged <- X # rearranged X
    opt.rs.ind <- which(X.rs == opt.rs.cur.col) # indices of rows which lead to optimal row sums
    X.rearranged.opt.row <- colMeans(X.rearranged[opt.rs.ind,, drop = FALSE]) # average over multiple rows if optimum is not unique

    ## Return
    list(bound = optim.fun(X.rs), # computed bound (\underline{s}_N or \overline{s}_N)
         tol = tol., # tolerance for the computed bound
         converged = tol.reached, # indicating whether converged
         opt.row.sums = opt.row.sums[1:iter], # the computed optimal row sums after each block rearrangement
         X.rearranged = X.rearranged, # the block-rearranged matrix
         X.rearranged.opt.row = X.rearranged.opt.row) # (averaged) rows of X.rearranged leading to the final optimal row sum
}


### 1.2 Rearrangement Algorithm ################################################

##' @title Computing lower/upper bounds for the worst/best VaR or best ES with the RA
##' @param level confidence level
##' @param qF d-list of marginal quantile functions
##' @param N number of discretization points
##' @param abstol absolute convergence tolerance (to determine convergence)
##' @param n.lookback number of column rearrangements to look back for deciding
##'        about 'convergence' (must be a number in {1,.., max.ra-1}, typically
##'        around d).
##' @param max.ra maximal number of column rearrangements
##' @param method character indicating which risk measure is approximated
##' @param sample logical indicating whether each column of the two working
##'        matrices is randomly permuted before the rearrangements begin
##' @return list containing the
##'         1) Computed lower and upper bound for (worst or best) risk measure
##'         2) The relative rearrangement gap
##'            "|(upper bound - lower bound) / upper bound|"
##'         3) Individual absolute tolerances reached (for each bound)
##'         4) 2-vector of logicals indicating whether the individual bounds reached
##'            the desired tolerances (=> convergence)
##'         5) Number of columns considered for rearrangement
##'         6) Vectors of optimal (minimal for worst VaR; maximal for best VaR;
##'            ES_alpha for best ES) row sums after each considered column rearrangement
##'         7) List of (N, d) input matrices X (for each bound)
##'         8) List of rearranged Xs (for each bound)
##'         9) List of (averaged) rows of the rearranged X leading to the optimal
##'            row sum (for each bound)
##' @author Marius Hofert
##' @note Notation is from p. 2757 in Embrechts, Puccetti, Rueschendorf (2013);
##'       variables are named according to the 'worst' VaR case.
RA <- function(level, qF, N, abstol = 0, n.lookback = length(qF),
               max.ra = Inf, method = c("worst.VaR", "best.VaR", "best.ES"),
               sample = TRUE)
{
    ## Checks and Step 1 (get N, abstol)
    stopifnot(0 < level, level < 1, is.null(abstol) || abstol >= 0,
              length(N) >= 1, N >= 2, is.logical(sample),
              is.list(qF), sapply(qF, is.function), (d <- length(qF)) >= 2, max.ra > d)
    method <- match.arg(method)

    ## Compute lower bound

    ## Step 2 (build \underline{X}^\alpha)
    p <- switch(method, # N-vector of prob. in *increasing* order
    "worst.VaR" = { # discretize the [alpha, 1) tail
        level + (1-level)*(0:(N-1))/N
    },
    "best.VaR" = { # discretize the [0, alpha) tail
        level*(0:(N-1))/N
    },
    "best.ES" = { # discretize the space [0, 1) of probabilities
        (0:(N-1))/N
    },
    stop("Wrong 'method'"))
    X.low <- sapply(qF, function(qF) qF(p))
    ## Adjust those that are -Inf (for method = "best")
    ## use level*((0+1)/2 / N) = level/(2N) instead of 0 quantile
    if(method == "best.VaR" || method == "best.ES")
        X.low[1,] <- sapply(1:d, function(j)
        if(is.infinite(X.low[1,j])) qF[[j]](level/(2*N)) else X.low[1,j])

    ## Steps 3--7 (determine \underline{X}^*)
    res.low <- if(method == "best.ES") {
        rearrange(X.low, tol = abstol, tol.type = "absolute",
                  n.lookback = n.lookback, max.ra = max.ra,
                  method = method, sample = sample,
                  is.sorted = TRUE, level = level) # need to pass 'level'
    } else {
        rearrange(X.low, tol = abstol, tol.type = "absolute",
                  n.lookback = n.lookback, max.ra = max.ra,
                  method = method, sample = sample,
                  is.sorted = TRUE)
    }

    ## Compute upper bound

    ## Step 2 (build \overline{X}^\alpha)
    p <- switch(method, # N-vector of prob. in *increasing* order
    "worst.VaR" = { # discretize the (alpha, 1] tail
        level + (1-level)*(1:N)/N
    },
    "best.VaR" = { # discretize the (0, alpha] tail
        level*(1:N)/N
    },
    "best.ES" = { # discretize the space (0, 1] of probabilities
        (1:N)/N
    }, stop("Wrong 'method'"))
    X.up <- sapply(qF, function(qF) qF(p))
    ## Adjust those that are Inf (for method = "worst")
    ## use level+(1-level)*(N-1+N)/(2*N) = level+(1-level)*(1-1/(2*N)) instead of 1 quantile
    if(method == "worst.VaR" || method == "best.ES")
        X.up[N,] <- sapply(1:d, function(j)
        if(is.infinite(X.up[N,j])) qF[[j]](level+(1-level)*(1-1/(2*N))) else X.up[N,j])

    ## Step 3--7 (determine \overline{X}^*)
    res.up <- if(method == "best.ES") {
        rearrange(X.up, tol = abstol, tol.type = "absolute",
                  n.lookback = n.lookback, max.ra = max.ra,
                  method = method, sample = sample,
                  is.sorted = TRUE, level = level) # need to pass 'level'
    } else {
        rearrange(X.up, tol = abstol, tol.type = "absolute",
                  n.lookback = n.lookback, max.ra = max.ra,
                  method = method, sample = sample,
                  is.sorted = TRUE)
    }

    ## Return
    list(bounds = c(low = res.low$bound, up = res.up$bound), # (\underline{s}_N, \overline{s}_N)
         rel.ra.gap = abs((res.up$bound-res.low$bound)/res.up$bound), # relative RA gap
         ind.abs.tol = c(low = res.low$tol, up = res.up$tol), # individual absolute tolerances
         converged = c(low = res.low$converged, up = res.up$converged), # converged?
         num.ra = c(low = length(res.low$opt.row.sums), up = length(res.up$opt.row.sums)), # number of considered column rearrangements (low, up)
         opt.row.sums = list(low = res.low$opt.row.sum, up = res.up$opt.row.sums), # optimal row sums (low, up)
         X = list(low = X.low, up = X.up), # input matrices X (low, up)
         X.rearranged = list(low = res.low$X.rearranged, up = res.up$X.rearranged), # rearranged Xs (low, up)
         X.rearranged.opt.row = # (averaged) rows of X.rearranged leading to the final optimal row sum
             list(low = res.low$X.rearranged.opt.row, up = res.up$X.rearranged.opt.row))
}


### 1.3 Adaptive Rearrangement Algorithm #######################################

##' @title Computing lower/upper bounds for the worst/best VaR or best ES with the ARA
##' @param level confidence level
##' @param qF d-list of marginal quantile functions
##' @param N.exp vector of exponents of 2 used as discretization points
##' @param reltol vector of length 2 containing the relative convergence tolerances
##'        for determining the individual and joint convergence; can also be of
##'        length 1 in which case it provides the relative joint convergence
##'        tolerance and the relative individual convergence tolerance is taken
##'        as NULL in this case.
##' @param n.lookback number of column rearrangements to look back for deciding
##'        about 'convergence' (must be a number in {1,.., max.ra-1}, typically
##'        around d).
##' @param max.ra maximal number of column rearrangements per N
##' @param method character indicating which risk measure is approximated
##' @param sample logical indicating whether each column of the two working
##'        matrices is randomly permuted before the rearrangements begin
##' @return list containing the
##'         1) Computed lower and upper bound for (worst or best) risk measure
##'         2) The relative rearrangement gap
##'            "|(upper bound - lower bound) / upper bound|"
##'         3) Relative individual tolerances reached (for each bound) and
##'            relative joint tolerance between the bounds
##'         4) 3-vector of logicals indicating whether the individual bounds and
##'            the two bounds jointly reached the desired tolerances (=> convergence)
##'         5) The number of discretization points used
##'         6) Number of columns considered for rearrangement
##'         7) Vectors of optimal (minimal for worst VaR; maximal for best VaR;
##'            ES_alpha for best ES) row sums after each considered column rearrangement
##'         8) List of (N, d) input matrices X (for each bound)
##'         9) List of rearranged Xs (for each bound)
##'        10) List of (averaged) rows of the rearranged X leading to the optimal
##'            row sum (for each bound)
##' @author Marius Hofert
ARA <- function(level, qF, N.exp = seq(8, 19, by = 1), reltol = c(0, 0.01),
                n.lookback = length(qF), max.ra = 10*length(qF),
                method = c("worst.VaR", "best.VaR", "best.ES"),
                sample = TRUE)
{
    ## Checks and Step 1 (get N, reltol)
    lreltol <- length(reltol)
    stopifnot(0 < level, level < 1, lreltol == 1 || lreltol == 2, reltol >= 0,
              length(N.exp) >= 1, N.exp >= 1, is.logical(sample),
              is.list(qF), sapply(qF, is.function), (d <- length(qF)) >= 2, max.ra > d)
    method <- match.arg(method)

    ## Determine tolerances
    itol <- if(lreltol == 2) reltol[1] else NULL # individual tolerance
    jtol <- if(lreltol == 2) reltol[2] else reltol[1] # joint tolerance

    ## Loop over N
    for(N in 2^N.exp) {

        ## Compute lower bound

        ## Step 2 (build \underline{X}^\alpha)
        p <- switch(method, # N-vector of prob. in *increasing* order
        "worst.VaR" = { # discretize the [alpha, 1) tail
            level + (1-level)*(0:(N-1))/N
        },
        "best.VaR" = { # discretize the [0, alpha) tail
            level*(0:(N-1))/N
        },
        "best.ES" = { # discretize the space [0, 1) of probabilities
            (0:(N-1))/N
        }, stop("Wrong 'method'"))
        X.low <- sapply(qF, function(qF) qF(p))
        ## Adjust those that are -Inf (for method = "best")
        ## use level*((0+1)/2 / N) = level/(2N) instead of 0 quantile
        if(method == "best.VaR" || method == "best.ES")
            X.low[1,] <- sapply(1:d, function(j)
            if(is.infinite(X.low[1,j])) qF[[j]](level/(2*N)) else X.low[1,j])

        ## Steps 3--7 (determine \underline{X}^*)
        res.low <- if(method == "best.ES") {
            rearrange(X.low, tol = itol, tol.type = "relative",
                      n.lookback = n.lookback, max.ra = max.ra,
                      method = method, sample = sample,
                      is.sorted = TRUE, level = level) # need to pass 'level'
        } else {
            rearrange(X.low, tol = itol, tol.type = "relative",
                      n.lookback = n.lookback, max.ra = max.ra,
                      method = method, sample = sample,
                      is.sorted = TRUE)
        }

        ## Compute upper bound

        ## Step 2 (build \overline{X}^\alpha)
        p <- switch(method, # N-vector of prob. in *increasing* order
        "worst.VaR" = { # discretize the (alpha, 1] tail
            level + (1-level)*(1:N)/N
        },
        "best.VaR" = { # discretize the (0, alpha] tail
            level*(1:N)/N
        },
        "best.ES" = { # discretize the space (0, 1] of probabilities
            (1:N)/N
        }, stop("Wrong 'method'"))
        X.up <- sapply(qF, function(qF) qF(p))
        ## Adjust those that are Inf (for method = "worst")
        ## use level+(1-level)*(N-1+N)/(2*N) = level+(1-level)*(1-1/(2*N)) instead of 1 quantile
        if(method == "worst.VaR" || method == "best.ES")
            X.up[N,] <- sapply(1:d, function(j)
            if(is.infinite(X.up[N,j])) qF[[j]](level+(1-level)*(1-1/(2*N))) else X.up[N,j])

        ## Step 3--7 (determine \overline{X}^*)
        res.up <- if(method == "best.ES") {
            rearrange(X.up, tol = itol, tol.type = "relative",
                      n.lookback = n.lookback, max.ra = max.ra,
                      method = method, sample = sample,
                      is.sorted = TRUE, level = level) # need to pass 'level'
        } else {
            rearrange(X.up, tol = itol, tol.type = "relative",
                      n.lookback = n.lookback, max.ra = max.ra,
                      method = method, sample = sample,
                      is.sorted = TRUE)
        }

        ## Determine (individual and joint) convergence
        joint.tol <- abs((res.low$bound - res.up$bound) / res.up$bound)
        joint.tol.reached <- joint.tol <= jtol
        if(res.low$converged && res.up$converged && joint.tol.reached) break

    }

    ## Return
    list(bounds = c(low = res.low$bound, up = res.up$bound), # (\underline{s}_N, \overline{s}_N)
         rel.ra.gap = abs((res.up$bound - res.low$bound) / res.up$bound), # relative RA gap
         tol = c(low = res.low$tol, up = res.up$tol, joint = joint.tol), # individual and joint relative tolerances
         converged = c(low = res.low$converged, up = res.up$converged,
                       joint = joint.tol.reached), # converged?
         N.used = N, # number of discretization points used
         num.ra = c(low = length(res.low$opt.row.sums), up = length(res.up$opt.row.sums)), # number of considered column rearrangements (low, up)
         opt.row.sums = list(low = res.low$opt.row.sums,
                             up = res.up$opt.row.sums), # optimal row sums (low, up) for the N used
         X = list(low = X.low, up = X.up), # input matrices X (low, up)
         X.rearranged = list(low = res.low$X.rearranged, up = res.up$X.rearranged), # rearranged Xs (low, up)
         X.rearranged.opt.row = # (averaged) rows of X.rearranged leading to the final optimal row sum
             list(low = res.low$X.rearranged.opt.row, up = res.up$X.rearranged.opt.row))
}


### 1.4 Adaptive Block Rearrangement Algorithm #################################

##' @title Computing lower/upper bounds for the worst/best VaR or best ES with the ABRA
##' @param level confidence level
##' @param qF d-list of marginal quantile functions
##' @param N.exp vector of exponents of 2 used as discretization points
##' @param absreltol vector of length 2 containing the absolute and relative
##'        convergence tolerances for determining the individual and joint
##'        convergences, respectively; can also be of length 1 in which case
##'        it provides the relative joint convergence tolerance and the
##'        absolute individual convergence tolerance is taken as 0 in this case.
##' @param n.lookback number of column rearrangements to look back for deciding
##'        about 'convergence' (must be a number in {1,.., max.ra-1}, typically
##'        around d) according to the row sum variance.
##' @param max.ra maximal number of column rearrangements per N
##' @param method character indicating which risk measure is approximated
##' @param sample logical indicating whether each column of the two working
##'        matrices is randomly permuted before the rearrangements begin
##' @return list containing the
##'         1) Computed lower and upper bound for (worst or best) risk measure
##'         2) The relative rearrangement gap
##'            "|(upper bound - lower bound) / upper bound|"
##'         3) Absolute tolerances reached (for each bound) and
##'            relative joint tolerance between the bounds
##'         4) 3-vector of logicals indicating whether the individual bounds and
##'            the two bounds jointly reached the desired tolerances (=> convergence)
##'         5) The number of discretization points used
##'         6) Number of blocks considered for rearrangement
##'         7) Vectors of optimal (minimal for worst VaR; maximal for best VaR;
##'            ES_alpha for best ES) row sums after each considered block rearrangement
##'         8) List of (N, d) input matrices X (for each bound)
##'         9) List of rearranged Xs (for each bound)
##'        10) List of (averaged) rows of the rearranged X leading to the optimal
##'            row sum (for each bound)
##' @author Marius Hofert and Martin Stefanik
ABRA <- function(level, qF, N.exp = seq(8, 19, by = 1), absreltol = c(0, 0.01),
                 n.lookback = NULL, max.ra = Inf,
                 method = c("worst.VaR", "best.VaR", "best.ES"),
                 sample = TRUE)
{
    ## Checks and Step 1 (get N, abstol)
    labsreltol <- length(absreltol)
    stopifnot(0 < level, level < 1, labsreltol == 1 || labsreltol == 2, absreltol >= 0,
              length(N.exp) >= 1, N.exp >= 1, is.logical(sample),
              is.list(qF), sapply(qF, is.function), (d <- length(qF)) >= 2)
    method <- match.arg(method)

    ## Determine tolerances
    itol <- if(labsreltol == 2) absreltol[1] else 0 # individual tolerance
    jtol <- if(labsreltol == 2) absreltol[2] else absreltol[1] # joint tolerance

    ## Function to compute the optimal n.lookback
    m.opt <- function(N, d, col.var) {
        b0 <-  1.6482635640
        b1 <- -0.0013742363
        b2 <-  0.0112121293
        b3 <-  0.0001273265
        ceiling(exp(b0 + b1 * N + b2 * d + b3 * col.var))
    }

    ## Loop over N
    for (N in 2^N.exp) {

        ## Compute lower bound

        ## Step 2 (build \underline{X}^\alpha)
        p <- switch(method, # N-vector of prob. in *increasing* order
        "worst.VaR" = { # discretize the [alpha, 1) tail
            level + (1-level)*(0:(N-1))/N
        },
        "best.VaR" = { # discretize the [0, alpha) tail
            level*(0:(N-1))/N
        },
        "best.ES" = { # discretize the space [0, 1) of probabilities
            (0:(N-1))/N
        }, stop("Wrong 'method'"))
        X.low <- sapply(qF, function(qF) qF(p))
        ## Adjust those that are -Inf (for method = "best")
        ## use level*((0+1)/2 / N) = level/(2N) instead of 0 quantile
        if(method == "best.VaR" || method == "best.ES")
            X.low[1,] <- sapply(1:d, function(j)
            if(is.infinite(X.low[1,j])) qF[[j]](level/(2*N)) else X.low[1,j])

        ## Steps 3--7 (determine \underline{X}^*)
        n.lookback.cur <- if(is.null(n.lookback))
                              m.opt(N, d = d, col.var = mean(apply(X.low, 2, var))) else n.lookback
        res.low <- if(method == "best.ES") {
            block_rearrange(X.low, tol = itol, tol.type = "absolute",
                            n.lookback = n.lookback.cur, max.ra = max.ra,
                            method = method, sample = sample,
                            level = level) # need to pass 'level'
        } else {
            block_rearrange(X.low, tol = itol, tol.type = "absolute",
                            n.lookback = n.lookback.cur, max.ra = max.ra,
                            method = method, sample = sample)
        }

        ## Computer upper bound

        ## Step 2 (build \overline{X}^\alpha)
        p <- switch(method, # N-vector of prob. in *increasing* order
        "worst.VaR" = { # discretize the (alpha, 1] tail
            level + (1-level)*(1:N)/N
        },
        "best.VaR" = { # discretize the (0, alpha] tail
            level*(1:N)/N
        },
        "best.ES" = { # discretize the space (0, 1] of probabilities
            (1:N)/N
        }, stop("Wrong 'method'"))
        X.up <- sapply(qF, function(qF) qF(p))
        ## Adjust those that are Inf (for method = "worst")
        ## use level+(1-level)*(N-1+N)/(2*N) = level+(1-level)*(1-1/(2*N)) instead of 1 quantile
        if(method == "worst.VaR" || method == "best.ES")
            X.up[N,] <- sapply(1:d, function(j)
            if(is.infinite(X.up[N,j])) qF[[j]](level+(1-level)*(1-1/(2*N))) else X.up[N,j])

        ## Step 3--7 (determine \overline{X}^*)
        n.lookback.cur <- if(is.null(n.lookback))
                              m.opt(N, d = d, col.var = mean(apply(X.up, 2, var))) else n.lookback
        res.up <- if(method == "best.ES") {
            block_rearrange(X.up, tol = itol, tol.type = "absolute",
                            n.lookback = n.lookback.cur, max.ra = max.ra,
                            method = method, sample = sample,
                            level = level) # need to pass 'level'
        } else {
            block_rearrange(X.up, tol = itol, tol.type = "absolute",
                            n.lookback = n.lookback.cur, max.ra = max.ra,
                            method = method, sample = sample)
        }

        ## Determine (individual and joint) convergence
        joint.tol <- abs((res.low$bound - res.up$bound) / res.up$bound)
        joint.tol.reached <- joint.tol <= jtol
        if(res.low$converged && res.up$converged && joint.tol.reached) break

    }

    ## Return
    list(bounds = c(low = res.low$bound, up = res.up$bound), # (\underline{s}_N, \overline{s}_N)
         rel.ra.gap = abs((res.up$bound - res.low$bound) / res.up$bound), # relative RA gap
         tol = c(low = res.low$tol, up = res.up$tol, joint = joint.tol), # individual and joint relative tolerances
         converged = c(low = res.low$converged, up = res.up$converged,
                       joint = joint.tol.reached), # converged?
         N.used = N, # number of discretization points used
         num.ra = c(low = length(res.low$opt.row.sums), up = length(res.up$opt.row.sums)), # number of considered column rearrangements (low, up)
         opt.row.sums = list(low = res.low$opt.row.sums,
                             up = res.up$opt.row.sums), # optimal row sums (low, up) for the N used
         X = list(low = X.low, up = X.up), # input matrices X (low, up)
         X.rearranged = list(low = res.low$X.rearranged, up = res.up$X.rearranged), # rearranged Xs (low, up)
         X.rearranged.opt.row = # (averaged) rows of X.rearranged leading to the final optimal row sum
             list(low = res.low$X.rearranged.opt.row, up = res.up$X.rearranged.opt.row))
}
