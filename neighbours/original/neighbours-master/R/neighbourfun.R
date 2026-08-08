neighbourfun <- function(min = 0,
                         max = 1,
                         kmin = NULL,
                         kmax = NULL,
                         stepsize,
                         sum = TRUE,
                         random = TRUE,
                         update = FALSE,
                         type = "numeric",
                         active = TRUE,
                         length = NULL,
                         A = NULL, ...) {

    wmin <- min
    wmax <- max
    budget <- sum

    .sub <- function(e, env) {
        ## replace a subexpression within e by another
        ## expression, bound in env. Examples:
        ##    .sub(quote(a + b), list(a = function(x) x + 1))
        ##    ## => function(x) x + 1 + b
        ##    .sub(quote(a + b), list(a = quote((function(x) x + 1))))
        ##    ## => (function(x) x + 1) + b

        expr <- substitute(substitute(e2, env),
                           env = list(e2 = e))
        eval(expr)
    }

    if (type == "numeric") {

        if (isTRUE(budget) || is.numeric(budget) && length(budget) == 1L) {

            ## NOTE 'budget' is a scalar, default is 1. Two
            ##      elements are selected: one is increased,
            ##      the other is decreased.

            .body <- quote({

                ## NOTE If wmin/wmax are -Inf/Inf, then all
                ##      elements could be selected.
                decrease <- which((x > wmin)[active])
                increase  <- which((x < wmax)[active])

                ## NOTE More than two elements could be selected.
                ##      It might also be necessary to check if
                ##      enough elements are available:
                ##      min(size, length(decrease), length(icnrease))
                i <- decrease[sample.int(length(decrease), size = 1L)]
                j <- increase[sample.int(length(increase), size = 1L)]
                stepsize <- .stepsize
                stepsize <- min(x[i] - wmin[i], wmax[j] - x[j], stepsize)
                x[i] <- x[i] - stepsize
                x[j] <- x[j] + stepsize
                x
            })
            names(.body)[10] <- c("return value")

            if (length(wmin) == 1) {

                ## if 'wmin' is a scalar, there is no
                ## need for subsetting

                ## .body[[7]]
                ## ==>  stepsize <- min(x[i] - wmin[i], wmax[j] - x[j], stepsize)
                ## .body[[7]][[3]]
                ## ==>  min(x[i] - wmin[i], wmax[j] - x[j], stepsize)
                ## .body[[7]][[3]][[2]]
                ## ==>  x[i] - wmin[i]
                ## .body[[7]][[3]][[2]][[3]]
                ## ==>  wmin[i]
                .body[[7]][[3]][[2]][[3]] <- wmin
            }

            if (length(wmax) == 1) {

                ## if 'wmax' is a scalar, there is no
                ## need for subsetting

                ## .body[[7]]
                ## ==>  stepsize <- min(x[i] - wmin[i], wmax[j] - x[j], stepsize)
                ## .body[[7]][[3]]
                ## ==>  min(x[i] - wmin[i], wmax[j] - x[j], stepsize)
                ## .body[[7]][[3]][[3]]
                ## ==>  wmax[j] - x[j]
                ## .body[[7]][[3]][[3]][[2]]
                ## ==>  wmax[j]
                .body[[7]][[3]][[3]][[2]] <- wmax
            }

        } else if (isFALSE(budget) || is.numeric(budget) && length(budget) == 2L) {

            ## budget is of length 2 or not defined

            .body <- quote({
                stepsize <- sample(c(-1, 1), size = 1) * .stepsize
                if (stepsize < 0) {
                    decrease <- which((x > wmin)[active])
                    i <- decrease[sample.int(length(decrease), size = 1L)]
                    stepsize <- max(wmin - x[i], stepsize, budget[1L] - sum(x))
                    ##                                     ^^^^^^^^^^^^^^^^^^^
                } else {
                    increase  <- which((x < wmax)[active])
                    i <- increase[sample.int(length(increase),  size = 1L)]
                    stepsize <- min(wmax - x[i], stepsize, budget[2L] - sum(x))
                    ##                                     ^^^^^^^^^^^^^^^^^^^
                }
                x[i] <- x[i] + stepsize
                x
            })
            if (isFALSE(budget)) {
                .body[[3]][[3]][[4]][[3]][[4]] <- NULL
                .body[[3]][[4]][[4]][[3]][[4]] <- NULL
            }
        } else {
            stop("budget must be logical or numeric")
        }


        ## [random]
        .body <- .sub(.body,
                      list(.stepsize =
                               if (random)
                                   quote(stepsize * runif(1L))
                               else
                                   stepsize))


        ## [wmin/wmax]
        ## if (length(wmin) > 1L || length(wmax) > 1L) {

        ##     ## wmin or wmax or both have length > 1
        ##     if (length(wmin) == 1L)
        ##         wmin <- rep(wmin, length(wmax))
        ##     if (length(wmax) == 1L)
        ##         wmax <- rep(wmax, length(wmin))

        ##     if (!isTRUE(active))
        ##         .body <- .sub(.body, list(wmin = quote(wmin[active]),
        ##                                   wmax = quote(wmax[active])))

        ## } else if (!isFALSE(budget) && length(budget) == 1L) {

        ##     ## wmin and wmax have length 1: no subsetting
        ##     .body[[7L]] <- quote(
        ##         stepsize <- min(x[i] - wmin, wmax - x[j], stepsize))
        ## }

        ## ## [active]
        ## if (!isTRUE(active)) {
        ##     .body <- .sub(.body, list(x = quote(x[active])))
        ##     .body[[length(.body)]] <- quote(x)
        ## }


        ## [update]
        if (is.character(update) && update == "Ax") {
            .body[[10L]] <- quote(
                attr(x, "Ax") <- attr(x, "Ax") +
                    A[, c(i, j)] %*% c(-stepsize, stepsize))
            .body[[11L]] <- quote(x)
        }



        ans <- function(x, ...) {}
        body(ans) <- .body
        return(ans)
    }

    if (type == "permute") {

        if (missing(stepsize) || stepsize == 1L)
            stepsize <- 2L

        if (stepsize == 2L) {
            ans <- function(x, ...) {
                ii <- sample.int(length(x), stepsize)
                x[ii] <- x[ii[2:1]]
                x
            }
        } else {
            ans <- function(x, ...) {
                ii <- sample.int(length(x), size = stepsize)
                x[ii] <- x[sample(ii)]
                x
            }
        }
        return(ans)
    }

    if (type == "logical") {

        if (missing(stepsize))
            stepsize <- 1L

        .body <- quote({
            i <- sample.int(length, stepsize)
            x[i] <- !x[i]
            x
        })
        if (is.null(kmin) && is.null(kmax)) {

            if (is.null(length)) {
                .body[[2]] <- .sub(.body[[2]],
                                   list(length = quote(length(x)),
                                        stepsize = stepsize))
            } else if (!isTRUE(active)) {
                .body[[2]] <- .sub(.body[[2]],
                                   list(length = sum(active)))
                .body <- .sub(.body,
                              list(x = quote(x[active])))
            }
            ans <- function(x, ...) {}
            body(ans) <- .body
            return(ans)


        } else if (!is.null(kmin) && !is.null(kmax) && kmin == kmax) {

            ## logical with constant number of TRUE values

            if (!is.null(active)) {

                function(x, ...) {
                    xx <- x[active]
                    true  <- which( xx)
                    false <- which(!xx)
                    xx[true [sample.int(length( true), size = stepsize)]] <- FALSE
                    xx[false[sample.int(length(false), size = stepsize)]] <- TRUE
                    x[active] <- xx
                    x
                }


            } else {
                function(x, ...) {
                    true  <- which( x)
                    false <- which(!x)
                    x[true [sample.int(length( true), size = stepsize)]] <- FALSE
                    x[false[sample.int(length(false), size = stepsize)]] <- TRUE
                    x
                }
            }

        } else if (!is.null(kmin) && !is.null(kmax) && kmin < kmax) {

            function(x, ...) {
                true  <- which( x)
                false <- which(!x)
                n.true <- length(true)
                if (n.true == kmax) {
                    x[true[sample.int(length(true), size = stepsize)]] <- FALSE
                } else if (n.true > kmin) {
                    i <- sample.int(length(x), size = stepsize)
                    x[i] <- !x[i]
                } else {
                    x[false[sample.int(length(false), size = stepsize)]] <- TRUE
                }
                x
            }

        }

    } else if (type == "5/10/40") {
        wmax  <- 0.05
        wmax2 <- 0.1
        max.sumL <- 0.4
        if (is.null(kmax))
            kmax <- 33
        function(x, ...) {
            k <- sum(abs(x) > 0)
            eps <- runif(1)*0.5/100

            to_sell <- x > 0
            to_buy  <- if (k == kmax)
                           x > 0 & x < wmax2
                       else
                           x < wmax2
            to_sell <- which(to_sell)
            to_buy  <- which(to_buy)
            sumL <- sum(x[x > wmax])

            i <- to_sell[sample.int(length(to_sell), size = 1L)]
            j <- to_buy [sample.int(length(to_buy),  size = 1L)]
            eps <- if (x[j] < wmax)
                       min(eps, wmax  - x[j], x[i])
                   else if (x[j] == wmax)
                       min(eps, wmax2 - x[j], x[i], max(0, max.sumL - sumL - x[j]))
                   else
                       min(eps, wmax2 - x[j], x[i], max(0, max.sumL - sumL))
            x[i] <- x[i] - eps
            x[j] <- x[j] + eps
            x
        }
    } else
        stop("no matches")
}

neighborfun <- neighbourfun

compare_vectors <- function(...,
                            sep = "",
                            diff.char = "|") {

    ## TODO make arguments
    compare1 <- TRUE  ## compare all solutions with the 1st
    rows <- TRUE      ## print rows
    FALSE.TRUE <- c("0", "1")


    vecs <- list(...)
    if (!is.logical(vecs[[1]]))
        stop("currently only supported for logical vectors")
    len.x <- length(vecs)
    if (length(unique(lengths(vecs))) != 1L)
        stop("vectors have different lengths")
    if (mode(vecs[[1L]]) == "logical") {
        d <- numeric(length(vecs) - 1L)
        cat(as.integer(vecs[[1]]), "\n", sep = "")
        if (len.x > 1L) {
            for (i in 2:length(vecs)) {
                if (nchar(diff.char))
                    cat(ifelse(vecs[[i - 1L]] == vecs[[i]], " ", diff.char),
                        "\n", sep = "")
                cat(as.integer(vecs[[i]]), "\n", sep = "")
                d[i - 1L] <- sum(vecs[[i - 1L]] != vecs[[i]])
            }
            if (len.x == 2L)
                message("The vectors differ in  ", d, "  place",
                        if (d != 1) "s", ".")
        }
    }
    invisible(d)
}

random_vector <- function(length,
                          min = 0,
                          max = 1,
                          kmin = NULL,
                          kmax = NULL,
                          sum = NULL,
                          type = "numeric",
                          n = 1,
                          ...) {
    ans <- NULL
    if (type == "logical") {
        if (missing(length))
            stop(sQuote("length"), " missing")
        if (is.null(kmin))
            kmin <- 0
        if (is.null(kmax))
            kmax <- length
        stopifnot(kmin <= kmax)
        ans <- array(logical(length*n), dim = c(length, n))
        for (j in seq_len(n)) {

            if (kmin == kmax)
                k <- kmin
            else
                k <- sample(seq(from = kmin, to = kmax), size = 1)

            ## ans <- logical(length)
            i <- sample(length, size = k)
            ans[i, j] <- TRUE
        }
        if (n == 1L)
            dim(ans) <- NULL

    } else if (type == "numeric") {
        if (missing(length))
            stop(sQuote("length"), " missing")
        stopifnot(min <= max)

        if (is.null(kmin) && is.null(kmax)) {
            if (n == 1) {
                ans <- runif(length, min = min, max = max)
            } else {
                ans <- runif(length*n, min = min, max = max)
                dim(ans) <- c(length, n)
            }
        } else {
            if (is.null(kmin))
                kmin <- 0
            if (is.null(kmax))
                kmax <- length
            ans <- runif(length*n, min = min, max = max)
            dim(ans) <- c(length, n)
            for (j in seq_len(n)) {
                if (kmin == kmax)
                    k <- length - kmin
                else
                    k <- length - sample(seq(from = kmin, to = kmax), size = 1)
                i <- sample(length, size = k)
                ans[i, j] <- 0
            }
        }
    }
    ans
}
