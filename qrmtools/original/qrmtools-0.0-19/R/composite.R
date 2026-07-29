### Composite distributions ####################################################

##' @title Check a Composite Distribution 'distr'
##' @param distr composite distribution object
##' @return logical indicating whether the provided composite distribution
##'         is valid
##' @author Marius Hofert
check_distr <- function(distr)
{
    stopifnot(is.list(distr))
    l <- length(distr)
    if(!(l >= 2)) stop("'distr' has to be a list of length >= 2")
    vapply(distr, function(x) {
        ## Case 1: the component is a list (has to be of length 2 with either
        ##         a name string and parameter (vector) or two functions
        ##         (either [p*(),d*()] or [p*(),q*()])
        if(is.list(x)) {
            if(length(x) != 2) stop("All 'distr' sub-lists have to be of length 2")
            ## Case 1.1: first component is a name => second component is the parameter (vector)
            if(is.character(x[[1]])) {
                if(!is.numeric(x[[2]]))
                    stop("Parameters have to be provided as numeric (vectors)") # check parameter (vector)
                existsFunction(paste0("p", x[[1]]), NA) # check if p*() exists
            } else
            ## Case 1.2: first component is the function p*(), second is d*() or q*()
                is.function(x[[1]]) && is.function(x[[2]]) # check p*() and d*()
        } else {
        ## Case 2: the component is a function (in this case interpreted as p*(),
        ##         e.g., an empirical df; a density is assumed to *not* exist in
        ##         this case)
            is.function(x)
        }
    }, NA)
}

##' @title Evaluate <d/p/q/r>fun(x)
##' @param x arguments of fun()
##' @param fun distribution base name
##' @param param distributional parameters
##' @param prefix one of d/p/q/r
##' @return function values
##' @author Marius Hofert
##' @note - This even works when 'fun' has no parameters
##'       - Inspired by asCall() from 'copula'
eval_named_distr <- function(x, fun, param, prefix)
{
    fun <- paste0(prefix, fun)
    expr <- if(length(param) == 0) quote(FUN(x)) else
            as.call(c(quote(FUN), c(quote(x.), as.expression(param))))
    expr[[1]] <- as.name(fun)
    eval(expr, list(x. = x))
}

##' @title d/p/q-Functions for Composite Distributions
##' @param x evaluation points
##' @param cuts right-end cut-off points determining the partition elements
##'        (length 1 less than distr and weights)
##' @param distr list of length equal to length(cuts)+1 containing either lists
##'        or functions. If distr[i] is a...
##'        1) ... list, it either contains...
##'        1.1) ... a name string and a parameter (vector)
##'        1.2) ... the functions p*() and then d*() or q*()
##'        2) ... function, it is p*() (d*() is assumed to *not* exist in this
##'           case => NA)
##' @param weights probability weights (have to add up to 1)
##' @param method method (d/p/q)
##' @return d/p/q values of the composite distribution
##' @author Marius Hofert
composite <- function(x, cuts, distr, weights, method = c("d","p","q"))
{
    l <- length(distr)
    if(l != length(cuts) + 1)
        stop("The length of 'cuts' must be equal to length('distr')-1")
    if(any(!is.finite(cuts))) stop("'cuts' have to be finite")
    stopifnot((lx <- length(x)) > 0, l > 1, is.numeric(cuts), check_distr(distr),
              weights > 0, l == length(weights), sum(weights) == 1)
    method <- match.arg(method)
    if(method == "q") stopifnot(0 <= x, x <= 1)
    if(method == "p" || method == "q") s <- c(0, cumsum(weights))

    ## Loop over the specified buckets
    mass <- numeric(l) # probability mass in each bucket
    res <- numeric(lx)
    for(i in seq_len(l)) {

        ## Determine the mass in bucket i (P(X in (a,b]) = F(b)-F(a))
        endpts <- if(i == 1) cuts[i] else if(i == l) cuts[i-1] else c(cuts[i-1], cuts[i])
        F <- if(is.list(distr[[i]])) {
            fun <- distr[[i]][[1]]
            if(is.character(fun)) # case 1.1)
                eval_named_distr(endpts, fun = fun, param = distr[[i]][[2]], prefix = "p")
            else fun(endpts) # case 1.2)
        } else distr[[i]](endpts) # case 2)
        mass[i] <- if(i == 1) F[1] else if(i == l) 1-F[1] else F[2]-F[1]
        if(mass[i] <= 0) stop("No mass in bucket", i)

        ## evaluate d/p/q on bucket i
        ind <- if(method == "d" || method == "p") {
            if(i == 1) x <= cuts[i] else if(i == l) cuts[i-1] < x else
            cuts[i-1] < x & x <= cuts[i] # logical indicating whether x is in bucket i
        } else  # method == "q"
            (if(i == 1) s[i] <= x else s[i] < x) & x <= s[i+1] # logical indicating whether F(x) is in F(<bucket i>)
        x. <- x[ind]
        res[ind] <-
            switch(method,
                   "d" = {
                       dx <- if(is.function(distr[[i]])) rep(NA, sum(ind)) else {
                           fun <- distr[[i]][[1]]
                           if(is.character(fun))
                               eval_named_distr(x., fun = fun,
                                                param = distr[[i]][[2]], prefix = "d")
                           else distr[[i]][[2]](x.)
                       }
                       dx * weights[i] / mass[i]
                   },
                   "p" = {
                       Fx <- if(is.function(distr[[i]])) distr[[i]](x.) else {
                           fun <- distr[[i]][[1]]
                           if(is.character(fun))
                               eval_named_distr(x., fun = fun, param = distr[[i]][[2]],
                                                prefix = "p")
                           else distr[[i]][[1]](x.)
                       }
                       Flow <- if(i == 1) 0 else F[i-1]
                       s[i] + (Fx - Flow) * weights[i] / mass[i]
                   },
                   "q" = {
                       Flow <- if(i == 1) 0 else F[1]
                       p. <- Flow + mass[i] * (x.-s[i]) / weights[i]
                       if(is.function(distr[[i]]))
                           stop("No quantile function provided for distr[[",i,"]]")
                       fun <- distr[[i]][[1]]
                       if(is.character(fun))
                           eval_named_distr(p., fun = fun, param = distr[[i]][[2]],
                                            prefix = "q")
                       else distr[[i]][[2]](p.)
                   },
                   stop("Wrong 'method'"))
    }
    res
}

dcomposite <- function(x, cuts, distr, weights)
    composite(x, cuts = cuts, distr = distr, weights = weights, method = "d")

pcomposite <- function(q, cuts, distr, weights)
    composite(q, cuts = cuts, distr = distr, weights = weights, method = "p")

qcomposite <- function(p, cuts, distr, weights)
    composite(p, cuts = cuts, distr = distr, weights = weights, method = "q")

rcomposite <- function(n, cuts, distr, weights)
    composite(runif(n), cuts = cuts, distr = distr, weights = weights, method = "q")
