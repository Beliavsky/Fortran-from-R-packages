
mona <- function(x, trace.lev = 0)
{
    ## check type of input matrix
    if(!(iM <- is.matrix(x)) && !is.data.frame(x))
        stop("x must be a matrix or data frame.")
    if(!all(vapply(lapply(as.data.frame(x),
			  function(y) levels(as.factor(y))),
		   length, 1) == 2))
        stop("All variables must be binary (e.g., a factor with 2 levels, both present).")
    n <- nrow(x)
    p <- ncol(x)
    if(p < 2)
	stop("mona() needs at least p >= 2 variables (in current implementation)")
    dnx <- dimnames(x)
    ## Change levels of input matrix to {0,1, NA=2}:
    iF <- function(.) as.integer(as.factor(.))
    x <- (if(iM) apply(x, 2, iF) else vapply(x, iF, integer(n))) - 1L
    hasNA <- anyNA(x)
    if(hasNA) x[is.na(x)] <- 2L
## was
##     x <- apply(as.matrix(x), 2, factor)
##     x[x == "1"] <- "0"
##     x[x == "2"] <- "1"
##     x[is.na(x)] <- "2"
##     storage.mode(x) <- "integer"

    ## call C routine (from original fortran) -- ../src/mona.c
    res <- .C(cl_mona,
                    as.integer(n),
                    as.integer(p),
                    x = x,
                    error = as.integer(trace.lev),
                    nban = integer(n),
                    ner = integer(n),
                    integer(n),
                    lava = integer(n), # => variable numbers in every step; 0: no variable
                    integer(p))

    ## stop with a message when two many missing values:
    if(res$error != 0) {
        ## NB: Need "full simple strings below, to keep it translatable":
	switch(res$error
	       ## 1 :
	       , stop("No clustering performed, an object was found with all values missing.")
	       ## 2 :
	       , stop("No clustering performed, found variable with more than half values missing.")
	       ## 3 : never triggers because of binary check above
	       , stop("No clustering performed, a variable was found with all non missing values identical.")
	       ## 4 :
	       , stop("No clustering performed, all variables have at least one missing value.")
	       ## 5: -- cannot trigger here: already handled above
	       , stop("mona() needs at least p >= 2 variables (in current implementation)")
	       )
    }
    ##O res$x <- matrix(as.numeric(substring(res$x,
    ##O                                      1:nchar(res$x), 1:nchar(res$x))),
    ##O                      n, p)
    ## storage.mode(res$x) <- "integer" # keeping dim()
    dimnames(res$x) <- dnx
    ## add labels to Fortran output
    if(length(dnx[[2]]) != 0) {
        lava <- as.character(res$lava)
        lava[lava != "0"] <- dnx[[2]][res$lava]
        lava[lava == "0"] <- "NULL"
        res$lava <- lava
    }
    ## construct "mona" object
    structure(class = "mona",
              list(data = res$x, hasNA = hasNA, order = res$ner,
                   variable = res$lava[-1], step = res$nban[-1],
                   order.lab = if(length(dnx[[1]]) != 0) dnx[[1]][res$ner],
                   call = match.call()))
}

print.mona <- function(x, ...)
{
    ## FIXME: 1) Printing this is non-sense in the case where the data is unchanged
    ##        2) If it was changed, mona(), i.e. 'x' here should contain the info!
    d <- dim(x$data) # TODO: maybe *not* keep 'data', but keep 'dim'
    cat("mona(x, ..) fit;  x of dimension ", d[1],"x",d[2],"\n", sep="")
    if(x$hasNA) {
        cat("Because of NA's, revised data:\n")
        print(x$data, quote = FALSE, ...)
    }
    cat("Order of objects:\n")
    print(if (length(x$order.lab) != 0) x$order.lab else x$order,
          quote = FALSE, ...)
    cat("Variable used:\n")
    print(x$variable, quote = FALSE, ...)
    cat("Separation step:\n")
    print(x$step, ...)
    cat("\nAvailable components:\n")
    print(names(x), ...)
    invisible(x)
}

## FIXME: print(summary(.)) should differ from print()

summary.mona <- function(object, ...)
{
    class(object) <- "summary.mona"
    object
}

print.summary.mona <- function(x, ...)
{
    print.mona(x, ...)
    invisible(x)
}

