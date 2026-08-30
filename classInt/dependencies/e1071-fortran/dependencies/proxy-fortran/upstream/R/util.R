
## wrappers for class dist
##
## note that all type checking and coercing
## is now done in C, as well as handling of
## attributes.
##
## fixme: create generic functions?
##
## ceeboo 2007

dim.dist <-
function(x)
    rep.int(attr(x, "Size"), 2)         # works with nrow and ncol

dimnames.dist <-
names.dist <-
function(x)
    rep(list(attr(x, "Labels")), 2)

"dimnames<-.dist" <-
"names<-.dist" <-
function(x, value)
{
    if (is.list(value))
        value <- value[[1]]
    if (is.null(value))
        attr(x, "Labels") <- NULL
    else {
        if (length(value) != attr(x, "Size"))
            stop("Length of labels does not conform Size.")
        attr(x, "Labels") <- as.character(value)
    }
    x
}

row.dist <-
function(x)
    .Call(R_row_dist, x, FALSE)

col.dist <-
function(x)
    .Call(R_row_dist, x, TRUE)

##

subset.dist <-
"[[.dist" <-
function(x, subset, ...)
{
    if (missing(subset))
        return(x)
    .Call(R_subset_dist, x, unique(subset))
}

##

rowSums.dist <-
colSums.dist <-
function(x, na.rm = FALSE)
    .Call(R_rowSums_dist, x, na.rm)

##

rowMeans.dist <-
colMeans.dist <-
function(x, na.rm = FALSE, diag = TRUE)
{
    if (!is.logical(diag))
        stop("'diag' not of type logical")
    s <- rowSums.dist(x, na.rm)
    if (na.rm) {
        x[!(is.na(x) | is.nan(x))] <- 1
        s / (rowSums.dist(x, na.rm) + (diag == TRUE))
    } else
        s / (length(s) - (diag == FALSE))
}

###

"[.dist" <-
function(x, i, j, drop = FALSE)
{
    na = nargs() - !missing(drop)
    if (na < 3L)
        as.vector(x)[i]
    else
        as.matrix(x)[i, j, drop = drop]
}
