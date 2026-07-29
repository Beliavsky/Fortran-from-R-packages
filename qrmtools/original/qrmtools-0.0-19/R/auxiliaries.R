### Auxiliary tools ############################################################

##' @title Constructing Sequences (in Log-Scale)
##' @param from see ?seq
##' @param to see ?seq
##' @param length.out see ?seq
##' @param log logical indicating whether a logarithmic sequence is to be
##'        constructed
##' @param base base of the logarithm
##' @param ... see ?seq
##' @return numeric sequence
##' @author Marius Hofert
seq2 <- function(from, to, length.out = NULL, log = FALSE, base = 10, ...)
{
    stopifnot(is.logical(log), base > 0) # everything else checked in seq.default()
    if(log) {
        stopifnot(from > 0, to > 0)
        base^seq(from = log(from, base = base), to = log(to, base = base),
                 length.out = length.out, ...)
    } else seq(from = from, to = to, length.out = length.out, ...)
}
