#' Test different methods for solving with \code{nleqslv}
#'
#' The function tests different methods and global strategies for solving a
#' system of nonlinear equations with \code{nleqslv}
#'
#' The function solves the function \code{fn} with \code{\link{nleqslv}} for
#' the specified methods and global strategies. When argument \code{Nrep} has
#' been set to a number greater than or equal to 1, repetitions of the solving
#' process are performed and the used CPU time in seconds is recorded.
#'
#' If checking a user supplied jacobian is enabled, then \code{testnslv} will
#' stop immediately when a possibly incorrect jacobian is detected.
#'
#' @param x A numeric vector with an initial guess of the root.
#' @param fn A function of \code{x} returning the function values.
#' @param jac A function to return the Jacobian for the \code{fn} function.
#' For a vector valued function \code{fn} the Jacobian must be a numeric matrix
#' of the correct dimensions.  For a scalar valued function \code{fn} the
#' \code{jac} function may return a scalar.  If not supplied numerical
#' derivatives will be used.
#' @param \dots Further arguments to be passed to \code{fn} and \code{jac} and
#' \code{\link{nleqslv}}.
#' @param method The methods to use for finding a solution.
#' @param global The global strategies to test. The argument may consist of
#' several possibly abbreviated items.
#' @param Nrep Number of repetitions to apply. Default is no repetitions.
#' @param title a description of this experiment.
#' @return \code{testnslv} returns an object of class \code{"test.nleqslv"}
#' which is a list containing the following elements \describe{
#' \item{\code{call}}{the matched call} \item{\code{out}}{ a dataframe
#' containing the results with the following columns \describe{
#' \item{\code{Method}}{method used.} \item{\code{Global}}{global strategy
#' used.} \item{\code{termcd}}{termination code of \code{nleqslv}.}
#' \item{\code{Fcnt}}{number of function evaluations used by the method and
#' global strategy.  This excludes function evaluations made when computing a
#' numerical Jacobian.} \item{\code{Jcnt}}{number of Jacobian evaluations.}
#' \item{\code{Iter}}{number of outer iterations used by the algorithm.}
#' \item{\code{Message}}{a string describing the termination code in an
#' abbreviated form.} \item{\code{Fnorm}}{square of the euclidean norm of the
#' vector of function results divided by 2.} \item{\code{cpusecs}}{CPU seconds
#' used by the requested number of repetitions (only present when argument
#' \code{Nrep} is not 0).} } } \item{\code{title}}{the description if
#' specified} } The abbreviated strings are \describe{
#' \item{\code{Fcrit}}{Convergence of function values has been achieved.}
#' \item{\code{Xcrit}}{This means that the relative distance between two
#' consecutive x-values is smaller than \code{xtol}.}
#' \item{\code{Stalled}}{The algorithm cannot find an acceptable new point.}
#' \item{\code{Maxiter}}{Iteration limit \code{maxit} exceeded.}
#' \item{\code{Illcond}}{Jacobian is too ill-conditioned.}
#' \item{\code{Singular}}{Jacobian is singular.}
#' \item{\code{BadJac}}{Jacobian is unusable.}
#' \item{\code{ERROR}}{\code{nleqslv} stopped because of a fatal error.} }
#' @section Warning: Any \code{nleqslv} error message will be displayed
#' immediately and an error for the particular combination of method and global
#' strategy will be recorded in the final dataframe.
#' @keywords nonlinear optimize
#' @examples
#'
#' dslnex <- function(x) {
#'     y <- numeric(2)
#'     y[1] <- x[1]^2 + x[2]^2 - 2
#'     y[2] <- exp(x[1]-1) + x[2]^3 - 2
#'     y
#' }
#' xstart <- c(0.5,0.5)
#' fstart <- dslnex(xstart)
#' testnslv(xstart,dslnex)
#' # this will encounter an error
#' xstart <- c(2.0,0.5)
#' fstart <- dslnex(xstart)
#' testnslv(xstart,dslnex)
#'
#' @export testnslv
testnslv <- function(x, fn, jac=NULL, ...,
                            method=c("Newton", "Broyden"),
                            global=c("cline", "qline", "gline", "pwldog", "dbldog", "hook", "none"),
                            Nrep=0L, title=NULL
                          )
{
    # utility functions
    catmsg <- function(m,g,res) {
        cat(sprintf("Error (method=%s global=%s): %s\n",m,g,attr(res,"condition")$message))
    }

    makeerrlist <- function(m,g,cpusecs=NULL) {
        if(is.null(cpusecs)) {
            z <- list(Method=m, Global=g, termcd=NA, Fcnt=NA, Jcnt=NA, Iter=NA, Message="ERROR",Fnorm=NA)
        } else {
            z <- list(Method=m, Global=g, termcd=NA, Fcnt=NA, Jcnt=NA, Iter=NA, Message="ERROR",Fnorm=NA,
                             cpusecs=cpusecs)
        }
        z
    }

    makereslist <- function(m,g,res,cpusecs=NULL) {
        fnorm <- sum(res$fvec^2)/2
        # necessary to test for termcd<0 and >6 otherwise R errors later in output
        # see R-help about switch
        if(res$termcd < 0 ) stop("User supplied jacobian most likely incorrect: cannot continue") else
        if(res$termcd > 7 ) message <- "BADCD" else
            message <- switch(res$termcd, "Fcrit", "Xcrit", "Stalled", "Maxiter", "Illcond", "Singular", "BadJac")

        if(is.null(cpusecs)) {
           z <- list(Method=m, Global=g, termcd=res$termcd, Fcnt=res$nfcnt, Jcnt=res$njcnt,
                            Iter=res$iter, Message=message, Fnorm=fnorm)
        } else {
           z <- list(Method=m, Global=g, termcd=res$termcd, Fcnt=res$nfcnt, Jcnt=res$njcnt,
                            Iter=res$iter, Message=message, Fnorm=fnorm, cpusecs=cpusecs)
        }
        z
    }

    methods <- match.arg(method, c("Newton", "Broyden"), several.ok=TRUE)
    globals <- match.arg(global, c("cline", "qline", "gline", "pwldog", "dbldog", "hook", "none"), several.ok=TRUE)

    my.call <- match.call()
    reslist <- vector("list", length(methods)*length(globals))

    # use try to avoid process stopping for Jacobian with non-finite values
    # if that happens, go to next method/global
    # avoidable fatal user errors will also lead to useless next method/global
    idx <- 1
    for(m in methods)
        for(g in globals) {
            if( Nrep >= 1) {
                mytime <- system.time( for(k in seq_len(Nrep)) {
                                res <- try(nleqslv(x, fn, jac, ..., method=m, global=g), silent=TRUE)
                                if(inherits(res,"try-error")) break
                            }, gcFirst = FALSE)
                cpus <- mytime[3]
            } else {
                res <- try(nleqslv(x, fn, jac, ..., method=m, global=g),silent=TRUE)
                cpus <- NULL
            }
            if(inherits(res,"try-error")) {
                catmsg(m,g,res)
                z <- makeerrlist(m,g,cpus)
            } else {
                z <- makereslist(m,g,res,cpus)
            }
            reslist[[idx]] <- z
            idx <- idx+1
        }

# from http://stackoverflow.com/questions/4512465/what-is-the-most-efficient-way-to-cast-a-list-as-a-data-frame?rq=1

    ## @Martin Morgan's Map() sapply() solution:
    f <- function(x) function(i) sapply(x, `[[`, i)
    z <- as.data.frame(Map(f(reslist), names(reslist[[1]])), stringsAsFactors=FALSE)

    res <- list()
    res$out <- z
    res$call <- my.call
    res$title <- title
    class(res) <- "test.nleqslv"
    res
}



#' Printing the result of \code{testnslv}
#'
#' Print method for \code{test.nleqslv} objects.
#'
#' This is the \code{print} method for objects inheriting from class
#' \code{test.nleqslv}. It prints the call to \code{testnslv} followed by the
#' description of the experiment (if the \code{title} argument was specified in
#' the call to \code{testnslv}) and the dataframe containing the results of
#' \code{testnslv}.
#'
#' @aliases print print.test.nleqslv
#' @param x a \code{test.nleqslv} object
#' @param digits specifies the minimum number of significant digits to be
#' printed in values.
#' @param width.cutoff integer passed to \code{\link{deparse}} which sets the
#' cutoff at which line-breaking is tried.
#' @param \dots additional arguments to \code{print}.
#' @return It returns the object \code{x} invisibly.
#' @keywords print
#' @export
#' @examples
#'
#' dslnex <- function(x) {
#'     y <- numeric(2)
#'     y[1] <- x[1]^2 + x[2]^2 - 2
#'     y[2] <- exp(x[1]-1) + x[2]^3 - 2
#'     y
#' }
#' xstart <- c(1.5,0.5)
#' fstart <- dslnex(xstart)
#' z <- testnslv(xstart,dslnex)
#' print(z)
#'
print.test.nleqslv <- function(x, digits=4, width.cutoff=45L, ...) {
    if(!inherits(x, "test.nleqslv"))
        stop("method is only for test.nleqslv objects")

    # calculate total number of function evaluations if numeric jacobian used

    cat("Call:\n",paste0(deparse(x$call, width.cutoff=width.cutoff), collapse = "\n"), "\n\n", sep = "")
    if(is.null(x$title)) cat("Results:\n") else cat("Results: ",x$title,"\n", sep="")
    print(x$out, digits=digits,...)
    invisible(x)
}
