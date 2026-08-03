#' ecd: A package for the stable lambda distribution family.
#'
#' The ecd package provides the core classes and functions for the stable lambda distribution family.
#' The stable lambda distribution is implemented in \code{\link{dsl}} section.
#' The lambda distribution uses the \code{ecld} namespace. SGED is considered part of ecld.
#' (See \code{\link{ecld-class}} for definition.)
#' The original elliptic lambda distribution uses the generic methods or \code{ecd} namespace. 
#' (See \code{\link{ecd-class}} for definition.)
#' The option pricing API uses the \code{ecop} namespace.
#' (See \code{\link{ecop-class}} for definition.)
#' Most helper utilities are named under either \code{ecd} or \code{ecld}.
#' 
#' @author Stephen H-T. Lihn
#'
#' @docType package
#' @name ecd-package
#' @import xts methods polynom graphics moments stabledist parallel yaml RSQLite
#' 
#' @seealso The two main classes are \code{\link{ecd-class}} and \code{\link{ecld-class}}
#' 
NULL

# Some areas of this package require multi-core capability
cores <- switch( Sys.info()[['sysname']],
    Windows = 1,
    Linux   = parallel::detectCores(),
    Darwin  = parallel::detectCores(),
    parallel::detectCores()
    )

if (is.null(getOption("mc.cores"))) {
    options("mc.cores"=cores)
}

# MPFR default settings

if (is.null(getOption("ecd.precBits"))) {
    options("ecd.precBits"=120L)
}

# MPFR default Inf conversion, number of sigma as replacement for +/- Inf
# for integrateR and imgf

.ecd.mpfr.N.sigma <- 300

# end

