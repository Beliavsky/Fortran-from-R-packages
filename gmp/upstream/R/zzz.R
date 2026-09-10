## we "need" S4 methods for dispatch on both (x,y)  .noGenerics <- TRUE
.conflicts.OK <- TRUE

## gmp-ify base function(s):
environment(outer) <- environment()# i.e. asNamespace("gmp")

.gmpVersion <- function() .Call(R_gmp_get_version)
gmpVersion <- function()
    numeric_version(sub("^([0-9]+\\.[0-9]+\\.[0-9]+).*","\\1", .gmpVersion()))

.onLoad <- function(libname, pkgname) {
    options("gmp:warnModMismatch" = TRUE, ## see ../man/biginteger.Rd
            "gmp:warnNoInv" = TRUE) ## ../man/add.biginteger.Rd | ../src/bigmod.cc

    ## as.big[zq]() need package dynloaded :
    gmpEnv <- parent.env(environment())
    gmpEnv$ NA_bigz_ <- as.bigz(NA)
    gmpEnv$ NA_bigq_ <- as.bigq(NA)
    invisible()
}


