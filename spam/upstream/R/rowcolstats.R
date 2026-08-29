# HEADER ####################################################
# This is file spam/R/rowcolstats.R.                        #
# It is part of the R package spam,                         #
#  --> https://CRAN.R-project.org/package=spam              #
#  --> https://CRAN.R-project.org/package=spam64            #
#  --> https://git.math.uzh.ch/reinhard.furrer/spam         #
# by Reinhard Furrer [aut, cre], Florian Gerber [aut],      #
#    Roman Flury [aut], Daniel Gerber [ctb],                #
#    Kaspar Moesinger [ctb]                                 #
# HEADER END ################################################



rowSums.spam <- function(x,...) {
    ## print("1")
    if(  getOption("spam.force64") )
        SS <- .format64()
    else
        SS <- .format.spam(x)
    
    return(.C64("rowsums",
                SIGNATURE=c("double", SS$signature, SS$signature,
                    "double"),

                x@entries,
                x@rowpointers,
                x@dimension[1],

                rs = vector_dc("double", x@dimension[1]),

                INTENT=c("r", "r", "r",
                    "w"),
                NAOK = getOption("spam.NAOK"),
                PACKAGE = SS$package)$rs )
}

colSums.spam <- function(x,...) {
    ## print("2")
    if(  getOption("spam.force64") )
        SS <- .format64()
    else
        SS <- .format.spam(x)
        
    return(.C64("colsums",
                SIGNATURE=c("double", SS$signature, SS$signature, SS$signature,
                    "double"),

                x@entries,
                x@colindices,
                x@rowpointers,
                x@dimension[1],

                cs = vector_dc("double", x@dimension[2]),

                INTENT=c("r", "r", "r","r",
                    "w"),
                NAOK = getOption("spam.NAOK"),
                PACKAGE = SS$package)$cs )
}

rowMeans.spam <- function(x,...) {
        ## print("3")
    if(  getOption("spam.force64") )
        SS <- .format64()
    else
        SS <- .format.spam(x)

    return(.C64("rowmeans",
                SIGNATURE=c("double", SS$signature, SS$signature,
                    SS$signature, SS$signature, 
                    "double"),

                x@entries,
                x@rowpointers,
                x@dimension[1],
                x@dimension[2],
                getOption("spam.structurebased"),

                rm = vector_dc("double", x@dimension[1]),

                INTENT=c("r", "r", "r",
                    "r", "r",
                    "rw"),
                NAOK = getOption("spam.NAOK"),
                PACKAGE = SS$package)$rm )
}

colMeans.spam <- function(x,...) {
            ## print("4")
    if(  getOption("spam.force64") )
        SS <- .format64()
    else
        SS <- .format.spam(x)
    
    return(.C64("colmeans",
                SIGNATURE=c("double", SS$signature, SS$signature,
                    SS$signature, SS$signature, SS$signature, 
                    "double", SS$signature),

                x@entries,
                x@colindices,
                x@rowpointers,
                x@dimension[1],
                x@dimension[2],
                getOption("spam.structurebased"),

                cm = vector_dc("double", x@dimension[2]),
                vector_dc(SS$type, x@dimension[2]),

                INTENT=c("r", "r", "r",
                    "r", "r", "r",
                    "rw", "rw"),
                NAOK = getOption("spam.NAOK"),
                PACKAGE = SS$package)$cm )
}



setMethod("rowSums","spam",rowSums.spam)
setMethod("colSums","spam",colSums.spam)
setMethod("rowMeans","spam",rowMeans.spam)
setMethod("colMeans","spam",colMeans.spam)
