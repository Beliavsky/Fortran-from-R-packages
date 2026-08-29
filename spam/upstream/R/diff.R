# HEADER ####################################################
# This is file spam/R/diff.R.                               #
# It is part of the R package spam,                         #
#  --> https://CRAN.R-project.org/package=spam              #
#  --> https://CRAN.R-project.org/package=spam64            #
#  --> https://git.math.uzh.ch/reinhard.furrer/spam         #
# by Reinhard Furrer [aut, cre], Florian Gerber [aut],      #
#    Roman Flury [aut], Daniel Gerber [ctb],                #
#    Kaspar Moesinger [ctb]                                 #
# HEADER END ################################################

     


########################################################################
diff.spam <- 
function (x, lag = 1, differences = 1, ...) 
{
    xlen <-   dim(x)[1L]
    if (length(lag) > 1L || length(differences) > 1L || lag < 
        1L || differences < 1L) 
        stop("'lag' and 'differences' must be integers >= 1")
    if (lag * differences >= xlen) {        
        return( as.spam( matrix( 0, 0L, dim(x)[2L])))
#       return( numeric( 0))        
# https://bugs.r-project.org/show_bug.cgi?id=18972
# -    return(x[0L]) # empty, but of proper mode
# +    return(if (ismat) x[0L, , drop = FALSE] else x[0L]) # empty, but of proper mode
# https://bugs.r-project.org/attachment.cgi?id=3569&action=edit
    }    
    for (i in 1L:differences){
      x <- x[(1L+lag):xlen,, drop = FALSE] - x[1L:(xlen-lag),, drop = FALSE] 
      xlen <- xlen - lag
    }
    return( x)
}

setMethod("diff","spam",diff.spam)
