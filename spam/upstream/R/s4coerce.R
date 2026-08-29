# HEADER ####################################################
# This is file spam/R/s4coerce.R.                           #
# It is part of the R package spam,                         #
#  --> https://CRAN.R-project.org/package=spam              #
#  --> https://CRAN.R-project.org/package=spam64            #
#  --> https://git.math.uzh.ch/reinhard.furrer/spam         #
# by Reinhard Furrer [aut, cre], Florian Gerber [aut],      #
#    Roman Flury [aut], Daniel Gerber [ctb],                #
#    Kaspar Moesinger [ctb]                                 #
# HEADER END ################################################



# a few coercions that make sense...

#      showMethods(coerce)

setAs("spam","logical", def=function(from) {
    if(getOption("spam.structurebased")) {
        return( as.logical(from@entries))     
    }else{
        inefficiencywarning( gettextf("This operation may be inefficient"), prod(dim(from)))
        return( as.logical(as.matrix(from)))
    }})

setAs("spam","vector", def=function(from) {
    if(getOption("spam.structurebased")) {
        return( as.vector(from@entries))     
    }else{
        inefficiencywarning( gettextf("This operation may be inefficient"), prod(dim(from)))
        return( as.vector(as.matrix(from)))
    }})

setAs("spam","integer", def=function(from) {
    if(getOption("spam.structurebased")) {
        return( as.integer(from@entries))     
    }else{
        inefficiencywarning( gettextf("This operation may be inefficient"), prod(dim(from)))
        return( as.integer(as.matrix(from)))
    }})

setAs("spam","matrix", def=function(from) {
    inefficiencywarning( gettextf("This operation may be inefficient"), prod(dim(from)))
    return( as.logical(as.matrix(from)))
})
setAs("spam","list", def=function(from) {
    return( triplet(from))
})
