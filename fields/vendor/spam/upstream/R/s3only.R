# HEADER ####################################################
# This is file spam/R/s3only.R.                             #
# It is part of the R package spam,                         #
#  --> https://CRAN.R-project.org/package=spam              #
#  --> https://CRAN.R-project.org/package=spam64            #
#  --> https://git.math.uzh.ch/reinhard.furrer/spam         #
# by Reinhard Furrer [aut, cre], Florian Gerber [aut],      #
#    Roman Flury [aut], Daniel Gerber [ctb],                #
#    Kaspar Moesinger [ctb]                                 #
# HEADER END ################################################

var.spam <- function(x, ...) {
    inefficiencywarning( "This 'var' operation may be inefficient", prod(dim(x)))
    var(as.matrix(x), ...)
}
