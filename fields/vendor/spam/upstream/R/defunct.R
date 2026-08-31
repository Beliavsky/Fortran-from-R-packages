# HEADER ####################################################
# This is file spam/R/defunct.R.                            #
# It is part of the R package spam,                         #
#  --> https://CRAN.R-project.org/package=spam              #
#  --> https://CRAN.R-project.org/package=spam64            #
#  --> https://git.math.uzh.ch/reinhard.furrer/spam         #
# by Reinhard Furrer [aut, cre], Florian Gerber [aut],      #
#    Roman Flury [aut], Daniel Gerber [ctb],                #
#    Kaspar Moesinger [ctb]                                 #
# HEADER END ################################################


validspamobject <- function(...) {
  .Defunct(new = 'validate_spam', package = 'spam', msg = "validspamobject() is defunct. Use validate_spam()")
}


