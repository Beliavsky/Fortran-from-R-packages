# HEADER ####################################################
# This is file spam/R/gmult.R.                              #
# It is part of the R package spam,                         #
#  --> https://CRAN.R-project.org/package=spam              #
#  --> https://CRAN.R-project.org/package=spam64            #
#  --> https://git.math.uzh.ch/reinhard.furrer/spam         #
# by Reinhard Furrer [aut, cre], Florian Gerber [aut],      #
#    Roman Flury [aut], Daniel Gerber [ctb],                #
#    Kaspar Moesinger [ctb]                                 #
# HEADER END ################################################

gmult <- function(x, splits, fact) {
  if(!is.spam(x)) { # preserve zero entries
    x <- as.spam(x) }

  splits <- sort.int(unique(as.integer(splits)))
  if(splits[length(splits)] < max(dim(x))+1) {
    splits <- c(splits, max(dim(x))+1) }

  fact <- as.matrix(fact)
  stopifnot(dim(fact) == (length(splits)-1))

  ll <- length(x@entries)

  if(.format.spam(x)$package == "spam64" || getOption("spam.force64")) {
    SS <- .format64()
  } else {
    SS <- .format32
  }

  x@entries <- dotCall64::.C64("gmult_f",
                               SIGNATURE = c("double", SS$signature, SS$signature, SS$signature, SS$signature,
                                             "double", SS$signature, "double"),
                               a = x@entries,
                               ia = x@colindices,
                               ja = x@rowpointers,
                               na = ll,
                               splits = splits,
                               fact = fact,
                               nfact = dim(fact)[2],
                               out = dotCall64::numeric_dc(ll),
                               INTENT = c( "r", "r", "r", "r",
                                           "r", "r", "r", "rw"),
                               NAOK=getOption("spam.NAOK"),
                               PACKAGE = SS$package)$out

  return(x)
}


