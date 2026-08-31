# HEADER ####################################################
# This is file spam/R/rep_len64.R.                          #
# It is part of the R package spam,                         #
#  --> https://CRAN.R-project.org/package=spam              #
#  --> https://CRAN.R-project.org/package=spam64            #
#  --> https://git.math.uzh.ch/reinhard.furrer/spam         #
# by Reinhard Furrer [aut, cre], Florian Gerber [aut],      #
#    Roman Flury [aut], Daniel Gerber [ctb],                #
#    Kaspar Moesinger [ctb]                                 #
# HEADER END ################################################

rep_len64 <- function(x, length.out, NAOK = getOption("spam.NAOK")){
    if(getOption("spam.force64") || length.out > 2147483647){
        .format64()
        return(.C64("rep_len64_c",
                    SIGNATURE = c("double", "int64", "int64", "double"),
                    
                    x = x,
                    lx = length(x),
                    length = length.out,
                    
                    out = numeric_dc(length.out),
                    
                    INTENT = c("r","r","r","w"),
                    NAOK = NAOK,
                    PACKAGE = "spam64"
                    )$out)
    }
    return(rep_len(x = x, length.out = length.out))
}
