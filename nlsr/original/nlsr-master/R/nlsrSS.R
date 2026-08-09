nlsrSS <- function(formula, data){
    # wrapper to allow nlxb to self-start
    strt<-getInitial(formula, data)
    cat("suggested start=\n")
    print(strt)
    fit <- nlxb(formula, data, strt, control=list(japprox="SSJac"))
    fit
}
