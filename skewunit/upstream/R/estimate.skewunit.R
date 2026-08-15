estimate.skewunit<-function (x, family1 = "asin", family2 = "asin", est.var = TRUE) 
{
    if (any(x < 0 | x > 1)) 
        stop("x should be between 0 and 1")
    if (!any(family1 == c("asin", "sbeta", "Uquad", "triang", 
        "JSB"))) 
        stop("distribution is not recognized")
    if (!any(family2 == c("asin", "sbeta", "Uquad", "triang", 
        "JSB")) & !is.null(family2)) 
        stop("distribution is not recognized")
    if(!is.null(family2))
	{
    if ((family1 != "JSB" & family1 != "sbeta") & (family2 != 
        "JSB" & family2 != "sbeta")) {
        case = 1
    }
    if (((family1 == "JSB" | family1 == "sbeta") & (family2 != 
        "JSB" & family2 != "sbeta")) | ((family1 != "JSB" & family1 != 
        "sbeta") & (family2 == "JSB" | family2 == "sbeta"))) {
        case = 2
    }
    if ((family1 == "JSB" | family1 == "sbeta") & (family2 == 
        "JSB" | family2 == "sbeta")) {
        case = 3
    }
	}
    if(any(family1!=c("sbeta","JSB")) & is.null(family2)){case=4}
    if(any(family1==c("sbeta","JSB")) & is.null(family2)){case=5}
    llike.skewunit <- function(theta, x, family1, family2, mod.param = FALSE, 
        case = 1) {
	if(family1=="Uquad")
	{
		if(length(which(x==0.5))>0) x=x[-which(x==0.5)]
	}
        if (case == 1) {
            lambda <- theta[1]/sqrt(1 + theta[1]^2)
            if (!mod.param) 
                lambda <- theta[1]
            logf <- dskewunit(x, lambda = lambda, family1 = family1, 
                family2 = family2, log = TRUE)
        }
        if (case == 2) {
            lambda <- theta[1]/sqrt(1 + theta[1]^2)
            delta <- exp(theta[2])
            if (!mod.param) {
                lambda <- theta[1]
                delta <- theta[2]
            }
            logf <- dskewunit(x, lambda = lambda, delta = delta, 
                family1 = family1, family2 = family2, log = TRUE)
        }
        if (case == 3) {
            lambda <- theta[1]/sqrt(1 + theta[1]^2)
            delta1 <- exp(theta[2])
            delta2 <- exp(theta[3])
            if (!mod.param) {
                lambda <- theta[1]
                delta1 <- theta[2]
                delta2 <- theta[3]
            }
            logf <- dskewunit(x, lambda = lambda, delta = delta1, 
                delta2 = delta2, family1 = family1, family2 = family2, 
                log = TRUE)
        }
	if(case==4){
	logf=dskewunit(x, family1 = family1, family2 = family2, log = TRUE)
	}
	if(case==5){
	delta <- exp(theta[1])
	if(!mod.param) delta=theta[1]
	logf=dskewunit(x, delta=delta, family1 = family1, family2 = family2, log = TRUE)
	}
        -sum(logf)
    }
	if(case==4)
	{
	aux1=aux2=list(convergence=0, value=llike.skewunit(1, x, family1, family2, mod.param = FALSE, case=case))
	}
    if (case == 1) {
        aux1 = suppressWarnings(optim(rep(-0.2, case), llike.skewunit, 
            x = x, family1 = family1, family2 = family2, mod.param = TRUE, 
            case = case, method = "Brent", control = list(maxit = 1e+05), 
            lower = -100, upper = 100))
        aux2 = suppressWarnings(optim(rep(0.2, case), llike.skewunit, 
            x = x, family1 = family1, family2 = family2, mod.param = TRUE, 
            case = case, method = "Brent", control = list(maxit = 1e+05), 
            lower = -100, upper = 100))
    }
    if (case == 2 | case==3) {
        aux1 = suppressWarnings(optim(rep(-0.2, case), llike.skewunit, 
            x = x, family1 = family1, family2 = family2, mod.param = TRUE, 
            case = case, method = "Nelder-Mead", control = list(maxit = 1e+05)))
        aux2 = suppressWarnings(optim(rep(0.2, case), llike.skewunit, 
            x = x, family1 = family1, family2 = family2, mod.param = TRUE, 
            case = case, method = "Nelder-Mead", control = list(maxit = 1e+05)))
    }
    if (case == 5) {
        aux1 = suppressWarnings(optim(0, llike.skewunit, 
            x = x, family1 = family1, family2 = family2, mod.param = TRUE, 
            case = case, method = "Brent", control = list(maxit = 1e+05),  lower = -100, upper = 100))
        aux2 = suppressWarnings(optim(rep(0.2, case), llike.skewunit, 
            x = x, family1 = family1, family2 = family2, mod.param = TRUE, 
            case = case, method = "Nelder-Mead", control = list(maxit = 1e+05)))
    }
	param=NULL
    if (aux1$convergence == 0) 
        aux = aux1
    if (aux2$convergence == 0 & -aux2$value > -aux1$value) 
        aux = aux2
    if (case == 1) 
        param = aux$par[1]/sqrt(1 + aux$par[1]^2)
    if (case ==2 | case==3) 
        param = c(aux$par[1]/sqrt(1 + aux$par[1]^2), exp(aux$par[-1]))
    if (case == 5) 
        param = exp(aux1$par[1])
	if(case!=4)	param=matrix(param, ncol=1)
    if (est.var & case!=4) {
    param = matrix(param, ncol = 1)
    colnames(param) = c("estimate")
    se = c()
        test = suppressWarnings(try(solve(hessian(llike.skewunit, 
            x0 = param, x = x, family1 = family1, family2 = family2, 
            mod.param = FALSE, case = case)), silent = TRUE))
	if(!grepl("Error",test)[1])
	{
        if (is.numeric(test) & min(diag(test)) > 0) {
            se = sqrt(diag(test))
            param <- cbind(param, se)
            colnames(param) <- c("estimate", "s.e.")
        }
	}
    }
    if (case == 1) 
        rownames(param) = c("lambda")
    if (case == 2) 
        rownames(param) = c("lambda", "delta")
    if (case == 3) 
        rownames(param) = c("lambda", "delta1", "delta2")
    if (case == 5) 
        rownames(param) = c("delta")
	if(case==4) case=0
	if(case==5) case=1
    logvero <- -aux$value
    AIC = -2 * logvero + 2 * case
    BIC = -2 * logvero + log(length(x)) * case
    object.out <- list(coefficients = param, logLik = logvero, 
        AIC = AIC, BIC = BIC, family1 = family1, family2 = family2)
    class(object.out) <- "skewunit"
    object.out
}
