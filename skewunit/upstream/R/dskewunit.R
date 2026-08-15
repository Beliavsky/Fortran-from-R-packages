dskewunit<-function (x, lambda = 0, delta = 1, delta2 = 1, family1 = "asin", 
    family2 = "asin", log = FALSE) {
    if (!any(family1 == c("asin", "sbeta", "Uquad", "triang", 
        "JSB"))) 
        stop("distribution is not recognized")
    if (!any(family2 == c("asin", "sbeta", "Uquad", "triang", 
        "JSB")) & !is.null(family2)) 
        stop("distribution is not recognized")
    f <- get(paste("d", family1, sep = ""))
	if(!is.null(family2)) G <- get(paste("p", family2, sep = ""))
    if (delta <= 0) 
        stop("delta should be positive")
    if (lambda < -1 | lambda > 1) 
        stop("lambda should be between -1 and 1")
	log.lik=f(x, log = TRUE)
  if(any(family1 == c("sbeta", "JSB"))) log.lik=f(x, delta=delta, log=TRUE)
   if(!is.null(family2))
   {
    if ((family1 != "JSB" | family1 != "sbeta") & (family2 != 
        "JSB" | family2 != "sbeta")) 
        log.lik <- log(2) + G(lambda * (x - 1/2) + 1/2, log.p = TRUE) + 
            f(x, log = TRUE)
    if ((family1 == "JSB" | family1 == "sbeta") & (family2 != 
        "JSB" | family2 != "sbeta")) 
        log.lik <- log(2) + G(lambda * (x - 1/2) + 1/2, log.p = TRUE) + 
            f(x, log = TRUE, delta = delta)
    if ((family1 != "JSB" | family1 != "sbeta") & (family2 == 
        "JSB" | family2 == "sbeta")) 
        log.lik <- log(2) + G(lambda * (x - 1/2) + 1/2, log.p = TRUE, 
            delta = delta) + f(x, log = TRUE)
    if ((family1 == "JSB" | family1 == "sbeta") & (family2 == 
        "JSB" | family2 == "sbeta")) 
        log.lik <- log(2) + G(lambda * (x - 1/2) + 1/2, log.p = TRUE, 
            delta = delta2) + f(x, log = TRUE, delta = delta)
    }
    if (!log) 
        fy <- exp(log.lik)
    else fy <- log.lik
    fy <- ifelse(x <= 0 | x >= 1, 0, fy)
    fy
}
