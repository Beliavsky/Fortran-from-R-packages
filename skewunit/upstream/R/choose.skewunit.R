choose.skewunit = function (x, criteria = "AIC") 
{
    fit <- c()
    if (criteria == "AIC") {
        fit[1] = estimate.skewunit(x, family1 = "asin", family2 = "asin", 
            est.var = FALSE)$AIC
        fit[2] = estimate.skewunit(x, family1 = "asin", family2 = "Uquad", 
            est.var = FALSE)$AIC
        fit[3] = estimate.skewunit(x, family1 = "asin", family2 = "triang", 
            est.var = FALSE)$AIC
        fit[4] = estimate.skewunit(x, family1 = "asin", family2 = "JSB", 
            est.var = FALSE)$AIC
        fit[5] = estimate.skewunit(x, family1 = "asin", family2 = "sbeta", 
            est.var = FALSE)$AIC
        fit[6] = estimate.skewunit(x, family1 = "Uquad", family2 = "asin", 
            est.var = FALSE)$AIC
        fit[7] = estimate.skewunit(x, family1 = "Uquad", family2 = "Uquad", 
            est.var = FALSE)$AIC
        fit[8] = estimate.skewunit(x, family1 = "Uquad", family2 = "triang", 
            est.var = FALSE)$AIC
        fit[9] = estimate.skewunit(x, family1 = "Uquad", family2 = "JSB", 
            est.var = FALSE)$AIC
        fit[10] = estimate.skewunit(x, family1 = "Uquad", family2 = "sbeta", 
            est.var = FALSE)$AIC
        fit[11] = estimate.skewunit(x, family1 = "triang", family2 = "asin", 
            est.var = FALSE)$AIC
        fit[12] = estimate.skewunit(x, family1 = "triang", family2 = "Uquad", 
            est.var = FALSE)$AIC
        fit[13] = estimate.skewunit(x, family1 = "triang", family2 = "triang", 
            est.var = FALSE)$AIC
        fit[14] = estimate.skewunit(x, family1 = "triang", family2 = "JSB", 
            est.var = FALSE)$AIC
        fit[15] = estimate.skewunit(x, family1 = "triang", family2 = "sbeta", 
            est.var = FALSE)$AIC
        fit[16] = estimate.skewunit(x, family1 = "JSB", family2 = "asin", 
            est.var = FALSE)$AIC
        fit[17] = estimate.skewunit(x, family1 = "JSB", family2 = "Uquad", 
            est.var = FALSE)$AIC
        fit[18] = estimate.skewunit(x, family1 = "JSB", family2 = "triang", 
            est.var = FALSE)$AIC
        fit[19] = estimate.skewunit(x, family1 = "JSB", family2 = "JSB", 
            est.var = FALSE)$AIC
        fit[20] = estimate.skewunit(x, family1 = "JSB", family2 = "sbeta", 
            est.var = FALSE)$AIC
        fit[21] = estimate.skewunit(x, family1 = "sbeta", family2 = "asin", 
            est.var = FALSE)$AIC
        fit[22] = estimate.skewunit(x, family1 = "sbeta", family2 = "Uquad", 
            est.var = FALSE)$AIC
        fit[23] = estimate.skewunit(x, family1 = "sbeta", family2 = "triang", 
            est.var = FALSE)$AIC
        fit[24] = estimate.skewunit(x, family1 = "sbeta", family2 = "JSB", 
            est.var = FALSE)$AIC
        fit[25] = estimate.skewunit(x, family1 = "sbeta", family2 = "sbeta", 
            est.var = FALSE)$AIC
        fit[26] = estimate.skewunit(x, family1 = "asin", family2 = NULL, 
            est.var = FALSE)$AIC
        fit[27] = estimate.skewunit(x, family1 = "Uquad", family2 = NULL, 
            est.var = FALSE)$AIC
        fit[28] = estimate.skewunit(x, family1 = "triang", family2 = NULL, 
            est.var = FALSE)$AIC
        fit[29] = estimate.skewunit(x, family1 = "JSB", family2 = NULL, 
            est.var = FALSE)$AIC
        fit[30] = estimate.skewunit(x, family1 = "sbeta", family2 = NULL, 
            est.var = FALSE)$AIC
    }
    if (criteria == "BIC") {
        fit[1] = estimate.skewunit(x, family1 = "asin", family2 = "asin", 
            est.var = FALSE)$BIC
        fit[2] = estimate.skewunit(x, family1 = "asin", family2 = "Uquad", 
            est.var = FALSE)$BIC
        fit[3] = estimate.skewunit(x, family1 = "asin", family2 = "triang", 
            est.var = FALSE)$BIC
        fit[4] = estimate.skewunit(x, family1 = "asin", family2 = "JSB", 
            est.var = FALSE)$BIC
        fit[5] = estimate.skewunit(x, family1 = "asin", family2 = "sbeta", 
            est.var = FALSE)$BIC
        fit[6] = estimate.skewunit(x, family1 = "Uquad", family2 = "asin", 
            est.var = FALSE)$BIC
        fit[7] = estimate.skewunit(x, family1 = "Uquad", family2 = "Uquad", 
            est.var = FALSE)$BIC
        fit[8] = estimate.skewunit(x, family1 = "Uquad", family2 = "triang", 
            est.var = FALSE)$BIC
        fit[9] = estimate.skewunit(x, family1 = "Uquad", family2 = "JSB", 
            est.var = FALSE)$BIC
        fit[10] = estimate.skewunit(x, family1 = "Uquad", family2 = "sbeta", 
            est.var = FALSE)$BIC
        fit[11] = estimate.skewunit(x, family1 = "triang", family2 = "asin", 
            est.var = FALSE)$BIC
        fit[12] = estimate.skewunit(x, family1 = "triang", family2 = "Uquad", 
            est.var = FALSE)$BIC
        fit[13] = estimate.skewunit(x, family1 = "triang", family2 = "triang", 
            est.var = FALSE)$BIC
        fit[14] = estimate.skewunit(x, family1 = "triang", family2 = "JSB", 
            est.var = FALSE)$BIC
        fit[15] = estimate.skewunit(x, family1 = "triang", family2 = "sbeta", 
            est.var = FALSE)$BIC
        fit[16] = estimate.skewunit(x, family1 = "JSB", family2 = "asin", 
            est.var = FALSE)$BIC
        fit[17] = estimate.skewunit(x, family1 = "JSB", family2 = "Uquad", 
            est.var = FALSE)$BIC
        fit[18] = estimate.skewunit(x, family1 = "JSB", family2 = "triang", 
            est.var = FALSE)$BIC
        fit[19] = estimate.skewunit(x, family1 = "JSB", family2 = "JSB", 
            est.var = FALSE)$BIC
        fit[20] = estimate.skewunit(x, family1 = "JSB", family2 = "sbeta", 
            est.var = FALSE)$BIC
        fit[21] = estimate.skewunit(x, family1 = "sbeta", family2 = "asin", 
            est.var = FALSE)$BIC
        fit[22] = estimate.skewunit(x, family1 = "sbeta", family2 = "Uquad", 
            est.var = FALSE)$BIC
        fit[23] = estimate.skewunit(x, family1 = "sbeta", family2 = "triang", 
            est.var = FALSE)$BIC
        fit[24] = estimate.skewunit(x, family1 = "sbeta", family2 = "JSB", 
            est.var = FALSE)$BIC
        fit[25] = estimate.skewunit(x, family1 = "sbeta", family2 = "sbeta", 
            est.var = FALSE)$BIC
        fit[26] = estimate.skewunit(x, family1 = "asin", family2 = NULL, 
            est.var = FALSE)$BIC
        fit[27] = estimate.skewunit(x, family1 = "Uquad", family2 = NULL, 
            est.var = FALSE)$BIC
        fit[28] = estimate.skewunit(x, family1 = "triang", family2 = NULL, 
            est.var = FALSE)$BIC
        fit[29] = estimate.skewunit(x, family1 = "JSB", family2 = NULL, 
            est.var = FALSE)$BIC
        fit[30] = estimate.skewunit(x, family1 = "sbeta", family2 = NULL, 
            est.var = FALSE)$BIC
    }
    sfit <- sort(fit)
	ind=order(fit)
    dist1 <- c(rep(c("asin", "Uquad", "triang", "JSB", "sbeta"), 
        c(5, 5, 5, 5, 5)),c("asin", "Uquad", "triang", "JSB", "sbeta"))
    dist2 <- c(rep(c("asin", "Uquad", "triang", "JSB", "sbeta"), 5), rep("",5))
    param <- ifelse(dist1 == "JSB" | dist1 == "sbeta", 1, 0) + 
        ifelse(dist2 == "JSB" | dist2 == "sbeta", 1, 0) + 1-ifelse(dist2=="",1,0)
    resu <- data.frame(dist1, dist2, round(fit, 2), param)[ind, 
        ]
    rownames(resu) <- 1:nrow(resu)
    colnames(resu) <- c("f", "G", criteria, "parameters")
    final <- estimate.skewunit(x, family1 = resu[1, 1], family2 = resu[1, 
        2])
    object.out = list(summary = resu, coefficients = final$coefficients, 
        logLik = final$logLik, AIC = final$AIC, BIC = final$BIC, 
        family1 = resu[1, 1], family2 = resu[1, 2])
    class(object.out) <- "skewunit"
    object.out
}
