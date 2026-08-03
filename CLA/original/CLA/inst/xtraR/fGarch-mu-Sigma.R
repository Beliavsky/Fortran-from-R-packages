### ---> ../../R/

if(FALSE) { ## ==> in another file -- TODO

    data(SP500, package="FRAPO")
    system.time(muS.sp500 <- muSigma(SP500))   #   26 sec. (lynne, 2017)

    data(NASDAQ, package="FRAPO")
    system.time(muS.nasdaq <- muSigma(NASDAQ)) #  122.3 sec. (lynne, 2017)

    ## Prove they are the same as Yanhao's version in her thesis:
    if(exists("trans.sp")) {# from her (*rds file)
        cat("comparing Yanhao Shi's  'trans.sp' with our muSigma() result:\n")
        stopifnot(all.equal(as.vector(muS.sp500$mu),
                            as.vector(trans.sp$mu), tol = 1e-15))
        stopifnot(all.equal(muS.sp500$covmat,
                            trans.sp $covmat, tol = 1e-15))
    }
    if(exists("trans.nasdaq")) {# from her (*rds file)
        cat("comparing Yanhao Shi's  'trans.nasdaq' with our muSigma() result:\n")
        stopifnot(all.equal(as.vector(  muS.nasdaq$mu),
                            as.vector(trans.nasdaq$mu), tol = 1e-15))
        stopifnot(all.equal(  muS.nasdaq $covmat,
                            trans.nasdaq $covmat, tol = 1e-15))
    }

    object.size(muS.sp500 ) ##  1843784 bytes :  1.8 MB ; ok for CRAN pkg
    object.size(muS.nasdaq) ## 38720584 bytes : 39 Mega too large for CRAN package
    ##

    if(FALSE) { ## as MM: do once
        save   (muS.sp500, file = "~/R/Pkgs/CLA/data/muSig_sp500.rda", compress = "xz")
        saveRDS(muS.sp500, file = "~/R/Pkgs/CLA/data/muSig_sp500.rds", compress = "xz")
    }
}
