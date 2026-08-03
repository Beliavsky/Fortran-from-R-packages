\name{bus}
\alias{bus}
\docType{data}
\title{ Automatic vehicle recognition data}
\description{
    The data set bus (Hettich and Bay, 1999) corresponds to a study in 
    automatic vehicle recognition (see Maronna et al. 2006, page 213, Example 6.3)). 
    This data set from the 
    Turing Institute, Glasgow, Scotland, contains measures of shape 
    features extracted from vehicle silhouettes. The images were acquired 
    by a camera looking downward at the model vehicle from a fixed angle 
    of elevation. Each of the 218 rows corresponds to a view of a bus 
    silhouette, and contains 18 attributes of the image. 
}
\usage{data(bus)}
\format{
    A data frame with 218 observations on the following 18 variables: 

    \describe{
        \item{\code{V1}}{compactness}
        \item{\code{V2}}{circularity}
        \item{\code{V3}}{distance circularity}
        \item{\code{V4}}{radius ratio}
        \item{\code{V5}}{principal axis aspect ratio}
        \item{\code{V6}}{maximum length aspect ratio}
        \item{\code{V7}}{scatter ratio}
        \item{\code{V8}}{elongatedness}
        \item{\code{V9}}{principal axis rectangularity}
        \item{\code{V10}}{maximum length rectangularity}
        \item{\code{V11}}{scaled variance along major axis}
        \item{\code{V12}}{scaled variance along minor axis}
        \item{\code{V13}}{scaled radius of gyration}
        \item{\code{V14}}{skewness about major axis}
        \item{\code{V15}}{skewness about minor axis}
        \item{\code{V16}}{kurtosis about minor axis}
        \item{\code{V17}}{kurtosis about major axis}
        \item{\code{V18}}{hollows ratio}
    }
}
%%\details{}
\source{
    Hettich, S. and Bay, S.D. (1999), The UCI KDD Archive,
    Irvine, CA:University of California, Department of Information and Computer Science,
%%    \url{http://kdd.ics.uci.edu}
    'http://kdd.ics.uci.edu'
}
\references{
    Maronna, R., Martin, D. and Yohai, V., (2006). Robust Statistics: Theory and Methods. Wiley, New York.
}
\examples{
    ## Reproduce Table 6.3 from Maronna et al. (2006), page 213
    data(bus)
    bus <- as.matrix(bus)
    
    ## calculate MADN for each variable
    xmad <- apply(bus, 2, mad)      
    cat("\nMin, Max of MADN: ", min(xmad), max(xmad), "\n")


    ## MADN vary between 0 (for variable 9) and 34. Therefore exclude 
    ##  variable 9 and divide the remaining variables by their MADNs.
    bus1 <- bus[, -9]
    madbus <- apply(bus1, 2, mad)
    bus2 <- sweep(bus1, 2, madbus, "/", check.margin = FALSE)

    ## Compute classical and robust PCA (Spherical/Locantore, Hubert, MCD and OGK)    
    pca  <- PcaClassic(bus2)
    rpca <- PcaLocantore(bus2)
    pcaHubert <- PcaHubert(bus2, k=17, kmax=17, mcd=FALSE)
    pcamcd <- PcaCov(bus2, cov.control=CovControlMcd())
    pcaogk <- PcaCov(bus2, cov.control=CovControlOgk())

    ev    <- getEigenvalues(pca)
    evrob <- getEigenvalues(rpca)
    evhub <- getEigenvalues(pcaHubert)
    evmcd <- getEigenvalues(pcamcd)
    evogk <- getEigenvalues(pcaogk)

    uvar    <- matrix(nrow=6, ncol=6)
    svar    <- sum(ev)
    svarrob <- sum(evrob)
    svarhub <- sum(evhub)
    svarmcd <- sum(evmcd)
    svarogk <- sum(evogk)
    for(i in 1:6){
        uvar[i,1] <- i
        uvar[i,2] <- round((svar - sum(ev[1:i]))/svar, 3)
        uvar[i,3] <- round((svarrob - sum(evrob[1:i]))/svarrob, 3)
        uvar[i,4] <- round((svarhub - sum(evhub[1:i]))/svarhub, 3)
        uvar[i,5] <- round((svarmcd - sum(evmcd[1:i]))/svarmcd, 3)
        uvar[i,6] <- round((svarogk - sum(evogk[1:i]))/svarogk, 3)
    }
    uvar <- as.data.frame(uvar)
    names(uvar) <- c("q", "Classical","Spherical", "Hubert", "MCD", "OGK")
    cat("\nBus data: proportion of unexplained variability for q components\n")
    print(uvar)
 
    ## Reproduce Table 6.4 from Maronna et al. (2006), page 214
    ##
    ## Compute classical and robust PCA extracting only the first 3 components
    ## and take the squared orthogonal distances to the 3-dimensional hyperplane
    ##
    pca3 <- PcaClassic(bus2, k=3)               # classical
    rpca3 <- PcaLocantore(bus2, k=3)            # spherical (Locantore, 1999)
    hpca3 <- PcaHubert(bus2, k=3)               # Hubert
    dist <- pca3@od^2
    rdist <- rpca3@od^2
    hdist <- hpca3@od^2

    ## calculate the quantiles of the distances to the 3-dimensional hyperplane
    qclass  <- round(quantile(dist, probs = seq(0, 1, 0.1)[-c(1,11)]), 1)
    qspc <- round(quantile(rdist, probs = seq(0, 1, 0.1)[-c(1,11)]), 1)
    qhubert <- round(quantile(hdist, probs = seq(0, 1, 0.1)[-c(1,11)]), 1)
    qq <- cbind(rbind(qclass, qspc, qhubert), round(c(max(dist), max(rdist), max(hdist)), 0))
    colnames(qq)[10] <- "Max"
    rownames(qq) <- c("Classical", "Spherical", "Hubert")
    cat("\nBus data: quantiles of distances to hiperplane\n")
    print(qq)

    ## 
    ## Reproduce Fig 6.1 from Maronna et al. (2006), page 214
    ## 
    cat("\nBus data: Q-Q plot of logs of distances to hyperplane (k=3) 
    \nfrom classical and robust estimates. The line is the identity diagonal\n")
    plot(sort(log(dist)), sort(log(rdist)), xlab="classical", ylab="robust")
    lines(sort(log(dist)), sort(log(dist)))
   
    
}
\keyword{datasets}
