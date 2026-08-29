




### rgrf and suchlike
rgrf <- function(  # we can also include dgrf...
 n,          # number of samples
 locs,       # either pass locs or set nx for regular/perturbed grid
 nx, ny=nx,  # at least the first needs to be set.
 xlim=c(0,1), ylim=c(0,1), # size of domain for regular/perturbed grid
 tau=0,      # perturbation if between (0,1/2). 
             # Separation of locations is at least 1-2tau 
 
 Covariance, # similar as with mle
 theta,      # (or theta0)

 beta=0,     # mean if length one, regression coefficients if longer
 X,          # 'design' matrix for mean Xbeta

 method=c('chol'),  # based on Choleski factorization
 method.args=list(sparse=FALSE), # list of arguments that can be passed to the corresponding approach
                     # for chol it can be RStruct, chol.args, cov.args
 eps = getOption("spam.eps"),
 drop = TRUE,
 attributes=TRUE,  # should these been passed back?
 ...)  {

    call <- match.call()
    # construct locations:
    if (missing(locs)) {
        dx <- diff(xlim)/nx
        dy <- diff(ylim)/ny
        locs <- expand.grid(x=seq(from=xlim[1]+dx/2, to=xlim[2]-dx/2, length=nx),
                            y=seq(from=ylim[1]+dx/2, to=ylim[2]-dx/2, length=ny))
        
        Nlocs <- nx * ny
        stopifnot(tau >= 0,  tau < 1/2)
         etau <- tau*c(dx,dy)
       
        if ( tau>eps)  locs <- locs + cbind(runif(Nlocs, -etau,etau), runif(Nlocs, -etau,etau))
         
    } else {
        
        locs <- as.matrix(locs)
        Nlocs <- dim(locs)[1]
        etau <- c(0,0)
    }
    
    # construct mean:
    if (length(beta)>1) {
        X <- as.matrix(X)
        stopifnot( dim(X)[1] == Nlocs,  dim(X)[2] == length(beta))
        mu <- c( X %*% beta ) 
    } else if (length(beta)==1) {
       mu <- if(missing(X)) beta else X * beta
    }

    
    if (Covariance == 'cov.nug') {  # slight exception here...
        out <- sweep( as.matrix( ( array(rnorm(n*Nlocs),c(n,Nlocs)) * theta[1] )), 2, mu, "+")
    } else {
        compact <- if (method.args$sparse) 
            (Covariance %in% c("cov.sph","cov.wend1","cov.wend2","cov.wu1","cov.wu3")) else FALSE
    
    
        if (method=='chol') {
            if (compact) {

                # construct distance and covariance matrix:
                distmat <- nearest.dist(locs, upper=NULL, delta=theta[1])
                Sigma <- do.call(Covariance, c(list(distmat, theta), method.args$cov.args))
               
                
                # now similar as in rmvnorm.spam
                if (is(method.args$Rstruct,"spam.chol.NgPeyton"))
                    cholS <- do.call("update.spam.chol.NgPeyton", list( method.args$Rstruct, Sigma, method.args$chol.args))
                else
                    cholS <- do.call("chol.spam", list( Sigma, method.args$chol.args))
                                        # cholS is the upper triangular part of the permutated matrix Sigma
                iord <- ordering(cholS, inv=TRUE)
                
                R <- as.spam(cholS)
                retval <- as.matrix( ( array(rnorm(n*Nlocs),c(n,Nlocs)) %*% R)[,iord,drop=F])
                                        # It is often better to order the sample than the matrix
                                        # R itself.
                out <- sweep(retval, 2, mu, "+")
            } else {


                if (requireNamespace("pkg", quietly = TRUE)) {
                    distmat <- fields::rdist(locs)
                } else {
                    distmat <- as.matrix(dist(locs))
                }

                Sigma <- do.call(Covariance, c(list(distmat, theta), method.args$cov.args))

                R <-   if (is.null(method.args$chol.args)) chol(Sigma) else
                              do.call("chol", list( Sigma, method.args$chol.args)) 
                N <- dim( Sigma)[1]
                
                out <- sweep( as.matrix( ( array(rnorm(n*Nlocs),c(n,Nlocs)) %*% R)), 2, mu, "+")
            }
            
        } else {
            stop("'method' not implemented yet.")
        }
    }
    if(attributes) {
        attr(out, "locs") <- locs
        attr(out, "mean") <- mu
        attr(out, "call") <- call
        attr(out, "effective.tau") <- etau
    }
    if(drop)   return(drop(out))
    return(out)
    
}
