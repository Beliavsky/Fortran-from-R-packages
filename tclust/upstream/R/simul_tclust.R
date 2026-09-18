#' Simulate contaminated data set for applying TCLUST
#' 
#' @description Simulate 10\% contaminated data set for applying TCLUST
#' 
#' @param n number of observations
#' @param p dimension (p>=2 and p>q)
#' @param k number of clusters (only k=3 and k=6 are allowed!!!)
#' @param type 1 (spherical for rest.fact=1) or 2 (elliptical for rest.fact=9^2)
#' @param balanced 1 (all clusters equal size) or 2 [proportions (25,30,35)\% if k=3 and (12.5,15,17.5,12.5,15,17.5)\% if k=6]
#' 
#' @return a list with the following items
#' \itemize{
#'     \item x - The generated dataset
#'     \item true - The true classification
#' }
#' 
#' @export
#' 
#' @examples 
#' res <- simula.tclust(n=400,k=3,p=8,type=2,balanced=1)
#' plot(res$x,col=res$true+1)
#'
simula.tclust <- function(n,p=4,k=3,type=2,balanced=1){
  
  if (k!=3 & k!=6) stop("Invalid k")
  if (p<2) stop("Invalid p")
  
  if(k==6){nn <- n/2}else{nn <- n}
  
  if (balanced==1){ nn1=nn2=nn3=ceiling(nn*0.3)}else{
              nn1=ceiling(nn*0.25)
              nn2=ceiling(nn*0.3)
              nn3=ceiling(nn*0.35)}
  
  nh <- nn1+nn2+nn3

  Y <- matrix(NA,ncol=p,nrow=n)
  
  if (type==1){
    ### type=1
    if(k==3){
        Y[1:nn1,] <- MASS::mvrnorm(nn1,c(4,4),matrix(c(1,0,0,1),nrow=2))
        Y[(nn1+1):(nn1+nn2),] <- MASS::mvrnorm(nn2,c(-4,4),matrix(c(1,0,0,1),nrow=2))
        Y[(nn1+nn2+1):(nn1+nn2+nn3),] <- MASS::mvrnorm(nn3,c(0,0),matrix(c(1,0,0,1),nrow=2))
        Y[(nh+1):n,] <- matrix(runif(p*(n-nh),-20,20),ncol=p)
        if(p>2){Y[1:n,3:p] <- MASS::mvrnorm(n,rep(0,p-2),diag(p-2))}
      }
      else
      {
      Y[1:nn1,] <- MASS::mvrnorm(nn1,c(4,4),matrix(c(1,0,0,1),nrow=2))
      Y[(nn1+1):(nn1+nn2),] <- MASS::mvrnorm(nn2,c(-4,4),matrix(c(1,0,0,1),nrow=2))
      Y[(nn1+nn2+1):(nn1+nn2+nn3),] <- MASS::mvrnorm(nn3,c(0,0),matrix(c(1,0,0,1),nrow=2))
      Y[(nh+1):(nh+nn1),] <- MASS::mvrnorm(nn1,c(24,0),matrix(c(1,0,0,1),nrow=2))
      Y[(nh+nn1+1):(nh+nn1+nn2),] <- MASS::mvrnorm(nn2,c(16,0),matrix(c(1,0,0,1),nrow=2))
      Y[(nh+nn1+nn2+1):(nh+nn1+nn2+nn3),] <- MASS::mvrnorm(nn3,c(20,4),matrix(c(1,0,0,1),nrow=2))
      Y[(2*nh+1):n,] <- matrix(runif((n-2*nh)*p,-20,40),ncol=p)
      if(p>2){Y[1:n,3:p] <- MASS::mvrnorm(n,rep(0,p-2),diag(p-2))}  
    }
  }
  
  if (type==2){
    ### type=2
    alp1 <- runif(1,0,2*pi)
    u1 <- matrix(c(cos(alp1),sin(alp1),-sin(alp1),cos(alp1)),ncol=2)
    alp2 <- runif(1,0,2*pi)
    u2 <- matrix(c(cos(alp2),sin(alp2),-sin(alp2),cos(alp2)),ncol=2)
  
    if(k==3){
        Y[1:nn1,] <- MASS::mvrnorm(nn1,c(20,20),t(u1)%*%diag(c(1,9^2))%*%u1)
        Y[(nn1+1):(nn1+nn2),] <- MASS::mvrnorm(nn2,c(-20,-20),t(u1)%*%diag(c(9^2,1))%*%u1)
        Y[(nn1+nn2+1):(nn1+nn2+nn3),] <- MASS::mvrnorm(nn3,c(0,0),t(u1)%*%diag(c(3^2,3^2))%*%u1)
        Y[(nn1+nn2+nn3+1):n,] <- matrix(runif(p*(n-nh),-50,50),ncol=p)
        if(p>2){Y[1:n,3:p] <- MASS::mvrnorm(n,rep(0,p-2),diag(p-2))}
      }
      else
      {
        Y[1:nn1,] <- MASS::mvrnorm(nn1,c(20,20),t(u1)%*%diag(c(1,9^2))%*%u1)
        Y[(nn1+1):(nn1+nn2),] <- MASS::mvrnorm(nn2,c(-20,-20),t(u1)%*%diag(c(9^2,1))%*%u1)
        Y[(nn1+nn2+1):(nn1+nn2+nn3),] <- MASS::mvrnorm(nn3,c(0,0),t(u1)%*%diag(c(3^2,3^2))%*%u1)
        Y[(nh+1):(nh+nn1),] <-  MASS::mvrnorm(nn1,c(40,20),t(u2)%*%diag(c(1,9^2))%*%u2)
        Y[(nh+nn1+1):(nh+nn1+nn2),] <- MASS::mvrnorm(nn2,c(0,-20),t(u2)%*%diag(c(9^2,1))%*%u2)
        Y[(nh+nn1+nn2+1):(nh+nn1+nn2+nn3),] <- MASS::mvrnorm(nn3,c(20,0),t(u2)%*%diag(c(3^2,3^2))%*%u2)
        Y[(2*nh+1):n,] <- matrix(runif((n-2*nh)*p,-60,80),ncol=p)
        if(p>2){Y[1:n,3:p] <- MASS::mvrnorm(n,rep(0,p-2),diag(p-2))}  
      }
    }
  
  # true
  true <- rep(0,n)
  if(k==3){true[1:nh] <- rep(1:3,c(nn1,nn2,nn3))}else{true[1:(2*nh)] <- rep(1:6,c(nn1,nn2,nn3,nn1,nn2,nn3))}

  # Return
  return(list(x=Y,true=true))
}
