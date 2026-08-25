# legacy vgQ functions 

###########################################
###########################################
###
###               OBLIMIN
###
###########################################
###########################################

vgQOLD.oblimin <- function(L, gam=0){
  X <- L^2 %*% (!diag(TRUE,ncol(L))) 
  if (0 != gam) {
     p <- nrow(L)
     X <- (diag(1,p) - matrix(gam/p,p,p)) %*% X
     }
  list(Gq=L*X,
       f=sum(L^2 * X)/4,
       Method=  if (gam == 0)  "Oblimin Quartimin" else
				if (gam == .5) "Oblimin Biquartimin" else
		         paste("Oblimin g=", gam,sep="")  )
}

# original
# vgQOLD.oblimin <- function(L, gam=0){
#   Method <- paste("Oblimin g=",gam,sep="")
#   if (gam == 0) Method <- "Oblimin Quartimin"
#   if (gam == .5) Method <- "Oblimin Biquartimin"
#   if (gam == 1) Method <- "Oblimin Covarimin"
#   k <- ncol(L)
#   p <- nrow(L)
#   N <- matrix(1,k,k)-diag(k)
#   f <- sum(L^2 * (diag(p)-gam*matrix(1/p,p,p)) %*% L^2 %*% N)/4
#   Gq <- L * ((diag(p)-gam*matrix(1/p,p,p)) %*% L^2 %*% N)
#   return(list(Gq=Gq,f=f,Method=Method))
# }

###########################################
###########################################
###
###              QUARTIMIN
###
###########################################
###########################################

vgQOLD.quartimin <- function(L){
  X <-  L^2 %*% (!diag(TRUE,ncol(L))) 
  list(Gq= L*X,
       f= sum(L^2 * X)/4,  
       Method=  "Quartimin" )
  }

#original
#vgQOLD.quartimin <- function(L){
#  Method="Quartimin"
#  L2 <- L^2
#  k <- ncol(L)
#  M <- matrix(1,k,k)-diag(k)
#  f <- sum(L2 * (L2 %*% M))/4
#  Gq <- L * (L2 %*% M)
#  return(list(Gq=Gq,f=f,Method=Method))
#} 


###########################################
###########################################
###
###           TARGET ROTATION
###      required argument is a 
###       target matrix 'Target'
###
###########################################
###########################################


vgQOLD.target <- function(L, Target=NULL){
   if(is.null(Target)) stop("argument Target must be specified.")
   #   e.g.  Target <- matrix(c(rep(NA,4),rep(0,8),rep(NA,4)),8) 
   #  approximates  Michael Brown approach
   Gq <-  2 * (L - Target)
   Gq[is.na(Gq)] <- 0 #missing elements in target do not affect the first derivative 
   list(Gq=Gq,
        f=sum((L-Target)^2, na.rm=TRUE),  
	Method="Target rotation")  #The target rotation ? option in Michael Browne's algorithm should be NA
   }

###########################################
###########################################
###
###   PARTIALLY SPECIFIED TARGET ROTATION
###      required arguments are a 
###      target matrix 'Target' and 
###          weight matrix 'W'
###
###########################################
###########################################


vgQOLD.pst <- function(L, W=NULL, Target=NULL){
   if(is.null(W))      stop("argument W must be specified.")
   if(is.null(Target)) stop("argument Target must be specified.")
   # Needs weight matrix W with 1's at specified values, 0 otherwise
   # e.g. W = matrix(c(rep(1,4),rep(0,8),rep(1,4)),8). 
   # When W has only 1's this is procrustes rotation
   # Needs a Target matrix Target with hypothesized factor loadings.
   # e.g. Target = matrix(0,8,2)
   Btilde <- W * Target
   list(Gq= 2*(W*L-Btilde), 
        f = sum((W*L-Btilde)^2),
        Method="Partially specified target")
}

###########################################
###########################################
###
###               OBLIMAX
###
###########################################
###########################################

vgQOLD.oblimax <- function(L){
  list(Gq= -(4*L^3/(sum(L^4))-4*L/(sum(L^2))),
       f= -(log(sum(L^4))-2*log(sum(L^2))),
       Method="Oblimax")
}


###########################################
###########################################
###
###          MINIMUM ENTROPY
###
###########################################
###########################################


vgQOLD.entropy <- function(L){
  list(Gq= -(L*log(L^2 + (L^2==0)) + L),
       f= -sum(L^2*log(L^2 + (L^2==0)))/2, 
       Method="Minimum entropy")
}

###########################################
###########################################
###
###               QUARTIMAX
###
###########################################
###########################################


vgQOLD.quartimax <- function(L){
  list(Gq= -L^3,
       f= -sum(diag(crossprod(L^2)))/4, 
       Method="Quartimax")
}


###########################################
###########################################
###
###               VARIMAX
###
###########################################
###########################################


vgQOLD.varimax <- function(L){
  QL <- sweep(L^2,2,colMeans(L^2),"-")
  list(Gq= -L * QL,
       f= -sqrt(sum(diag(crossprod(QL))))^2/4, 
       Method="varimax")
}

###########################################
###########################################
###
###               SIMPLIMAX
###  argument: # close to zero loadings 'k'
###
###########################################
###########################################


vgQOLD.simplimax <- function(L, k=nrow(L)){
  # k: Number of close to zero loadings
  Imat <- sign(L^2 <= sort(L^2)[k])
  list(Gq= 2*Imat*L,
       f= sum(Imat*L^2), 
       Method="Simplimax")
}

###########################################
###########################################
###
### BENTLER'S INVARIANT PATTERN SIMPLICITY
###
###########################################
###########################################

vgQOLD.bentler <- function(L){
  L2 <- L^2
  M <- crossprod(L2)
  D <- diag(diag(M))
  list(Gq= -L * (L2 %*% (solve(M)-solve(D))),
       f= -(log(det(M))-log(det(D)))/4,
       Method="Bentler's criterion")
}

###########################################
###########################################
###
###           TANDEM CRITERIA
###
###########################################
###########################################


#vgQOLD.tandemI <- function(L){  # Tandem Criterion, Comrey, 1967.
#  Method <- "Tandem I"
#  LL <- (L %*% t(L))
#  LL2 <- LL^2
#  f <- -sum(diag(crossprod(L^2, LL2 %*% L^2)))
#  Gq1 <- 4 * L *(LL2 %*% L^2)
#  Gq2 <- 4 * (LL * (L^2 %*% t(L^2))) %*% L
#  Gq <- -Gq1 - Gq2 
#  return(list(Gq=Gq,f=f,Method=Method))
#}

vgQOLD.tandemI <- function(L){  # Tandem Criterion, Comrey, 1967.
  LL <- (L %*% t(L))
  LL2 <- LL^2
  Gq1 <- 4 * L *(LL2 %*% L^2)
  Gq2 <- 4 * (LL * (L^2 %*% t(L^2))) %*% L
  Gq <- -Gq1 - Gq2 
  list(Gq=Gq,
       f= -sum(diag(crossprod(L^2, LL2 %*% L^2))), 
       Method="Tandem I")
  }

vgQOLD.tandemII <- function(L){  # Tandem Criterion, Comrey, 1967.
  LL <- (L %*% t(L))
  LL2 <- LL^2
  f <- sum(diag(crossprod(L^2, (1-LL2) %*% L^2)))
  Gq1 <- 4 * L *((1-LL2) %*% L^2)
  Gq2 <- 4 * (LL * (L^2 %*% t(L^2))) %*% L
  Gq <- Gq1 - Gq2 
  list(Gq=Gq,
       f=f, 
       Method="Tandem II")
  }

###########################################
###########################################
###
###               GEOMIN
###
###########################################
###########################################


vgQOLD.geomin <- function(L, delta=.01){
  k <- ncol(L)
  p <- nrow(L)
  L2 <- L^2 + delta
  pro <- exp(rowSums(log(L2))/k) 
  list(Gq=(2/k)*(L/L2)*matrix(rep(pro,k),p),
       f= sum(pro), 
       Method="Geomin")
  }
  
###########################################
###########################################
###
###               BI-GEOMIN
###
###########################################
###########################################

  
vgQOLD.bigeomin <- function(L, delta = 0.01){
   Lg <- L[ , -1, drop = FALSE]
   out <- vgQOLD.geomin(Lg, delta = delta) 
     list(Gq=cbind(0, out$Gq),
       f= out$f, 
       Method="Bi-Geomin")
  }

###########################################
###########################################
###
###          CRAWFORD FERGUSON FAMILY
###           needs kappa parameter
###
###             EQUAMAX   PARSIMAX
###
###########################################
###########################################


vgQOLD.cf <- function(L, kappa=0){
  k <- ncol(L)
  p <- nrow(L)
  # kappa <- 0 # quartimax 
  # kappa <- 1/p # varimax
  # kappa <- k/(2*p) # equamax
  # kappa <- (k-1)/(p+k-2) # parsimax
  # kappa <- 1 # factor parsimony
  N <- matrix(1,k,k)-diag(k)
  M <- matrix(1,p,p)-diag(p)
  L2 <- L^2
  f1 <- (1-kappa)*sum(diag(crossprod(L2,L2 %*% N)))/4
  f2 <- kappa*sum(diag(crossprod(L2,M %*% L2)))/4
  list(Gq= (1-kappa) * L * (L2 %*% N) + kappa * L * (M %*% L2),
       f= f1 + f2,
       Method=  if (kappa == 0)  "Crawford-Ferguson Quartimax/Quartimin" else
				if (kappa == 1/p) "Crawford-Ferguson Varimax" else
				if (kappa == k/(2*p)) "Equamax" else
				if (kappa == (k-1)/(p+k-2)) "Parsimax" else
				if (kappa == 1)  "Factor Parsimony"   else
				paste("Crawford-Ferguson:k=",kappa,sep=""))
}

###########################################
###########################################
###
###               INFOMAX
###
###########################################
###########################################



vgQOLD.infomax <- function(L){
  k <- ncol(L)
  p <- nrow(L)
  S <- L^2
  s <- sum(S)
  s1 <- rowSums(S)
  s2 <- colSums(S)
  E <- S/s
  e1 <- s1/s
  e2 <- s2/s
  Q0 <- sum(-E * log(E))
  Q1 <- sum(-e1 * log(e1))
  Q2 <- sum(-e2 * log(e2))
  f <- log(k) + Q0 - Q1 - Q2
  H <- -(log(E) + 1)
  alpha <- sum(S * H)/s^2
  G0 <- H/s - alpha * matrix(1, p, k)
  h1 <- -(log(e1) + 1)
  alpha1 <- s1 %*% h1/s^2
  G1 <- matrix(rep(h1,k), p)/s - as.vector(alpha1) * matrix(1, p, k)
  h2 <- -(log(e2) + 1)
  alpha2 <- h2 %*% s2/s^2
  G2 <- matrix(rep(h2,p), ncol=k, byrow=T)/s - as.vector(alpha2) * matrix(1, p, k)
  Gq <- 2 * L * (G0 - G1 - G2)
  list(f = f,
  		Gq=Gq,
  		Method="Infomax")
}

###########################################
###########################################
###
###           MCCAMMON ENTROPY
###
###########################################
###########################################


vgQOLD.mccammon <- function(L){
  k <- ncol(L)
  p <- nrow(L)
  S <- L^2
  M <- matrix(1,p,p)
  s2 <- colSums(S)
  P <- S / matrix(rep(s2,p),ncol=k,byrow=T)
  Q1 <- -sum(P * log(P))
  H <- -(log(P) + 1)
  R <- M %*% S
  G1 <- H/R - M %*% (S*H/R^2)
  s <- sum(S)
  p2 <- s2/s
  Q2 <- -sum(p2 * log(p2))
  h <- -(log(p2) + 1)
  alpha <- h %*% p2
  G2 <- rep(1,p) %*% t(h)/s - as.vector(alpha)*matrix(1,p,k)
  Gq <- 2*L*(G1/Q1 - G2/Q2)
  f <- log(Q1) - log(Q2)
  list(f = f,
		Gq = Gq, 
		Method = "McCammon entropy")
}

###########################################
###########################################
###
###           BIFACTOR BIQUARTIMIN
###
###########################################
###########################################


vgQOLD.bifactor <- function(L){
	k <- ncol(L)
	Lt <- L[,2:k]
	Lt2 <- Lt^2
	N <- matrix(1, nrow=k-1, ncol=k-1) - diag(k-1)
	f <- sum(Lt2 * (Lt2 %*% N))
	Gt <- 4 * Lt * (Lt2 %*% N)
	G <- cbind(0, Gt)
	list(f = f,
			Gq = G,
			Method = "Bifactor Biquartimin") 
}

###########################################
###########################################
###
###               VARIMIN
###
###########################################
###########################################


   
vgQOLD.varimin <- function (L){
    QL <- sweep(L^2, 2, colMeans(L^2), "-")
    list(Gq = L * QL, 
        f = sqrt(sum(diag(crossprod(QL))))^2/4,
        Method = "varimin")
}

###########################################
###########################################
###
###             Lp rotation
###
###########################################
###########################################

# for lpT and lpQ functions see the file lp.R

vgQOLD.lp.wls<-function(L, W) {
   list(
     Gq = 2 * W * L/nrow(L),                 
     f = sum(W * L * L)/nrow(L),             
     Method = "Weighted least squares for Lp rotation" # Method description
   )
 }

###########################################
###########################################
###
###               PROMAX 
###            (not in use)
###
###########################################
###########################################

# promax is already defined in the stats (previously mva) package
# 
#GPromax <- function(A,pow=3){
# method <- "Promax"
# # Initial rotation: Standardized Varimax
# require(statsa)
# xx <- promax(A,pow)
# Lh <- xx$loadings
# Th <- xx$rotmat
# orthogonal <- F
# Table <- NULL
#return(list(loadings=Lh,Th=Th,Table=NULL,method,orthogonal=orthogonal))
#}
