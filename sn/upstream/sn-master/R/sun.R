#  file sn/R/sun.R   
#  This file is a component of the R package 'sn' 
#  copyright (C) 1997-2021 Adelchi Azzalini
# 
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 2 or 3 of the License
#  (at your option).
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/
#-------------------------------------------
#
# Some support functions
#
all.numeric <- function(...)
{# check if all elements are numeric
  lst <- list(...)
  n <- length(lst)
  if(n == 0) return(NULL)
  m <- is.numeric(lst[[1]])
  if(n == 1) return(m)
  for(k in 2:n) m <- m & is.numeric(lst[[k]])
  return(m)    
} 

blockDiag <- function(...)
{# create a block-diagonal matrix from a set of matrices 
  lst <- list(...)
  n <- length(lst)
  if(n == 0) return(NULL)
  m <- as.matrix(lst[[1]])
  if(n == 1) return(m)
  for(k in 2:n) {
    mk <- as.matrix(lst[[k]])
    m <- rbind(cbind(m, matrix(0, nrow(m), ncol(mk))), 
             cbind(matrix(0, nrow(mk), ncol(m)), mk))
  }           
  return(m)    
} 

tr <- function(x) 
{# trace of a numeric square matrix
  if(mode(x) != "numeric") stop("not a numeric argument")
  if(is.matrix(x)) {
    if(ncol(x) == nrow(x)) sum(diag(x)) else stop("not a square matrix")} 
  else if(length(x)==1) x else stop("not a square matrix") 
}
#-------------------------------------------

dsun <- function(x, xi, Omega, Delta, tau, Gamma, dp=NULL, log=FALSE, silent=FALSE, ...) 
{# SUN density function  
  if (!(missing(Delta) & missing(Omega)) && !is.null(dp)) 
        stop("You cannot set both component parameters and 'dp'")
    if (!is.null(dp)) {
        if (length(dp) != 5)  stop("wrong length of non-null 'dp'")
        xi <- drop(dp[[1]])
        Omega <- dp[[2]]
        Delta <- dp[[3]]
        tau <-  dp[[4]]
        Gamma <- dp[[5]]
    }
  if(!all.numeric(x, xi, Omega, Delta, tau, Gamma)) stop("non-numeric argument(s)")  
  d <- dim(Omega)[1]
  if(length(xi)!= d | dim(Omega)[2] != d)  stop("mismatch of dimensions")
  omega <- sqrt(diag(Omega))
  Omega.bar <- cov2cor(Omega)
  O.inv <- solve(Omega.bar)  
  m <- length(tau)
  if(m==1  & !silent) 
    warning("When m=1, functions for the SN/ESN distr'n are preferable")
  if(any(dim(Gamma) != c(m,m) | dim(Delta) != c(d,m)))
    stop("mismatch of dimensions")
  x <- if(is.vector(x)) matrix(x, 1, d) else data.matrix(x) 
  n <- nrow(x)
  if(m > 20) 
    {if(silent) return(rep(NA, n)) else stop("m exceeds the admissible size")}
  if (is.vector(xi)) xi <- outer(rep(1, n), as.vector(matrix(xi, 1, d)))
  tz <- t(x-xi)/omega
  D.Oinv <- t(Delta) %*% O.inv
  p1 <- pmnorm(t(tau + D.Oinv %*% tz), rep(0,m), Gamma - D.Oinv %*% Delta, ...)
  p2 <- pmnorm(tau, rep(0,m), Gamma, ...)
  if(n == 1) {
    if(any(c(attr(p1,"status"), attr(p2,"status")) != "normal completion"))
      warning("return status from pmnorm is not 'normal completion'")
    }
  pdfN <- dmnorm(x, xi, Omega, log=log)
  if(log) pdfN + logb(p1) - logb(p2)
  else pdfN * p1/p2
}

psun <- function(x, xi, Omega, Delta, tau, Gamma, dp=NULL, log=FALSE, silent=FALSE, ...) 
{# SUN distribution function    
  if (!(missing(Delta) & missing(Omega)) && !is.null(dp)) 
        stop("You cannot set both component parameters and 'dp'")
    if (!is.null(dp)) {
        if (length(dp) != 5)  stop("wrong length of non-null 'dp'")
        xi <- drop(dp[[1]])
        Omega <- dp[[2]]
        Delta <- dp[[3]]
        tau <-  dp[[4]]
        Gamma <- dp[[5]]
    }
  if(!all.numeric(x, xi, Omega, Delta, tau, Gamma)) stop("non-numeric argument(s)")    
  d <- dim(Omega)[1]
  if(dim(Omega)[2] != d)  stop("mismatch of dimensions")
  omega <- sqrt(diag(Omega))
  Omega.bar <- cov2cor(Omega)
  O.inv <- solve(Omega.bar)  
  m <- length(tau)
  if(m==1 & !silent)
    warning("When m=1, functions for the SN/ESN distribution are preferable")
  if(any(dim(Gamma) != c(m,m) | dim(Delta) != c(d,m)))
    stop("mismatch of dimensions")
  x <- if(is.vector(x)) matrix(x, 1, d) else data.matrix(x) 
  n <- nrow(x)
  if((d+m) > 20) 
    {if(silent) return(rep(NA, n)) else stop("(d+m) exceeds the admissible size")}
  if (is.vector(xi)) xi <- outer(rep(1, n), as.vector(matrix(xi, 1, d)))
  if(ncol(x) != ncol(xi)) stop("mismatch of dimensions")
  tz <- t(x-xi)/omega  
  y <- cbind(t(tz), outer(rep(1, n), tau))
  Omega.starNeg <- rbind(cbind(Omega.bar, -Delta),  cbind(t(-Delta), Gamma))
  p1 <- pmnorm(y, mean=rep(0, m+d), varcov=Omega.starNeg, ...)
  p2 <- pmnorm(tau, rep(0,m), Gamma, ...)
  if(n==1) {
    if(any(c(attr(p1,"status"), attr(p2,"status")) != "normal completion"))
    warning("return status from pmnorm is not 'normal completion'")
    }
  as.numeric(pmin(1, pmax(0, p1/p2)))
}

rsun <- function(n=1, xi, Omega, Delta, tau, Gamma, dp=NULL, silent=FALSE) 
{# SUN random numbers, use (7.4) of SN book  
  if (!(missing(Delta) & missing(Omega)) && !is.null(dp)) 
        stop("You cannot set both component parameters and 'dp'")
    if (!is.null(dp)) {
        if (length(dp) != 5)  stop("wrong length of non-null 'dp'")
        xi <- drop(dp[[1]])
        Omega <- dp[[2]]
        Delta <- dp[[3]]
        tau <-  dp[[4]]
        Gamma <- dp[[5]]
    }
  d <- dim(Omega)[1]
  if(length(xi)!= d | dim(Omega)[2] != d)  stop("mismatch of dimensions")
  omega <- sqrt(diag(Omega))
  Omega.bar <- cov2cor(Omega)
  # O.inv <- solve(Omega.bar)  
  m <- length(tau)
  if(m==1 & !silent) 
    warning("When m=1, functions for the SN/ESN family are preferable")
  if(any(dim(Gamma) != c(m,m) | dim(Delta) != c(d,m)))
    stop("mismatch of dimensions")
  Delta_invGamma <- Delta %*% solve(Gamma)  
  Psi.bar <- Omega.bar - Delta_invGamma %*%  t(Delta)
  u0 <- mnormt::rmnorm(n, rep(0, d), Psi.bar)
  u1 <- mnormt::rmtruncnorm(n, rep(0, m), Gamma, -tau)
  tz <- t(u0) + Delta_invGamma %*% t(u1)
  t(xi + omega * tz)
}  

#-------------------------

sunMean <- function(xi, Omega, Delta, tau, Gamma, dp=NULL, silent=FALSE, ...) 
{# expected value of SUN distribution
  if (!(missing(Delta) & missing(Omega)) && !is.null(dp)) 
        stop("You cannot set both component parameters and 'dp'")
    if (!is.null(dp)) {
        if (length(dp) != 5)  stop("wrong length of non-null 'dp'")
        xi <- drop(dp[[1]])
        Omega <- dp[[2]]
        Delta <- dp[[3]]
        tau <-  dp[[4]]
        Gamma <- dp[[5]]
    }
  if(!all.numeric(xi, Omega, Delta, tau, Gamma)) stop("non-numeric argument(s)")    
  d <- dim(Omega)[1]
  if(length(xi)!= d | dim(Omega)[2] != d)  stop("mismatch of dimensions")
  omega <- sqrt(diag(Omega))
  Omega.bar <- cov2cor(Omega)
  O.inv <- solve(Omega.bar)  
  m <- length(tau)
  if(m==1 & !silent) 
    warning("When m=1, functions for the SN/ESN family are preferable")
  if(any(dim(Gamma) != c(m,m) | dim(Delta) != c(d,m)))
    stop("mismatch of dimensions") 
  if(m > 20) 
    {if(silent) return(NA) else stop("m exceeds the admissible size")}   
  prob <- mnormt::pmnorm(tau, rep(0, m), Gamma, ...)
  if(m > 3  &&  (attr(prob,"status") != "normal completion") & !silent)
    warning("return status from pmnorm is not 'normal completion'")
  deriv <- dnorm(tau)/prob
  if(m>1) for(k in 1:m) {
      Gk <- Gamma[-k,-k, drop=FALSE]
      gk <- Gamma[-k, k, drop=FALSE]
      Ec <- as.vector(gk * tau[k])
      Vc <- Gk - gk %*% t(gk)  
      deriv[k] <- deriv[k] * pmnorm(tau[-k], Ec, Vc, ...)
      }
  as.numeric(xi + omega*as.vector(Delta %*% deriv))
}
mean.SUNdistr <- function(x) sunMean(dp=slot(x, "dp"), silent=TRUE)

sunVcov <- function(xi, Omega, Delta, tau, Gamma, dp=NULL, silent=FALSE, ...) 
{# variance (matrix) of SUN  distribution, using Proposition1 of RAV&AA-2020
  if (!(missing(Delta) & missing(Omega)) && !is.null(dp)) 
        stop("You cannot set both component parameters and 'dp'")
    if (!is.null(dp)) {
        if (length(dp) != 5)  stop("wrong length of non-null 'dp'")
        xi <- drop(dp[[1]])
        Omega <- dp[[2]]
        Delta <- dp[[3]]
        tau <-  dp[[4]]
        Gamma <- dp[[5]]
    }
  if(!all.numeric(xi, Omega, Delta, tau, Gamma)) stop("non-numeric argument(s)")  
  d <- dim(Omega)[1]
  if(length(xi)!= d | dim(Omega)[2] != d)  stop("mismatch of dimensions")
  omega <- sqrt(diag(Omega))
  Omega.bar <- cov2cor(Omega)
  O.inv <- solve(Omega.bar)  
  m <- length(tau)
  if(m > 20) 
    {if(silent) return(NA) else stop("m exceeds the admissible size")}  
  if(m==1 & !silent) 
    warning("When m=1, functions for the SN/ESN family are preferable")
  if(any(dim(Gamma) != c(m,m) | dim(Delta) != c(d,m))) 
    stop("mismatch of dimensions") 
  mom.U <- mnormt::mom.mtruncnorm(2, mean=rep(0,m), Gamma, lower=-tau, ...)
  omega.Delta <- omega * Delta
  Gamma.inv <- solve(Gamma)
  A <- omega.Delta %*% Gamma.inv
  B.BT <- Omega - omega.Delta %*% Gamma.inv %*% t(omega.Delta)
  E.U <- if(m==1) mom.U$cum[1] else mom.U$cum1
  E.U2 <- if(m==1) mom.U$mom[3] else mom.U$order2$m2
  var.U <- if(m==1) mom.U$cum[2] else mom.U$order2$cum2
  return(Omega - A %*% (Gamma- var.U) %*% t(A))
}
vcov.SUNdistr <- function(object) sunVcov(dp=slot(object, "dp"), silent=TRUE)
#-------------------------------------------
# expand array to matrix (which are used by RAV&AA-2020)
array2mat <- function(x, d) 
  if(length(x)==d |  length(dim(x))==2) return(x)   else  
  {
    n <- length(dim(x))
    if(n > 4) stop("length(dim(x))>4 not allowed")
    out <- NULL
    for(k in 1:d) { 
      s1 <- if(n==3) paste("x[, , k]") else paste("x[, , k, 1]")
      m1 <- eval(str2expression(s1))
      out <- rbind(out, m1)
     }
    if(n==4)  for(j in 2:d)  out <- cbind(out, array2mat(x[,,,j], d))
    return(out) 
    }

sunMardia <- function(xi, Omega, Delta, tau, Gamma, dp=NULL, silent=FALSE, ...) 
{# Mardia measures of multivariate skewness and kurtosis for SUN distributions 
  if(!(missing(Delta) & missing(Omega)) && !is.null(dp)) 
     stop("You cannot set both component parameters and 'dp'")
  if(!is.null(dp)) {
    if (length(dp) != 5)  stop("wrong length of non-null 'dp'")
    xi <- drop(dp[[1]])
    Omega <- dp[[2]]
    Delta <- dp[[3]]
    tau <-  dp[[4]]
    Gamma <- dp[[5]]
    }
  if(!all.numeric(xi, Omega, Delta, tau, Gamma)) stop("non-numeric argument(s)")
  d <- length(xi)
  m <- length(tau)  
  compNames <- rownames(Omega)
  HcompNames <- rownames(Gamma)
  if(is.null(compNames)) compNames <- paste("V", 1:d, sep="")
  if(is.null(HcompNames)) HcompNames <- paste("H", 1:m, sep="")
  u <- sunValues(dp=dp, compNames, HcompNames, ...)
  return(u$mardia)
}

makeSUNdistr <- function(dp, name, compNames, HcompNames, drop=TRUE) 
{
  if(!is.list(dp)) stop("dp is not a list")
  if(length(dp) != 5) stop("length(dp) is not 5")
  xi <- dp[[1]] 
  Omega <- dp[[2]]
  Delta <- dp[[3]]
  tau <- dp[[4]]
  Gamma <- dp[[5]]
  if(!all.numeric(xi, Omega, Delta, tau, Gamma)) stop("non-numeric argument(s)")
  d <- length(xi)
  m <- length(tau)
  if(!all(dim(Omega) == c(d,d)))  stop("mismatch of dimensions")
  if(missing(compNames)) { compNames <-
    if(length(names(xi)) == d) names(xi) else
      as.vector(outer("V", as.character(1:d), paste,sep="")) }
  if(!is.matrix(Gamma) | m==1) {
    if(length(c(Gamma))>1) stop("Wrong dp$Gamma")
    if(c(Gamma) != 1) stop("Since m=1,  dp$Gamma must be 1, but it is not") 
    if(drop) {
    delta <- c(Delta)
    if(length(delta) != d) stop("wrong size of Delta")
    if(length(tau) != 1) stop("wrong length(tau)")
    Om.delta <- solve(cov2cor(Omega)) %*% delta
    delta.star.sq <- sum(delta %*% Om.delta)
    if(delta.star.sq >= 1 | delta.star.sq < 0) stop("unfeasible arguments")
    alpha <- as.vector(Om.delta)/sqrt(1 - delta.star.sq)
    if(missing(name)) name <- "Unknown_ESN"
    if(d==1) {
      dp.ESN <- c(xi=xi, omega=sqrt(Omega), alpha=alpha, tau=tau)
      obj <- new("SECdistrUv", dp=dp.ESN, family="ESN", name=name)
      } 
    else {
      dp.ESN <- list(xi=xi, Omega=Omega, alpha=alpha, tau=tau)
      obj <- new("SECdistrMv", dp=dp.ESN, family="ESN", name=name, compNames=compNames)
      }
    return(obj)  
    } } 
  if(any(dim(Gamma) != c(m,m)) | any(dim(Delta) != c(d,m)))
    stop("mismatch of dimensions")
  omega <- sqrt(diag(Omega))
  if(!all(diag(Gamma)==1)) stop("diag(Gamma) are not all 1's")
  big.Omega <- rbind(cbind(Omega, omega*Delta), cbind(t(omega*Delta), Gamma))
  if(max(abs(big.Omega -t(big.Omega))) > .Machine$double.eps)
    stop("(Omega, Delta, Gamma) do not make a symmetric matrix")
  big.Omega <- 0.5*(big.Omega + t(big.Omega))  
  eigenvalues <- eigen(big.Omega, symmetric=TRUE, only.values = TRUE)$values 
  if(any(eigenvalues <= 0)) 
    stop("(Omega, Delta, Gamma) do not make a positive definite matrix") 
  name <- if (!missing(name))  as.character(name)[1]
          else paste("Unnamed-SUN(d=", as.character(d), ",m=", 
                     as.character(m), ")", sep = "")
  names(dp) <- c("xi", "Omega", "Delta", "tau", "Gamma")
  if(missing(compNames)) compNames <-
    as.vector(outer("V", as.character(1:d), paste,sep="")) 
  if(missing(HcompNames)) HcompNames <-
    as.vector(outer("H", as.character(1:m), paste,sep="")) 
  names(xi) <- compNames
  dimnames(Omega) <- list(compNames, compNames)    
  dimnames(Delta) <- list(compNames, HcompNames)     
  names(tau) <- HcompNames
  dimnames(Gamma) <- list(HcompNames, HcompNames)      
  dp0 <- list(xi=xi, Omega=Omega, Delta=Delta, tau=tau, Gamma=Gamma)         
  obj <- new("SUNdistr", dp = dp0, name = name, compNames=compNames, 
             HcompNames=HcompNames)
  if(!is(obj, "SUNdistr") & drop==FALSE) 
    stop("Error. No SUNdistr object created")
  obj
}

marginalSUNdistr <- function(object, comp, name, drop=TRUE) 
{# builds from 'obj' the SUN marginal distribution identified by 'comp' 
  # class.obj <- class(object)
  if(!is(object, "SUNdistr") & !is(object, "SECdistrMv")) stop("object of wrong class")
  if(is(object, "SECdistrMv")) { 
     if(slot(object, "family") == "ESN") {
        message("This object is an ESN distribution, passed on to 'SECdistrMv'")
        return(marginalSECdistr(object, comp, name, drop)) }
        else stop("wrong 'family' type of 'SECdistrMv' object")
        } 
  dp <- slot(object, "dp")
  Omega <- dp[[2]]
  d <- dim(Omega)[1]
  if(!all(comp %in% 1L:d)) stop("some comp values not admissible")
  dp.m <- list(xi=dp[[1]][comp], Omega=Omega[comp, comp, drop=FALSE], 
               Delta=dp[[3]][comp,, drop=FALSE], tau=dp[[4]], Gamma=dp[[5]])
  if(missing(name)) {
    comp.c <- paste(as.character(comp), collapse=",")
    name <- paste(slot(object, "name"), "[", comp.c, "]", sep="")       
    }         
  compNames <- slot(object, "compNames")[comp]
  hnames <- slot(object, "HcompNames")
  obj.m <- makeSUNdistr(dp.m, name, compNames, hnames, drop=drop)         
  # if(class(obj.m) != "SUNdistr") stop("Error. No SUNdistr object created")
  obj.m
}

affineTransSUNdistr <- function(object, a, A, name, compNames, HcompNames, drop=TRUE)
{# distribution of affine transformation X=a+t(A)Y; see SN book, top of p.199 
  if(!is(object, "SUNdistr")) stop("wrong object class")
  dp <- slot(object, "dp")
  d <- length(dp$xi)
  if(!is.matrix(A) || nrow(A) != d) stop("A is not a matrix or wrong nrow(A)")
  h <- ncol(A)
  if(length(a) != h) stop("size mismatch of arguments 'a' and 'A'")
  if(missing(name)) name <- paste(deparse(substitute(a)), " + t(",  
    deparse(substitute(A)), ") %*% (", slot(object, "name"),")", sep="")
  else name <- as.character(name)[1]
  if(missing(compNames)) 
    compNames <- as.vector(outer("V",as.character(1:h), paste,sep=""))
  if(missing(HcompNames)) HcompNames <- slot(object, "HcompNames")
  Omega <- dp$Omega
  omega <- sqrt(diag(Omega))
  OmegaX <- t(A) %*% Omega %*% A
  OmegaX <- (OmegaX + t(OmegaX))/2
  eig <- eigen(OmegaX, symmetric=TRUE, only.values=TRUE)$values
  if(any(eig <= 0)) stop("singular transformation") 
  omegaX <- sqrt(diag(OmegaX))
  DeltaA <- (1/omegaX)*t(A) %*% (omega * dp$Delta)
  dpX <- list(xi=as.vector(a + t(A) %*% matrix(dp$xi, ncol=1)), Omega=OmegaX,
              Delta=DeltaA, tau=dp$tau, Gamma=dp$Gamma)
  obj <- makeSUNdistr(dp=dpX, name, compNames, HcompNames, drop=drop)
  return(obj)
}
# 
convolutionSUNdistr <- function(object1, object2,  name, compNames, HcompNames)
{# convolution of two SUN distributions; see SN book eq.(7.8) on p.199
  if(!is(object1, "SUNdistr") | !is(object2, "SUNdistr")) 
    stop("wrong object class")
  dp1 <- slot(object1, "dp")
  dp2 <- slot(object2, "dp")
  m1 <- length(dp1$tau)
  m2 <- length(dp2$tau)
  if(length(dp1$xi) != length(dp2$xi)) stop("objects with different dimensions")  
  name1 <- slot(object1, "name")
  name2 <- slot(object2, "name")
  if(missing(name)) name <- paste("(", name1, ")+(", name2, ")", sep="")
  if(missing(compNames)) compNames <- 
    as.vector(outer("V", as.character(1:length(dp1$xi)), paste, sep=""))
  Omega1 <- dp1$Omega
  omega1 <- sqrt(diag(Omega1))
  Omega2 <- dp2$Omega
  omega2 <- sqrt(diag(Omega2))
  omega <- sqrt(omega1^2+omega2^2)
  Delta <- cbind((omega1/omega)* dp1$Delta,  (omega2/omega)* dp2$Delta )
  if(missing(compNames)) compNames <-
    as.vector(outer("V", as.character(1:length(dp1$xi)), paste, sep=""))
  if(missing(HcompNames))  HcompNames <- 
    c(paste(name1, slot(object1, "HcompNames"), sep="."), 
      paste(name2, slot(object2, "HcompNames"), sep=".")) 
  dimnames(Delta) <- list(compNames, HcompNames)             
  xi <- dp1$xi + dp2$xi     
  names(xi) <- compNames
  Omega <- Omega1 + Omega2
  dimnames(Omega) <- list(compNames, compNames)    
  tau <- c(dp1$tau, dp2$tau)
  names(tau) <- HcompNames
  Gamma <- blockDiag(dp1$Gamma, dp2$Gamma)
  dimnames(Gamma) <- list(HcompNames, HcompNames)    
  dp <- list(xi=xi, Omega=Omega, Delta=Delta, tau=tau, Gamma=Gamma)
  obj <- makeSUNdistr(dp=dp, name, compNames, HcompNames)
  return(obj)
}
#
conditionalSUNdistr <- function(object, comp, values, eventType="=", name, drop=TRUE) 
{# Conditional distribution for the "=" case as given by eq.(7.7) of SN book, and
 # later amendment; the distribution for the ">" case is given by RAV&AA (2020).
  if(!is(object, "SUNdistr")) stop("wrong object class")
  type <- match.arg(eventType, c("=", ">"))
  if(!is.numeric(values)) stop("non-numeric 'values'")
  dp <- slot(object, "dp")
  xi <- dp$xi
  Omega <- dp$Omega
  Delta <- dp$Delta
  tau <- dp$tau
  Gamma <- dp$Gamma
  d <- length(xi)
  m <- length(tau)
  if(!all(comp %in% 1:d)) stop("some 'comp' terms outside range")
  if(length(comp) == d) stop("degenerate conditional distribution")
  if(length(comp) != length(values)) stop("mismatch of comp and values sizes") 
  omega <- sqrt(diag(Omega))
  Omega11 <- Omega[comp, comp, drop=FALSE]
  Omega22 <- Omega[-comp, -comp, drop=FALSE]
  Omega.bar <- cov2cor(Omega)
  if(type == "=") {
    O11.inv <- solve(Omega11)
    tmp1 <- Omega[-comp, comp, drop=FALSE] %*% O11.inv 
    values0 <- matrix(values - xi[comp], ncol=1)
    xi2.1 <- c(xi[-comp] + tmp1 %*% values0)
    O22.1 <- Omega22 - tmp1 %*% Omega[comp, -comp, drop=FALSE]
    tmp2 <- solve(Omega.bar[comp, comp, drop=FALSE])
    Delta1 <- Delta[comp, , drop=FALSE]
    Delta2 <- Delta[-comp, , drop=FALSE]
    tau2.1 <- c(tau + t(Delta1) %*% tmp2 %*% (values0/omega[comp]))
    Delta2.1 <- Delta2 - Omega.bar[-comp,comp] %*% tmp2 %*% Delta1
    Gamma2.1 <- Gamma - t(Delta1) %*% tmp2 %*% Delta1
    s <- sqrt(diag(Gamma2.1))
    sDelta <- Delta2.1 %*% diag(1/s, m, m) 
    stau <- tau2.1/s
    sGamma <- cov2cor(Gamma2.1)
    if(missing(name)) name <- paste(slot(object, "name"), "|comp[",
      paste(comp,collapse=","), "]=(", paste(format(values), collapse=","), 
      ")", sep="")
    names <- slot(object, "compNames")[-comp]  
    dp.c <- list(xi=xi2.1, Omega=O22.1, Delta=sDelta, tau=stau, Gamma=sGamma)
    hnames <- slot(object, "HcompNames")
    obj <- makeSUNdistr(dp=dp.c, name, names, hnames, drop=drop)
    }
  if(type == ">") {
    xi.c <- xi[-comp]
    Delta.c <- cbind(Delta[-comp,, drop=FALSE], Omega.bar[-comp,comp,drop=FALSE])
    tau.c <- c((xi[comp] + (-values))/omega[comp], tau)
    Gamma.c <- rbind(cbind(Omega.bar[comp, comp, drop=FALSE], Delta[comp,,drop=FALSE]), 
                     cbind(t(Delta[comp,, drop=FALSE]), Gamma))
    dp.c <- list(xi=xi.c, Omega=Omega22, Delta=Delta.c, tau=tau.c, Gamma=Gamma.c)
    if(missing(name)) name <- paste(slot(object, "name"), "|comp[",
      paste(comp, collapse=","), "]>(", paste(format(values), collapse=","), 
      ")", sep="")
    names <- slot(object, "compNames")[-comp] 
    hnames <- c(slot(object, "compNames")[comp], slot(object, "HcompNames"))
    obj <- makeSUNdistr(dp=dp.c, name, names, hnames, drop=drop)           
    }
  return(obj)     
}
#
joinSUNdistr <- function(object1, object2, name, compNames, HcompNames) 
{# join two SUN distributions assuming independence
  obj1 <- object1
  obj2 <- object2
  if(!is(obj1, "SUNdistr")) obj1 <- convertSN2SUNdistr(obj1, silent=TRUE)
  if(is.null(obj1)) stop("object1 is neither a SUNdistr object nor adjustable") 
  if(!is(obj2, "SUNdistr")) obj2 <- convertSN2SUNdistr(obj2, silent=TRUE)
  if(is.null(obj2)) stop("object2 is neither a SUNdistr object nor adjustable") 
  dp1 <- slot(obj1, "dp")
  dp2 <- slot(obj2, "dp")
  name1 <- slot(obj1, "name")
  name2 <- slot(obj2, "name")
  if(missing(name)) name <- paste("(",name1, ")x(", name2, ")", sep="")
  if(missing(compNames))  compNames <- 
    c(paste(name1, slot(obj1, "compNames"), sep="."), 
      paste(name2, slot(obj2, "compNames"), sep="."))
  if(missing(HcompNames))  HcompNames <- 
    c(paste(name1, slot(obj1, "HcompNames"), sep="."), 
      paste(name2, slot(obj2, "HcompNames"), sep="."))     
  dp <- list(xi=c(dp1$xi, dp2$xi), Omega=blockDiag(dp1$Omega, dp2$Omega),
            Delta=blockDiag(dp1$Delta, dp2$Delta), tau=c(dp1$tau, dp2$tau),
            Gamma=blockDiag(dp1$Gamma, dp2$Gamma))
  makeSUNdistr(dp, name, compNames, HcompNames)          
}

convertSN2SUNdistr <- function(object, HcompNames="h", silent=FALSE)
{# converts SN/ESN into a SUN distribution 
  # obj.cl <- class(object)
  if(!is(object, "SECdistrUv") & !is(object, "SECdistrMv"))
    if(silent) return(NULL) else stop("wrong class object")
  obj.fm <- slot(object, "family")    
  if(!(obj.fm %in% c("SN", "ESN"))) 
    if(silent) return(NULL) else stop("wrong family of distributions")
  dp <- slot(object, "dp")   
  if(is(object, "SECdistrUv")) {
    xi <- dp[1]
    Omega <- matrix(dp[2]^2, 1, 1)
    alpha <- dp[3]
    Delta <- matrix(alpha/sqrt(1+alpha^2), 1, 1)
    tau <- if(length(dp)>3) dp[4] else 0
    names <- slot(object, "name")
    }
  if(is(object, "SECdistrMv")) {  
    xi <- dp[[1]]
    Omega <- dp[[2]]
    alpha <- dp[[3]]
    etc <- delta.etc(alpha, Omega)
    Delta <- matrix(etc$delta, ncol=1)
    tau <- if(length(dp)>3) dp[[4]] else 0
    names <- slot(object, "compNames")
    }
  dp <- list(xi=xi, Omega=Omega, Delta=Delta, tau=tau, Gamma=matrix(1, 1, 1))
  makeSUNdistr(dp=dp, slot(object, "name"), names, HcompNames[1], drop=FALSE)        
}

convertCSN2SUNpar <- function(mu, Sigma, D, nu, Delta)
{# convert a set of CSN parameters to their SUN equivalents
  if(!all.numeric(mu, Sigma, D, nu, Delta)) stop("non-numeric argument(s)")
  if(any(eigen(Sigma, only.values=TRUE)$values <= 0)) stop("invalid Sigma")
  if(any(eigen(Delta, only.values=TRUE)$values <= 0)) stop("invalid Delta")
  p <- NCOL(Sigma)
  q <- NCOL(Delta)
  if(length(mu) != p) stop("mismatch of dimensions")
  if(length(nu) != q) stop("mismatch of dimensions")
  if(any(dim(D) != c(q, p))) stop("mismatch of dimensions")
  DS <- D %*% Sigma
  M <- rbind(cbind(Sigma, t(DS)),  cbind(DS, Delta + DS %*% t(D)))
  M0 <- cov2cor(M)
  Gamma <- M0[p + (1:q), p + (1:q), drop=FALSE]
  DeltaSUN <- M0[1:p, p+(1:q), drop=FALSE]
  list(xi=mu, Omega=matrix(Sigma, p, p), Delta=DeltaSUN, tau=-nu, Gamma=Gamma)
}

#----------------
summary.SUNdistr <- function(object, ...) 
{# 
  dp <- slot(object, "dp")
  name <- slot(object, "name")
  compNames <- slot(object, "compNames")
  HcompNames <- slot(object, "HcompNames")
  u <- sunValues(dp=dp, compNames, HcompNames, ...)
  new("summary.SUNdistr", dp=dp, name=name, compNames=compNames,
    HcompNames=HcompNames, mean=u$mean, var.cov=u$vcov, gamma1=u$gamma1, 
    cum3=u$cum3,  mardia=u$mardia)
}

sunValues <- function(dp, compNames, HcompNames, silent=FALSE, ...) 
{# Some moments and other characteristics values of a SUN distribution. 
 # Computations are based on Proposition 1 and 2 of RAV&AA-2020
 # This function is *not* exported in NAMESPACE.
  if (length(dp) != 5)  stop("wrong length of non-null 'dp'")
  xi <- drop(dp[[1]])
  Omega <- dp[[2]]
  Delta <- dp[[3]]
  tau <-  dp[[4]]
  Gamma <- dp[[5]]
  if(!all.numeric(xi, Omega, Delta, tau, Gamma)) stop("non-numeric argument(s)")   
  d <- length(xi)
  m <- length(tau)
  if(m > 20) 
    {if(silent) return(NA) else stop("m exceeds the admissible size")}  
  if(missing(compNames)) compNames <- paste("V", 1:d, sep="")
  if(missing(HcompNames)) HcompNames <- paste("H", 1:m, sep="")
  omega <- sqrt(diag(Omega))
  omega.Delta <- omega * Delta
  Gamma.inv <- solve(Gamma)
  A <- omega.Delta %*% Gamma.inv            # A=\Lambda in (17) of RAV&AA-2020
  mom.U <- mnormt::mom.mtruncnorm(4, mean=rep(0,m), Gamma, lower=-tau, ...)
  E.U <- if(m==1) mom.U$cum[1] else mom.U$cum1
  mu1.X <- A %*% E.U                        # see \mu_1(X) in Proposition 1
  Esun <- dp$xi + drop(mu1.X)           
  names(Esun) <- compNames
  # E.U2 <- if(m==1) mom.U$mom[3] else mom.U$order2$m2
  var.U <- if(m==1) mom.U$cum[2] else mom.U$order2$cum2
  Vsun <- Omega - A %*% (Gamma- var.U) %*% t(A)    # see var(X) in Prop.1
  dimnames(Vsun) <- list(compNames, compNames)    
  Sigma <- var.X <- Vsun
  sigma <- sqrt(diag(Sigma))
  Sigma.inv <- solve(Sigma)  
  mu2.X  <- var.X + mu1.X %*% t(mu1.X)
  #---
  # Calcolo cumulanti/momenti centrali del terzo ordine.
  # Partiamo da \mu_3(X) della Proposizione 1 di RAV&AA-2020.
  # Calcoliamo (I_{d^2}+K_d) utilizzando eqn.(4) e (7) a p.57
  # di Magnus & Neudecker (2007, 3^ ed) 
  D <- duplicationMatrix(d)
  Dplus <- solve(t(D) %*% D) %*% t(D)
  twiceD.Dplus_d <-  2*D %*% Dplus
  B.BT <- Omega - omega.Delta %*% Gamma.inv %*% t(omega.Delta)  # \Psi in (17)  
  mu3.U <- if(m==1) mom.U$mom[3+1] else array2mat(mom.U$order3$m3, m) 
  mu3.X <- ( (A %x% A) %*% mu3.U %*% t(A) + twiceD.Dplus_d %*% (mu1.X %x% B.BT)           
            + matrix(B.BT, ncol=1) %*% t(mu1.X) )       
  # Now apply shift \xi=-mu1.X in (A.7) of RAV&AA; first two terms cancel out
  shift <- (-mu1.X)
  cum3 <- (twiceD.Dplus_d %*% (shift %*% t(shift) %x% mu1.X + shift %x% mu2.X) 
           + matrix(mu2.X, ncol=1) %*% t(shift) + mu3.X)   
  cum3 <- array(cum3, dim=c(d,d,d))                # convert matrix into array
  gamma1 <- cum3[cbind(1:d, 1:d, 1:d)]/sigma^3  
  #---
  # Mardia measures of skewness and kurtosis; use Proposition 2 of RAV&AA-2020
  AA <- Gamma.inv %*% t(omega.Delta) %*% Sigma.inv %*% omega.Delta %*% Gamma.inv
  # AA =\tilde\Lambda^T\tilde\Lambda = \Lambda^T\Sigma\inv\Lambda in Prop.2
  vec.mu3 <- if(m==1) mom.U$centr.mom[3] else c(mom.U$order3$cum3)
  gamma1M <- beta1M <- if(m==1) drop(vec.mu3^2 *AA^3) else
                           drop(t(vec.mu3) %*% (AA %x% AA %x% AA) %*% vec.mu3)
  mu4.U <- if(m==1) mom.U$centr.mom[4] else 
   { cum4 <- array2mat(mom.U$order4$cum4, m) # conversione cum4 in matrice
    # Usiamo (2.8)-(2.9) di Kollo & Srivastava (2005, Comms.Stat-TM) per
    # passare da cumulanti a momenti centrali del quarto ordine, con correzione!
    D <- duplicationMatrix(m)
    Dplus <- solve(t(D) %*% D) %*% t(D)
    twiceD.Dplus_m <-  2*D %*% Dplus
    cmom4N <- twiceD.Dplus_m %*% (var.U %x% var.U) + c(var.U) %*% t(c(var.U))  
    cum4 + cmom4N  # matrice dei quarti momenti centrali
    }
  tmp1 <- Gamma.inv %*% t(omega.Delta) %*% Sigma.inv
  tmp2 <- B.BT %*% Sigma.inv
  beta2M <- ( tr((AA %x% AA) %*% mu4.U)  + 2* tr(var.U %*%  AA) * tr(tmp2) 
    + tr(tmp2)^2 + 4 * tr(var.U %*%  tmp1 %*% B.BT %*% t(tmp1))
    + 2 * tr(tmp2 %*% tmp2) )
  mardia <- c(gamma1M=gamma1M, gamma2M=(beta2M-d*(d+2)))
  list(mean=Esun, vcov=Vsun, gamma1=gamma1, cum3=cum3, mardia=mardia)
}  
#----------------
# plotting SUN densities
#      
plot.SUNdistr <- function(x, range, nlevels=8, levels, npt,  main, comp, 
  compLabs, gap = 0.5, ...) 
{# plot density of object of class "SUNdistr" 
  obj <- x
  if(slot(obj, "class") != "SUNdistr") stop("object of wrong class")
  dp <- slot(obj, "dp")
  d <- length(dp$xi)
  if(missing(comp)) comp <- seq(1, d) 
  if(!all(comp %in% seq(1,d))) stop("illegal 'comp' value(s)")
  pd <- length(comp) # actual plotting dimension
  if(missing(npt)) npt <- if(pd==1) 251 else rep(101, pd)
  pobj <- if(pd == d) obj else marginalSUNdistr(obj, comp=comp, drop=FALSE)  
  name.pobj <- slot(obj, "name")
  if(pd < d) name.pobj <- paste(name.pobj,"[", paste(comp, collapse=","), "]", sep="")
  if(missing(main))  { main <- if(pd == 1 | pd == 2) 
      paste("Density function of", name.pobj) else
      paste("Bivariate densities of", name.pobj) 
      }
  compNames <- slot(pobj, "compNames")
  if(missing(compLabs)) compLabs <- compNames    
  if(length(compLabs) != pd) stop("wrong length of 'compLabs' vector")
  if(missing(range)) {
    range <- matrix(NA, 2, pd)
    dp.pobj <- slot(pobj, "dp")
    m <- sunMean(dp=dp.pobj)
    v <- sunVcov(dp=dp.pobj)
    s <- sqrt(diag(v))
    range <- rbind(m -3*s, m + 3*s)
    }
  dots <- list(...)
  nmdots <- names(dots)    
  if(pd == 1)   out <- plot.SUNdistrUv(pobj, range, npt, main, ...)
  if(pd == 2) {
    p <- plot.SUNdistrBv(pobj, range, nlevels, levels, npt, compLabs, main, ...)
    out <- list(object=pobj, plot=p)
    }            
  if(pd > 2) {
    textPanel <- function(x = 0.5, y = 0.5, txt, cex, font) 
      text(x, y, txt, cex = cex, font = font)
    localAxis <- function(side, x, y, xpd, bg, main,  oma, ...) {
      if (side%%2 == 1) Axis(x, side = side, xpd = NA,  ...) else 
      Axis(y, side = side, xpd = NA, ...)
      }
    localPlot <- function(..., oma, font.main, cex.main) plot.SUNdistrBv(...)
    text.diag.panel <- compLabs
    oma <- if ("oma" %in% nmdots) dots$oma else NULL
    if (is.null(oma)) {
      oma <- c(4, 4, 4, 4)
      if (!is.null(main))  oma[3L] <- 6
      }    
    opar <- par(mfrow = c(length(comp), length(comp)), 
                mar = rep(c(gap,gap/2), each=2), oma=oma)
    on.exit(par(opar))
    out <- list(object=pobj)
    count <- 1
    for (i in comp) 
      for (j in comp) {
        count <- count + 1
        if(i == j) {
          plot(1, type="n", xlab="", ylab="", axes=FALSE)
          text(1, 1, text.diag.panel[i], cex=2)
          box()
          out[[count]] <- list()
          names(out)[count] <- paste("diagonal component", compNames[i])
        } else {
        ji <- c(j,i) 
        marg <- marginalSUNdistr(pobj, comp=ji, drop=FALSE) 
        out[[count]] <- localPlot(x=marg, range=range[,ji], nlevels, levels,
             npt=npt[ji], compNames= compNames[ji], compLabs=compLabs[ji],              
             main="", yaxt="n", xaxt="n", ...)   
        names(out)[count] <- paste("plot of components (", j, ",", i, ")")
        if(i==comp[1]) axis(3) ; if(j==length(comp)) axis(4)
        if(j==comp[1]) axis(2) ; if(i==length(comp)) axis(1)    
        box() }
      }
    par(new = FALSE)
    if (!is.null(main)) {
      font.main <- if ("font.main" %in% nmdots) 
         dots$font.main else par("font.main") 
      cex.main <- if ("cex.main" %in% nmdots) 
         dots$cex.main  else par("cex.main") 
      mtext(main, side=3, TRUE, line=5, outer = TRUE, at=NA, cex=cex.main, 
            font=font.main, adj=0.5)
      }}
  invisible(out)
}

plot.SUNdistrBv <- 
  function(x, range, nlevels=8, levels, npt, compLabs, main, ...)
{# plot BiVariate SUN distribution (hence d=2)
  obj <- x
  if(slot(obj, "class") != "SUNdistr") stop("object of wrong class")
  dp <- slot(obj, "dp")
  d <- length(dp[[1]])
  if(d != 2) stop("wrong dimensions, d=2 is required")
  if(missing(npt)) npt <- rep(51, d)
  n1 <- npt[1]
  n2 <- npt[2]
  x1 <- seq(min(range[,1]), max(range[,1]), length=n1)
  x2 <- seq(min(range[,2]), max(range[,2]), length=n2)
  x1.x2 <- cbind(rep(x1, n2), as.vector(matrix(x2, n1, n2, byrow=TRUE)))
  X <- matrix(x1.x2, n1 * n2, 2, byrow = FALSE)
  pdf <- matrix(dsun(X, dp=dp), n1, n2)
  oo <- options()
  options(warn=-1)
  compNames <- slot(obj ,"compNames")
  if(missing(levels)) levels <- pretty(range(pdf, finite=TRUE), nlevels)[-1]
  if(missing(compLabs)) compLabs <- compNames
  contour(x1, x2, pdf, levels=levels, labels=format(levels),
    main=main, xlab=compLabs[1], ylab=compLabs[2], ...)
  options(oo) 
  cL <- contourLines(x1, x2, pdf, levels=levels)
  for(j in 1:length(cL)) cL[[j]]$level <- levels[j]
  return(list(x=x1, y=x2, names=compNames, density=pdf, contourLines=cL))
}    

plot.SUNdistrUv <- function(x, range, npt=251, main, ...)
{# plot density of object "SUNdistr" when d=1
  obj <- x
  if(slot(obj, "class") != "SUNdistr") stop("object of wrong class")
  dp <- slot(obj, "dp")
  if(length(dp[[1]]) != 1) stop("SUN distribution of wrong dimension")
  dots <- list(...)
  nmdots <- names(dots)  
  topline <- if(obj@name == "") "" else	
    paste("Probability density of ", obj@name, "\n", sep="")
  if(missing(main))  main <- paste(topline, "\nunivariate SUN distribution")
  mar <- if ("mar" %in% nmdots) dots$mar else NULL
  if (is.null(mar)) {
    mar <- c(4.5, 4.5, 4, 2)
    if (is.null(main))  mar[3L] <- 2
    }     
  omar <- par()$mar
  on.exit(par(omar))
  par(mar=mar)    
  x <- seq(min(range), max(range), length=npt)
  pdf <- as.vector(dsun(matrix(x, ncol=1), dp=dp))
  xLab <- if("xlab" %in% nmdots) dots$xlab else slot(obj, "name")
  yLab <- if("ylab" %in% nmdots) dots$ylab else "probability density"
  yLim <- if("ylim" %in% nmdots) dots$ylim else c(0, max(pdf))
  plot(x, pdf, type="n", xlab=xLab, ylab=yLab, ylim=yLim)
  lines(x, pdf, ...)
  abline(h=0, lty=2, col="gray50")
  if (!is.null(main)) {
    font.m <- if("font.main" %in% nmdots) dots$font.main else par("font.main") 
    cex.m <- if("cex.main" %in% nmdots) dots$cex.main else par("cex.main") 
    title(main, line=2, cex.main=cex.m, font.main=font.m)
    }    
  invisible(list(object=obj, x=x, density=pdf))
}
#============================ classes and methods ============================

setClass("SUNdistr",
   representation(dp="list", name="character",  compNames="character", 
     HcompNames="character"),
   validity=function(object){     
     dp <- object@dp
     if(length(dp) != 5) return(FALSE)
     if(!all(names(dp) == c("xi", "Omega", "Delta", "tau", "Gamma"))) return(FALSE)
     if(mode(unlist(dp)) != "numeric") return(FALSE)
     if(!is.character(object@name)) return(FALSE)
     if(length(object@name) != 1) return(FALSE)
     if(length(object@compNames) != length(dp[[1]])) return(FALSE)
     if(length(object@HcompNames) != length(dp[[4]])) return(FALSE)
     # numeric checks are assumed to be handled by makeSUNdistr
     TRUE
     }
)

setMethod("show", "SUNdistr",
  function(object){
    if(!is(object, "SUNdistr")) stop("wrong object class")
    if(object@name != "")  
      cat("Probability distribution of variable '", object@name, "'\n", sep="")
    dp <- slot(object, "dp")
    d <- length(dp[[1]])
    m <- length(dp[[4]])
    compNames <- slot(object, "compNames")
    HcompNames <-  slot(object, "HcompNames")
    cat("This is a SUN distribution of dimension d=", d, ", involving m=",
         m, " hidden variables:", sep="")       
    cat("\n\nd-component parameters (xi, Omega):\n")
    out <- rbind(xi=dp$xi, Omega=dp$Omega)
    rownames(out) <- c("xi", paste("Omega[", compNames, ",", sep=""))
    colnames(out) <- compNames
    print(out)
    cat("\nm-component parameters (Delta, tau, Gamma):\n")
    out <- rbind(dp$Delta, dp$tau, dp$Gamma) 
    rownames(out) <- c( paste("Delta[", compNames, ",", sep=""), "tau", 
           paste("Gamma[", HcompNames, ",", sep=""))
    colnames(out) <- HcompNames       
    print(out)       
  }
)

setClass("summary.SUNdistr",
  representation(dp="list", name="character", 
    compNames="character", HcompNames="character",
    mean="vector", var.cov="matrix", gamma1="vector", cum3="array",
    mardia="vector"),
  validity=function(object) {
    dp <- slot(object, "dp")
    if(length(dp) != 5) return(FALSE)
    if(mode(unlist(dp)) != "numeric") return(FALSE)
    d <- length(dp[[1]])
    # m <- length(dp[[4]])
    if(length(slot(object, "mean")) != d) return(FALSE)
    if(any(dim(slot(object, "var.cov")) != c(d,d))) return(FALSE)
    if(length(slot(object, "gamma1")) != d) return(FALSE)
    if(any(dim(slot(object, "cum3")) != c(d,d,d))) return(FALSE)
    if(length(slot(object, "mardia")) != 2) return(FALSE)
    TRUE
  }
)     

setMethod("show", "summary.SUNdistr",
  function(object){
    obj <- object
    dp <- slot(obj, "dp")
    sun <- new("SUNdistr", dp=dp, name=slot(obj, "name"), 
      compNames=slot(obj, "compNames"), HcompNames=slot(obj, "HcompNames"))
    show(sun)
    cat("\nExpected value:\n")
    print(slot(obj, "mean"))
    cat("\nVariance matrix:\n")
    print(slot(obj, "var.cov"))
    cat("\nCoefficients of marginal skewness (gamma1):\n")
    print(slot(obj, "gamma1"))
    cat("\nMardia's measures of multivariate skewness and kurtosis:\n")
    print(slot(obj, "mardia"))
  }
)   

setMethod("plot", signature(x="SUNdistr", y="missing"), plot.SUNdistr)
setMethod("mean", signature(x="SUNdistr"), mean.SUNdistr)
setMethod("vcov", signature(object="SUNdistr"), vcov.SUNdistr)
setMethod("summary", signature(object="SUNdistr"), summary.SUNdistr)

