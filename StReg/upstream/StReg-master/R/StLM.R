StLM <-
function(y,X,Trend=1,v=1,maxiter=1000,meth="BFGS",hes="FALSE",init="na")
{
X <- cbind(X)
Data <- cbind(y,X) # Data matrix

# Error Handling: Error message if y and X do not have same number of rows.
if (length(y)!= nrow(cbind(X))) stop("y and X matrix must have same number of rows")

Z <- t(cbind(y,X)) ; T <- nrow(cbind(X)) ; l <- ncol(cbind(y,X))  ## Data matrix transposed; Sample size and number of variables in VAR.
if (mean(Trend)==1 & length(Trend)==1) {Trend <- matrix(1,T,1) ; colnames(Trend) <- "const"; TREND <- Trend ; c <- 1} # If 1, then matrix of 1s is created for intercept term.

if (mean(Trend)==0 & length(Trend)==1) {Trend <- matrix(0,T,1); TREND <- Trend ; c <- 0} # If 0, then intercept is excluded from the model.

if(ncol(cbind(Trend))>1) {TREND <- Trend ; c <- ncol(Trend)} # If more than one column is supplied in Trend, then use the whole matrix as mean (t varying intercept).
if(Trend[1]!=1) {if (is.null(colnames(Trend))=="TRUE") {colnames(Trend) <- paste("t",".",seq(1,ncol(Trend),1),sep="")}} # If Trend matrix does not have column names provided, then automatically created.

# Error Handling: Error message if Trend matrix and Data matrix do not have same number of rows.
if (Trend[1]!=1) {if (nrow(cbind(X))!= nrow((Trend))) stop("Data and Trend matrix must have same number of rows")}

# Error Handling: Error message if X matrix with at least one column is not included.
if (ncol(cbind(X))<1) stop("no. of explanatory variables in X must be at least 1")

# Error Handling: Error message if Trend matrix with at least one column is not provided
if (Trend[1]!=1) {if (length(Trend)==1) stop("The Trend matrix must have at least one column")}

# Error Handling: Error message if degrees of freedom is not a positive number
if (v<=0) stop("degrees of freedom must be a positive number")

# Error Handling: Error message if maximum iterations for numerical optimization is less than 10
if (maxiter<10) stop("Iteration must be at least 10")

# Error Handling: Error message if number of parameters is too large. Reduce the variables or increase the sample size.
if (2*T<3*((3)*(l+1)*l/2)+3*l+c) stop(list(((3)*(l+1)*l/2)+3*l+c,"Too many parameters for given sample size."))

##Likelihood Function
L<-function(a)
{
 S <- BlockTop(a[1:((3)*(l+1)*l/2)],l,lag=0)$S # Var-Cov matrix created using block Toipletz matrix from BlockTop function
 F <- solve(S)/v # VarCov and its Inverse

if(c!=0) {MU <- matrix(c(a[(((3)*(l+1)*l/2)+1):(((3)*(l+1)*l/2)+l*c)]),l,c)} # Mean if intercept is to be included in the model

if(c!=0) {M0 <- kronecker(diag(1,1),MU)} ## Kronecker product of the Mean.
if(c==0) {M0 <- kronecker(diag(1,1),0)} ## Kronecker product of the Mean.

if(c!=0) {M <- M0%*%t(TREND)} ; if (c==0) {M <- 0}

D <- 1 + diag(t(Z-M)%*%F%*%(Z-M)) # Quadratic form of the Student's t likelihood function

## Likelihood function
 LLn <- (T)*const - 0.5*(T)*log(det(S)) - 0.5*(v+(1)*l)*sum(log(D))
 neg.LLn <- -LLn # negative of likelihood function to be minimized.
 neg.LLn
}

if(init[1]=="na") int <- c(runif(((3)*(l+1)*l/2)+l*c,0,.1))  # Initialization of the parameters if not provided

if(init[1]!="na") int <- init # If initial values are provided
const <- log(gamma((v+(1)*l)/2)) - log(gamma(v/2)) - 0.5*l*(1)*log(pi*v) # constant part of the log likelihood function.

# Optimization starts
op <- optim(int,L,hessian=hes,control=list(trace=1,maxit=maxiter,reltol=1e-14),method=meth)
a <- op$par ; Like <- op$value ; hess <- op$hes # Gather the optimization results

################################### Parameters Estimated ##################
S <- BlockTop(a[1:((3)*(l+1)*l/2)],l,lag=0)$S
if(c!=0) {MU <- matrix(c(a[(((3)*(l+1)*l/2)+1):(((3)*(l+1)*l/2)+l*c)]),l,c)}

if(c!=0) {M0 <- kronecker(diag(1,1),MU)}
if(c==0) {M0 <- kronecker(diag(1,1),0)}

if(c!=0) {M <- M0%*%t(TREND)} ; if (c==0) {M <- 0}

PP <- Par.dlrm(a,S,MU,l,0,v,T,M,c,TREND) # Estimation of parameters of conditional distribution
Delta <- PP$Delta ; Delta0 <- PP$Delta0
B1 <- PP$B1 ; s2 <- PP$s2 ; Q <- PP$Q
var.coef <- PP$var.coef

beta <- cbind(Delta0, t(B1))
colnames(beta) <- c(colnames(Trend),colnames(Data)[-1]) # coefficients are given names


############# Fitted values/Residuals/Conditional Covariance ###############
q <- v/(v+l*0+l-2)
Ct <- vector(length=(nrow(cbind(X))))
U <- muy <- vector(length=(nrow(cbind(X))))
trend <- vector(length=ncol(cbind(X))-1)
for(i in 1:(nrow(cbind(X))))
{
  if(c!=0) {Ct[i] <- 1 + t(Z[2:ncol(S),i]-M[2:ncol(S)])%*%Q%*%(Z[2:ncol(S),i]-M[2:ncol(S)])}
  if(c==0) {Ct[i] <- 1 + t(Z[2:ncol(S),i])%*%Q%*%(Z[2:ncol(S),i])}
  if(c!=0) {trend[i] <- Delta0%*%(Trend[i,])} ; if(c==0) {trend[i] <- 0%*%(Trend[i,])}
  muy[i] <- trend[i] + t(B1)%*%Z[c(2:l,rep(1:l,0)),i]
  U[i] <- t(Z[1,i])- muy[i]
}

#############################################################
if(hes=="TRUE") VARth <- solve(op$hessian) # Inversion of the hessian matrix
############################ Coefficients/SEs/P-values########
################
##Jacobian##SE##
################

if(hes=="TRUE")
{
Jc <- Jacob.dlrm(a,0,l,v,c)$J #jacobian(J,a) ## Computing the Jacobian matrix
SE <- sqrt(diag(Jc%*%VARth%*%t(Jc)))
p_value <- 2*(1-pt(abs(c(Delta0,t(B1),vech(s2)))/SE,(T)))
COEF <- round(cbind(c(Delta0,t(B1),vech(s2)),c(SE[1:length(SE)]),c(p_value[1:length(SE)])),8)
colnames(COEF) <- c("coef.","std.err.","p-value")

Jv <- ConJacob.dlrm(a,0,l,v,c)$Jv ## Computing the Jacobian matrix
SEv <- sqrt(diag(Jv%*%VARth%*%t(Jv)))
p_valuev <- 2*(1-pt(abs(c(var.coef))/SEv,(T)))
VAR.COEF <- round(cbind(c(var.coef),c(SEv[1:length(SEv)]),c(p_valuev[1:length(SEv)])),8)
colnames(VAR.COEF) <- c("var.coef","std.err.","p-value")

B <- COEF[1:ncol(beta),] ; rownames(B) <- colnames(beta)
}


if(hes=="FALSE") {
COEF <- round(cbind(c(Delta0,B1,vech(s2))),8) ; colnames(COEF) <- c("coef.")
VAR.COEF <- round(cbind(var.coef),8) ; colnames(VAR.COEF) <- c("var.coef")
B <- beta
}

if(Trend[1]!= 0) beta <- cbind(Delta0,t(B1))
if(Trend[1] == 0) beta <- cbind(t(B1))

rownames(beta) <- paste(colnames(Data)[1])

head.x <- matrix(nrow=0,ncol=ncol(Data))
#if(lag!=0) for(i in 1:lag){head.x[i,] <- paste(colnames(Data),".",i,sep="")}

ifelse(mean(Trend[,1])==0 & ncol(Trend)==1, colnames(beta) <- c(colnames(X)), ifelse(mean(Trend)==1 & ncol(Trend)==1,colnames(beta) <- c("const.",colnames(X),t(head.x)),colnames(beta) <- c("const.",colnames(Trend)[2:ncol(Trend)],colnames(X))) )

head.var.coef <- matrix(nrow=(ncol(X)+1),ncol=0)
#if(lag!=0) for(i in 1:lag){head.var.coef[,i] <- c(paste(colnames(Data),".",i,sep=""))}

#if(lag!=0) coef.label <- c(colnames(X),head.var.coef)
coef.label <- c(colnames(X))

var.head <- matrix(nrow=length(coef.label),ncol=length(coef.label))
for(i in 1:length(coef.label))
{
for(j in 1:length(coef.label))
{
var.head[i,j] <- paste(coef.label[j],".",coef.label[i],sep="")
}
}

res <- cbind(U)

rownames(VAR.COEF) <- c("const.",vech(var.head))

ad <- Student(res,Data,Ct,s2,v,0,B1)$ad # Anderson Darling test for Student's t test

p <- ifelse(hes=="FALSE",ncol(beta),ncol(beta))

fitted <- muy
R2 <- sum((fitted-mean(Data[((1):nrow(Data))]))^2 )/sum((Data[((1):nrow(Data))]-mean(Data[((1):nrow(Data))]))^2) # R-squared
F <- (R2/(1-R2))*(nrow(Data)-p)/(p-1) # F-statistics

if(hes=="TRUE") result <- list(beta=B,var.coef=VAR.COEF,like=-Like,sigma=cbind("s^2"=s2, "stderr.sigma"=tail(COEF[,2],length(vech(s2)))),cvar=Ct,trend=trend,res=res,fitted=muy,init=a,hes=hess,S=S,R.squared=R2,F.stat=F,ad=ad)
if(hes=="FALSE") result <- list(beta=B,var.coef=VAR.COEF,like=-Like,sigma=s2,cvar=Ct,trend=trend,res=res,fitted=muy,init=a,S=S,R.squared=R2,F.stat=F,ad=ad)

return(result)
}
