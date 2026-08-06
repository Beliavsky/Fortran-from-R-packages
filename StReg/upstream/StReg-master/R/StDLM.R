StDLM <-
function(y,X,Trend=1,lag=1,v=1,maxiter=1000,meth="BFGS",hes="FALSE",init="na")
{
Data <- cbind(y,X) # Data matrix
if (lag-round(lag)!=0) stop("lag must be an integer") # Error Handling: Lag must be a positive integer
if (lag<=0) stop("lag number must be positive integer") # Error Handling: lag must be a positive integer

if (length(y)!= nrow(cbind(X))) stop("y and X matrix must have same number of rows")

Z <- t(embed(cbind(y,X),lag+1)) ; T <- nrow(cbind(X)) ; l <- ncol(cbind(y,X))  ## Number of variables in VAR
if (mean(Trend)==1 & length(Trend)==1) {Trend <- matrix(1,T,1); TREND <- embed(Trend,lag+1) ; c <- 1} # If 1 is then matrix of 1s is created for intercept term.

if (mean(Trend)==0 & length(Trend)==1) {Trend <- matrix(0,T,1); TREND <- embed(Trend,lag+1) ; c <- 0} # If 0, then intercept is excluded from the model.

if(ncol(cbind(Trend))>1) {TREND <- embed(Trend,lag+1) ; c <- ncol(Trend)}
if(Trend[1]!=1) {if (is.null(colnames(Trend))=="TRUE") {colnames(Trend) <- paste("t",".",seq(1,ncol(Trend),1),sep="")}}

# Error Handling: Error message if Trend matrix and Data matrix do not have same number of rows. 
if (Trend[1]!=1) {if (nrow(cbind(X))!= nrow((Trend))) stop("Data and Trend matrix must have same number of rows")}

# Error Handling: Number of variables must of at least one in X
if (ncol(cbind(X))<1) stop("number of variables in X must be at least 1")

# Error Handling: the trend matrix must have at least one column.
if (Trend[1]!=1) {if (ncol(Trend)<1) stop("The Trend matrix must have at least one column")}

# Error Handling: Error messaage if degrees of freedom is not a positive number
if (v<=0) stop("degrees of freedom must be a positive integer")

# Error Handling: Error messaage if maximum iterations for numerical optimization is less than 10
if (maxiter<10) stop("Iteration must be at least 10")

# Error Handling: Error messaage if number of parameters is too large. Reduce the variables or increase the sample size.
if (2*T<3*((lag+3)*(l+1)*l/2)+3*l+c) stop(list(((lag+3)*(l+1)*l/2)+3*l+c,"Too many parameters for given sample size. Reduce the number of lags."))

##Likelihood Function
L<-function(a)
{
 S <- BlockTop(a[1:((lag+3)*(l+1)*l/2)],l,lag)$S #Var-Cov Matrix created using block Toipletz matrix from BlockTop function
 F <- solve(S)/v ## VarCov and its Inverse

if(c!=0) {MU <- matrix(c(a[(((lag+3)*(l+1)*l/2)+1):(((lag+3)*(l+1)*l/2)+l*c)]),l,c)}

if(c!=0) {M0 <- kronecker(diag(1,lag+1),MU)} ## Kronecker product of the Mean.
if(c==0) {M0 <- kronecker(diag(1,lag+1),0)} ## Kronecker product of the Mean. 

if(c!=0) {M <- M0%*%t(TREND)} ; if (c==0) {M <- 0}

D <- 1 + diag(t(Z-M)%*%F%*%(Z-M)) # Quadratic form of the Student's t likelihood function

## Likelihood function
 LLn <- (T-lag)*const - 0.5*(T-lag)*log(det(S)) - 0.5*(v+(lag+1)*l)*sum(log(D))
 neg.LLn <- -LLn # negative of likelihood function to be minimized.
 neg.LLn
}

if(init[1]=="na") int <- c(runif(((lag+3)*(l+1)*l/2)+l*c,0,.1))  # Initialization of the parameters if not provided

if(init[1]!="na") int <- init
const <- log(gamma((v+(lag+1)*l)/2)) - log(gamma(v/2)) - 0.5*l*(lag+1)*log(pi*v) # constant part of the log likelihood function.

op <- optim(int,L,hessian=hes,control=list(trace=1,maxit=maxiter,reltol=1e-14),method=meth) 
a <- op$par ; Like <- op$value ; hess <- op$hes

################################### Parameters Estimated ##################
S <- BlockTop(a[1:((lag+3)*(l+1)*l/2)],l,lag)$S
if(c!=0) {MU <- matrix(c(a[(((lag+3)*(l+1)*l/2)+1):(((lag+3)*(l+1)*l/2)+l*c)]),l,c)}

if(c!=0) {M0 <- kronecker(diag(1,lag+1),MU)}
if(c==0) {M0 <- kronecker(diag(1,lag+1),0)}

if(c!=0) {M <- M0%*%t(TREND)} ; if (c==0) {M <- 0}

PP <- Par.dlrm(a,S,MU,l,lag,v,T,M,c,TREND) # Estimation of parameters of conditional distribution
Delta <- PP$Delta ; Delta0 <- PP$Delta0
B1 <- PP$B1 ; s2 <- PP$s2 ; Q <- PP$Q
var.coef <- PP$var.coef

if(Trend[1]!= 0) beta <- cbind(Delta0,t(B1))
if(Trend[1] == 0) beta <- cbind(t(B1))

rownames(beta) <- paste(colnames(Data)[1])  # coefficients are given names

head.x <- matrix(nrow=lag,ncol=ncol(Data)) 
if(lag!=0) for(i in 1:lag){head.x[i,] <- paste(colnames(Data),".",i,sep="")}

if(lag!=0) ifelse(mean(Trend[,1])==0 & ncol(Trend)==1, colnames(beta) <- c(colnames(X),t(head.x)), ifelse(mean(Trend)==1 & ncol(Trend)==1,colnames(beta) <- c("const.",colnames(X),t(head.x)),colnames(beta) <- c("const.",colnames(Trend)[2:ncol(Trend)],colnames(X),t(head.x))) )
if(lag==0) ifelse(mean(Trend[,1])==0 & ncol(Trend)==1, colnames(beta) <- c(colnames(X)), ifelse(mean(Trend)==1 & ncol(Trend)==1,colnames(beta) <- c("const.",colnames(X),t(head.x)),colnames(beta) <- c("const.",colnames(Trend)[2:ncol(Trend)],colnames(X))) )

###################################################################
#############Fitted values/Residuals/Con. Covariance###############
###################################################################
q <- v/(v+l*lag+l-2)
Ct <- vector(length=(nrow(cbind(X))-lag)) 
U <- muy <- vector(length=(nrow(cbind(X))-lag))
trend <- vector(length=ncol(cbind(X))-1)
for(i in 1:(nrow(cbind(X))-lag)) 
{
  if(c!=0) {Ct[i] <- 1 + t(Z[2:ncol(S),i]-M[2:ncol(S)])%*%Q%*%(Z[2:ncol(S),i]-M[2:ncol(S)])}
  if(c==0) {Ct[i] <- 1 + t(Z[2:ncol(S),i])%*%Q%*%(Z[2:ncol(S),i])}
  if(c!=0) {trend[i] <- Delta0%*%(Trend[i,])} ; if(c==0) {trend[i] <- 0%*%(Trend[i,])}
  muy[i] <- trend[i] + t(B1)%*%Z[c(2:l,rep(1:l,lag)),i] 
  U[i] <- t(Z[1,i])- muy[i]   
}

#############################################################
if(hes=="TRUE") VARth <- solve(op$hessian)
############################Coefficients/SEs/P-values########
################
##Jacobian##SE##
################

if(hes=="TRUE") 
{
Jc <- Jacob.dlrm(a,lag,l,v,c)$J #jacobian(J,a) # Computing the Jacobian matrix
SE <- sqrt(diag(Jc%*%VARth%*%t(Jc)))
p_value <- 2*(1-pt(abs(c(Delta0,t(B1),vech(s2)))/SE,(T-lag)))
COEF <- round(cbind(c(Delta0,t(B1),vech(s2)),c(SE[1:length(SE)]),c(p_value[1:length(SE)])),8)
colnames(COEF) <- c("coef.","std.err.","p-value")

Jv <- ConJacob.dlrm(a,lag,l,v,c)$Jv # Computing the Jacobian matrix
SEv <- sqrt(diag(Jv%*%VARth%*%t(Jv)))
p_valuev <- 2*(1-pt(abs(c(var.coef))/SEv,(T-lag)))
VAR.COEF <- round(cbind(c(var.coef),c(SEv[1:length(SEv)]),c(p_valuev[1:length(SEv)])),8)
colnames(VAR.COEF) <- c("var.coef","std.err.","p-value")
B <- COEF[1:ncol(beta),] ; rownames(B) <- colnames(beta)
rownames(B) <- colnames(beta)
}

if(hes=="FALSE") { 
COEF <- round(cbind(c(Delta0,t(B1),vech(s2))),8) ; colnames(COEF) <- c("coef.")
VAR.COEF <- round(cbind(var.coef),8) ; colnames(VAR.COEF) <- c("var.coef")

}

head.var.coef <- matrix(nrow=(ncol(X)+1),ncol=lag) 
if(lag!=0) for(i in 1:lag){head.var.coef[,i] <- c(paste(colnames(Data),".",i,sep=""))}

if(lag!=0) coef.label <- c(colnames(X),head.var.coef)
if(lag==0) coef.label <- c(colnames(X))

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

ad <- Student(res,Data,Ct,s2,v,lag,B1)$ad

p <- ifelse(hes=="FALSE",ncol(beta),ncol(beta))

fitted <- muy
R2 <- sum((fitted-mean(Data[((lag+1):nrow(Data))]))^2 )/sum((Data[((lag+1):nrow(Data))]-mean(Data[((lag+1):nrow(Data))]))^2)
F <- (R2/(1-R2))*(nrow(Data)-p)/(p-1)

if(hes=="TRUE") result <- list(beta=B,var.coef=VAR.COEF,like=-Like,sigma=cbind("s^2"=s2, "stderr.sigma"=tail(COEF[,2],length(vech(s2)))),cvar=Ct,trend=trend,res=res,fitted=muy,init=a,hes=hess,S=S,R.squared=R2,F.stat=F,ad=ad)
if(hes=="FALSE") result <- list(beta=beta,var.coef=VAR.COEF,like=-Like,sigma=s2,cvar=Ct,trend=trend,res=res,fitted=muy,init=a,S=S,R.squared=R2,F.stat=F,ad=ad)

return(result)
}