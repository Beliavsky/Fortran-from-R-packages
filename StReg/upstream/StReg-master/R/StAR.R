StAR <-
function(Data,Trend=1,lag=1,v=1,maxiter=1000,meth="BFGS",hes="FALSE",init="na")
{
Data <- cbind(Data) # Data matrix
if (ncol(Data)!=1) stop("data must be in a vector form") # Error Handling: Error message if the data is not in vector form.
if (lag-round(lag)!=0) stop("lag must be an integer") # Error Handling: Lag must be a positive integer
if (lag<=0) stop("lag number must be positive integer") # Error Handling: lag must be a positive integer

Z <- t(embed(Data,lag+1)) ; T <- nrow(Data) ; l <- 1  # Sample size and l=1 indicating a univariate data matrix

if (Trend[1]==1 & length(Trend)==1) {Trend <- matrix(1,T,1) ; colnames(Trend) <- "const"; TREND <- embed(Trend,lag+1) ; c <- 1} # If 1 is then matrix of 1s is created for intercept term.

if (Trend[1]==0 & length(Trend)==1) {Trend <- matrix(0,T,1) ; TREND <- embed(Trend,lag+1) ; c <- 0} # If 0, then intercept is excluded from the model.

if(ncol(cbind(Trend))>1) {TREND <- embed(Trend,lag+1) ; c <- ncol(Trend)} # If more than one column is supplied in Trend, then use the whole matrix as mean (t varying intercept).
if(Trend[1]!=1) {if (is.null(colnames(Trend))=="TRUE") {colnames(Trend) <- paste("t",".",seq(1,ncol(Trend),1),sep="")}} # If Trend matrix does not have column names provided, then automatically created.

# Error Handling: Error message if Trend matrix and Data matrix do not have same number of rows. 
if (Trend[1]!=1) {if (nrow(Data)!= nrow(Trend)) stop("Data and Trend matrix must have same number of rows")}

if (ncol(Data)>2) stop("number of variables must equal to 1") # Error Handling: Number of columns in the data matrix must be 1.
if (Trend[1]!=1) {if (ncol(Trend)<1) stop("The Trend matrix must have at least one column")} # Error Handling: The trend matrix must have at least one column


# Error Handling: Error messaage if degrees of freedom is not a positive number
if (v<=0) stop("degrees of freedom must be a positive number")

# Error Handling: Error messaage if maximum iterations for numerical optimization is less than 10
if (maxiter<10) stop("Iteration must be at least 10")

# Error Handling: Error messaage if number of parameters is too large. Reduce the variables or increase the sample size.
if (2*T<3*((lag+3)*(l+1)*l/2)+3*l+c) stop(list(((lag+3)*(l+1)*l/2)+3*l+c,"Too many parameters for given sample size. Reduce the number of lags."))

##Likelihood Function
L<-function(a)
{
S <- BlockTop(a[1:((lag+3)*(l+1)*l/2)],l,lag)$S ##Var-Cov Matrix created using block Toipletz matrix from BlockTop function
F <- solve(S)/v ## VarCov and its Inverse

if(c!=0) {MU <- matrix(c(a[(((lag+3)*(l+1)*l/2)+1):(((lag+3)*(l+1)*l/2)+l*c)]),l,c)} # Mean if intercept is to be included in the model

if(c!=0) {M0 <- kronecker(diag(1,lag+1),MU)} ## Kronecker product of the Mean.
if(c==0) {M0 <- kronecker(diag(1,lag+1),0)} ## Kronecker product of the Mean.

if(c!=0) {M <- M0%*%t(TREND)} ; if (c==0) {M <- 0}

D <- 1 + diag(t(Z-M)%*%F%*%(Z-M)) # Quadratic form of the Student's t likelihood function

## Likelihood function
 LLn <- (T-lag)*const - 0.5*(T-lag)*log(det(S)) - 0.5*(v+(lag+1)*l)*sum(log(D))
 neg.LLn <- -LLn # negative of likelihood function to be minimized.
 neg.LLn

}

if(init[1]=="na") int <- c(runif(((lag+3)*(l+1)*l/2)+l*c,0,.1))   # Initialization of the parameters if not provided

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

PP <- Par.star(a,S,MU,l,lag,v,T,M,c,TREND) # Estimation of parameters of conditional distribution
Delta <- PP$Delta ; Delta0 <- PP$Delta0
B1 <- PP$B1 ; s2 <- PP$s2 ; Q <- PP$Q
var.coef <- PP$var.coef
if(mean(Trend[,1]) ==0) beta <- cbind( t(B1))
if(mean(Trend[,1]) !=0) beta <- cbind(Delta0, t(B1))

rownames(beta) <- paste(colnames(Data)) # coefficients are given names
head.x <- matrix(nrow=lag,ncol=ncol(Data)) ; for(i in 1:lag){head.x[i,] <- paste(colnames(Data),".",i,sep="")}
ifelse(mean(Trend[,1])==0 & ncol(Trend)==1, colnames(beta) <- c(t(head.x)), ifelse(mean(Trend[,1])==1 & ncol(Trend)==1,colnames(beta) <- c("const.",t(head.x)),colnames(beta) <- c("const.",colnames(Trend)[2:ncol(Trend)],t(head.x))))

#############Fitted values/Residuals/Con. Covariance###############
q <- v/(v+l*lag-2)
Ct <- vector(length=(T-lag)) 
U <- muy <- matrix(nrow=l,ncol=(T-lag))
trend <- vector(length=ncol(cbind(Data))-1)

for(i in 1:(T-lag)) 
{
  if(c!=0) {Ct[i] <- 1 + t(Z[(l+1):((lag+1)*l),i]-M[(l+1):((lag+1)*l)])%*%Q%*%(Z[(l+1):((lag+1)*l),i]-M[(l+1):((lag+1)*l)])}
  if(c==0) {Ct[i] <- 1 + t(Z[(l+1):((lag+1)*l),i])%*%Q%*%(Z[(l+1):((lag+1)*l),i])}
  if(c!=0) {trend[i] <- Delta0%*%(Trend[i,])} ; if(c==0) {trend[i] <- 0%*%(Trend[i,])}
  muy[,i] <- trend[i] + t(B1)%*%Z[(l+1):(l*(lag+1)),i] 
  U[,i] <- t(Z[1:l,i])- muy[,i]   
}

res <- t(U)
#############################################################
if(hes=="TRUE") VARth <- solve(op$hessian)
############################Coefficients/SEs/P-values########
################
##Jacobian##SE##
################

if(hes=="TRUE") 
{
Jc <- Jacob.star(a,lag,l,v,c)$J 	# Computing the Jacobian matrix
SE <- sqrt(diag(Jc%*%VARth%*%t(Jc)))
p_value <- 2*(1-pt(abs(c(Delta0,t(B1),vech(s2)))/SE,(T-lag)))
COEF <- round(cbind(c(Delta0,t(B1),vech(s2)),c(SE[1:length(SE)]),c(p_value[1:length(SE)])),8)
colnames(COEF) <- c("coef.","std.err.","p-value")

Jv <- ConJacob.star(a,lag,l,v,c)$Jv # Computing the Jacobian matrix
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
B <- beta
}


head.var.coef <- matrix(nrow=l,ncol=lag) 
for(i in 1:lag)
{
head.var.coef[,i] <- paste(colnames(Data),".",i,sep="")
}

var.head <-matrix(nrow=lag*l,ncol=lag*l)
for(i in 1:(lag*l))
{
for(j in 1:(lag*l))
{
var.head[i,j] <- paste(c(head.var.coef)[j],".",c(head.var.coef)[i],sep="")
}
}

rownames(VAR.COEF) <- c("const.",vech(var.head))

if(Trend[1]!= 0) beta <- cbind(Delta0,t(B1))
if(Trend[1] == 0) beta <- cbind(t(B1))

fitted <- t(muy)
res <- cbind(U)

ad <- Student(res,Data,Ct,s2,v,lag,B1)$ad

p <- ifelse(hes=="FALSE",ncol(beta),ncol(beta))

R2 <- sum((fitted-mean(Data[((lag+1):nrow(Data))]))^2 )/sum((Data[((lag+1):nrow(Data))]-mean(Data[((lag+1):nrow(Data))]))^2)
F <- (R2/(1-R2))*(nrow(Data)-p)/(p-1) 

if(hes=="TRUE") result <- list(beta=B,var.coef=VAR.COEF,like=-Like,sigma=cbind("s^2"=s2, "stderr.sigma"=tail(COEF[,2],length(vech(s2)))),cvar=Ct,trend=trend,res=res,fitted=fitted,init=a,hes=hess,S=S,R.squared=R2,F.stat=F,ad=ad)
if(hes=="FALSE") result <- list(beta=B,var.coef=VAR.COEF,like=-Like,sigma=s2,cvar=Ct,trend=trend,res=res,fitted=fitted,init=a,S=S,R.squared=R2,F.stat=F,ad=ad)
return(result)
}
