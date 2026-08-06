StVAR <-
function(Data,Trend=1,lag=1,v=1,maxiter=1000,meth="BFGS",hes="FALSE",init="na")
{

if (lag-round(lag)!=0) stop("lag must be an integer")
if (lag<=0) stop("lag number must be positive integer")

Z <- t(embed(Data,lag+1)) ; T <- nrow(Data) ; l <- ncol(Data)  ## No. of variables in VAR
if (mean(Trend)==1 & length(Trend)==1) {Trend <- matrix(1,T,1); TREND <- embed(Trend,lag+1) ; c <- 1}

if (mean(Trend)==0 & length(Trend)==1) {Trend <- matrix(0,T,1); TREND <- embed(Trend,lag+1) ; c <- 0}

if(ncol(cbind(Trend))>1) {TREND <- embed(Trend,lag+1) ; c <- ncol(Trend)}
if(Trend[1]!=1) {if (is.null(colnames(Trend))=="TRUE") {colnames(Trend) <- paste("t",".",seq(1,ncol(Trend),1),sep="")}}

if (Trend[1]!=1) {if (nrow(cbind(Data))!= nrow((Trend))) stop("Data and Trend matrix must have same number of rows")}
if (ncol(cbind(Data))<1) stop("no. of variables must be at least 1")
if (Trend[1]!=1) {if (ncol(Trend)<1) stop("The Trend matrix must have at least one column")}
if (v<=0) stop("degrees of freedom must be a positive integer")
if (maxiter<10) stop("Iteration must be at least 10")
if (2*T<3*((lag+3)*(l+1)*l/2)+3*l+c) stop(list(((lag+3)*(l+1)*l/2)+3*l+c,"Too many parameters for given sample size. Reduce the number of lags."))


##Likelihood Function
L<-function(a)
{
 S <- BlockTop(a[1:((lag+3)*(l+1)*l/2)],l,lag)$S ##Var-Cov Matrix
 F <- solve(S)/v ## VarCov and its Inverse

if(c!=0) {MU <- matrix(c(a[(((lag+3)*(l+1)*l/2)+1):(((lag+3)*(l+1)*l/2)+l*c)]),l,c)}

if(c!=0) {M0 <- kronecker(diag(1,lag+1),MU)}
if(c==0) {M0 <- kronecker(diag(1,lag+1),0)}

if(c!=0) {M <- M0%*%t(TREND)} ; if (c==0) {M <- 0}

D <- 1 + diag(t(Z-M)%*%F%*%(Z-M)) ##Quadratic form

## Likelihood function
 LLn <- (T-lag)*const - 0.5*(T-lag)*log(det(S)) - 0.5*(v+(lag+1)*l)*sum(log(D))
 neg.LLn <- -LLn
 neg.LLn
}

if(init[1]=="na") int <- c(runif(((lag+3)*(l+1)*l/2)+l*c,0,.1))  ##Initialization

if(init[1]!="na") int<-init
const <- log(gamma((v+(lag+1)*l)/2)) - log(gamma(v/2)) - 0.5*l*(lag+1)*log(pi*v)

## Maximization of the likelihood function
op <- optim(int,L,hessian=hes,control=list(trace=1,maxit=maxiter,reltol=1e-14),method=meth) 
a <- op$par ; Like <- op$value ; hess <- op$hes

###################################Parameters##################
###################################Parameters##################
###################################Parameters##################
S <- BlockTop(a[1:((lag+3)*(l+1)*l/2)],l,lag)$S

if(c!=0) {MU <- matrix(c(a[(((lag+3)*(l+1)*l/2)+1):(((lag+3)*(l+1)*l/2)+l*c)]),l,c)}

if(c!=0) {M0 <- kronecker(diag(1,lag+1),MU)}
if(c==0) {M0 <- kronecker(diag(1,lag+1),0)}

if(c!=0) {M <- M0%*%t(TREND)} ; if (c==0) {M <- 0}

PP <- Par.stvar(a,S,MU,l,lag,v,T,M,c)
Delta <- PP$Delta ; Delta0 <- PP$Delta0
B1 <- PP$B1 ; s2 <- PP$s2 ; Q <- PP$Q
var.coef <- PP$var.coef

beta <- cbind(Delta0,t(B1))

rownames(beta) <- paste(colnames(Data))
head.x <- matrix(nrow=lag,ncol=ncol(Data)) ; for(i in 1:lag){head.x[i,] <- paste(colnames(Data),".",i,sep="")}
ifelse(mean(Trend)!=1, colnames(beta) <- c(colnames(Trend),t(head.x)), colnames(beta) <- c("const.",t(head.x)) )

###################################################################
#############Fitted values/Residuals/Con. Covariance###############
###################################################################
q <- v/(v+l*lag-2)
Ct <- Ctt <- vector(length=(T-lag)) 
U <- muy <- matrix(nrow=l,ncol=(T-lag))
for(i in 1:(T-lag)) 
{
  if(c!=0) {Ct[i] <- 1 + t(Z[(l+1):((lag+1)*l),i]-M[(l+1):((lag+1)*l)])%*%Q%*%(Z[(l+1):((lag+1)*l),i]-M[(l+1):((lag+1)*l)])}
  if(c==0) Ct[i] <- 1 + t(Z[(l+1):((lag+1)*l),i])%*%Q%*%(Z[(l+1):((lag+1)*l),i])
  if(c!=0) {muy[,i] <- Delta[i,] + t(B1)%*%Z[(l+1):(l*(lag+1)),i] }
  if(c==0) {muy[,i] <- Delta + t(B1)%*%Z[(l+1):(l*(lag+1)),i] }
  U[,i] <- t(Z[1:l,i])- muy[,i]   
}

#############################################################
if(hes=="TRUE") VARth <- solve(op$hessian)
############################Coefficients/SEs/P-values########
################
##Jacobian##SE##
################

if(hes=="TRUE") 
{Jc <- Jacob.stvar(a,lag,l,v,c)$J #jacobian(J,a)
SE <- sqrt(diag(Jc%*%VARth%*%t(Jc)))
p_value <- 2*(1-pt(abs(c(Delta0,t(B1),vech(s2)))/SE,(T-lag)))
COEF <- round(cbind(c(Delta0,t(B1),vech(s2)),c(SE[1:length(SE)]),c(Delta0,B1,vech(s2))/c(SE[1:length(SE)]),c(p_value[1:length(SE)])),8)
colnames(COEF) <- c("coef.","std.err.","t-value","p-value")

## Labeling the beta matrices
k <- ncol(beta)
B <- array(dim=c(k,4,l))
std.err <- matrix(COEF[1:(k*l),2],l,k)
t.value <- matrix(COEF[1:(k*l),3],l,k)
p.value <- matrix(COEF[1:(k*l),4],l,k)

for(i in 1:l)
{
B[,,i] <- cbind(beta[i,],std.err[i,],t.value[i,],p.value[i,]) 
}

dimnames(B)[[1]] <- colnames(beta)
dimnames(B)[[2]] <- colnames(COEF)
dimnames(B)[[3]] <- colnames(Data)


Jv <- ConJacob.stvar(a,lag,l,v,c)$Jv #jacobian(J,a)
SEv <- sqrt(diag(Jv%*%VARth%*%t(Jv)))
p_valuev <- 2*(1-pt(abs(c(var.coef))/SEv,(T-lag)))
VAR.COEF <- round(cbind(c(var.coef),c(SEv[1:length(SEv)]),c(var.coef)/c(SEv[1:length(SEv)]),c(p_valuev[1:length(SEv)])),8)
colnames(VAR.COEF) <- c("var.coef","std.err.","t-value","p-value")
}

if(hes=="FALSE") { 
COEF <- round(cbind(c(Delta0,t(B1),vech(s2))),8) ; colnames(COEF) <- c("coef.")
VAR.COEF <- round(cbind(var.coef),8) ; colnames(VAR.COEF) <- c("var.coef")
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


## Trends
trends <- Delta
if(c!=0) {head.tre <- paste(colnames(Data),".","tre",sep="") ; colnames(trends) <- head.tre}

## Residuals
res <- t(U)
head.res <- paste(rep("u",ncol(Data)),".",colnames(Data),sep="")
colnames(res) <- head.res

## Fitted Values
fitted <- t(muy)
head.fit <- paste(colnames(Data),".",rep("hat",ncol(Data)),sep="")
colnames(fitted) <- head.fit

ifelse(c!=0,COEF <- COEF, COEF <- COEF[2:nrow(COEF),])
ifelse(c!=0,beta <- beta, beta <- beta[,2:ncol(beta)])


ad <- Student(res,Data,Ct,s2,v,lag,B1)$ad

p <- ifelse(hes=="FALSE",ncol(beta),ncol(beta))

F <- R2 <- vector(length=ncol(Data))
for(k in 1:ncol(Data))
{
R2[k] <- sum((fitted[,k]-mean(Data[((lag+1):nrow(Data)),k]))^2 )/sum((Data[((lag+1):nrow(Data)),k]-mean(Data[((lag+1):nrow(Data)),k]))^2)
F[k] <- (R2[k]/(1-R2[k]))*(nrow(Data)-p)/(p-1)
}

if(hes=="TRUE") result <- list(beta=B,var.coef=VAR.COEF,like=-Like,sigma=cbind("s^2"=vech(s2), "stderr.sigma"=tail(COEF[,2],length(vech(s2)))),cvar=Ct,trends=trends,res=res,fitted=fitted,init=a,hes=hess,S=S,R.squared=R2,F.stat=F,ad=ad)
if(hes=="FALSE") result <- list(beta=beta,var.coef=VAR.COEF,like=-Like,sigma=s2,cvar=Ct,trends=trends,res=res,fitted=fitted,init=a,S=S,R.squared=R2,F.stat=F,ad=ad)

return(result)
}
