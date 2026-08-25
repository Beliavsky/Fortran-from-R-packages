uni.selection=function(t.vec, d.vec, X.mat, P.value=0.001,K=10,score=TRUE,d0=0,
                       randomize=FALSE,CC.plot=FALSE,permutation=FALSE,M=200){
  
  n=length(t.vec)
  if(randomize==TRUE){
    rand=sample(1:n)  ### randomize patient IDs before cross-validation ###
    t.vec=t.vec[rand]
    d.vec=d.vec[rand]
    X.mat=X.mat[rand,]
  } 
  
  p=ncol(X.mat)
  if(score==TRUE){ res=uni.score(t.vec, d.vec, X.mat,d0) }else{ 
    res=uni.Wald(t.vec, d.vec, X.mat) 
  }
  temp=res$P<P.value
  if(sum(temp)==0) stop("no feature is selected; please increase P.value")

    Beta=res$beta[temp]
    Z=res$Z[temp]
    P=res$P[temp]
    X.cut=as.matrix(X.mat[,temp])
    q=ncol(X.cut)

    if(score==TRUE){ w=Z }else{ w=Beta }
    CC=X.cut%*%w
    c.index0=1-concordance(  Surv(t.vec,d.vec)~CC  )$concordance
    
    ### Log-partial likelihood ###
    atr_t=(matrix(t.vec,n,n,byrow=TRUE)>=matrix(t.vec,n,n))
    l.func=function(g){
      l=sum( (d.vec*CC*g) )-sum( d.vec*log(atr_t%*%exp(CC*g)) )
      as.numeric( l )
    }
   
    ### Cross-validation ###
    CC.CV=CC.test=NULL
    CVL=RCVL1=RCVL2=0
    g_kk_vec=numeric(K)
    folds=cut(seq(1,n),breaks=K,labels=FALSE)
    
    for(k in 1:K){
      fold_k=which(folds==k)
      t_k=t.vec[-fold_k]
      d_k=d.vec[-fold_k]
      CC_k=CC[-fold_k]
      X_k=X.mat[-fold_k,]
      n_k=length(t_k)
      
      if(score==TRUE){ res=uni.score(t_k, d_k, X_k,d0) }else{ 
        res=uni.Wald(t_k, d_k, X_k)
      }
      temp_k=res$P<P.value
      
      if(sum(temp_k)==0){
        CC_kk=rep(0,n)
        warning("no feature is selected in a cross-validation fold; please increase P.value")
      }else{
        if(score==TRUE){ w_k=res$Z[temp_k] }else{ w_k=res$beta[temp_k] }
        CC_kk=as.matrix(X.mat[,temp_k])%*%w_k
      }
      CC.test=c(CC.test,CC_kk[fold_k])

      ##### Resubstitution CVL (without cross-validation) #####
      res_k=coxph(Surv(t_k,d_k)~CC_k)
      RCVL1=RCVL1+l.func(res_k$coef)-res_k$loglik[2] 

      ##### Resubstitution CVL (only coefficients are cross-validated) ### 
      if(score==TRUE){ w_CV=uni.score(t_k, d_k, X.cut[-fold_k,],d0)$Z }else{ 
        w_CV=uni.Wald(t_k, d_k, X.cut[-fold_k,])$beta
      }
      CC.CV=c(CC.CV,X.cut[fold_k,]%*%as.matrix(w_CV))
      CC.CV_k=X.cut[-fold_k,]%*%as.matrix(w_CV)
      res_CV_k=coxph(Surv(t_k,d_k)~CC.CV_k)
      RCVL2=RCVL2+l.func(res_CV_k$coef)-res_CV_k$loglik[2] 
     
      ### CVL (fullly cross-validated) ###
      res_kk=coxph(Surv(t_k,d_k)~CC_kk[-fold_k])
      l_kk.func=function(g){
        l=sum( (CC_kk*g)[d.vec==1] )-sum( (log(atr_t%*%exp(CC_kk*g)))[d.vec==1] )
        as.numeric( l )
      }
      CVL=CVL+l_kk.func(res_kk$coef)-as.numeric(res_kk$loglik[2])
    }
    
    ##### c-index (Cross-validating estimation of CC) 
    c.index1=1-concordance(  Surv(t.vec,d.vec)~CC.CV  )$concordance
    c.index2=1-concordance(  Surv(t.vec,d.vec)~CC.test  )$concordance
    c.index=c("No cross-validation"=c.index0,"Incomplete cross-validation"=c.index1,
              "Full cross-validation"=c.index2)
    
    if(CC.plot==TRUE){
    par(mfrow=c(1,2))
    plot(CC,CC.test,xlab="CC (Not cross-validated)",
         ylab="CC (Fully cross-validated)" )
    plot(CC.CV,CC.test,xlab="CC (Only estimation cross-validated)",
         ylab="CC (Fully cross-validated)" )
    }
    
    FDR.perm=f.perm=NULL
    if(permutation==TRUE){
      q.perm=numeric(M)
      for(i in 1:M){
        set.seed(i)
        if(score==TRUE){ res=uni.score(t.vec, d.vec, X.mat[sample(1:n),],d0) }else{ 
          res=uni.Wald(t.vec, d.vec, X.mat[sample(1:n),])
        }
        q.perm[i]=sum(res$P<P.value)
      }
      f.perm=mean(q.perm) ## No. of false discoveries
      FDR.perm=f.perm/q ## False discovery rate
    }
    
    list(beta=Beta[order(P)],Z=Z[order(P)],P=P[order(P)],
         CVL=c("CVL"=CVL,"RCVL1"=RCVL1,"RCVL2"=RCVL2),
         Genes=c("No. of genes"=p,"No. of selected genes"=q,
                 "No. of falsely selected genes"=f.perm),
         FDR=c("P.value * (No. of genes)"=P.value*p/q,"Permutation"=FDR.perm)
    )

}
