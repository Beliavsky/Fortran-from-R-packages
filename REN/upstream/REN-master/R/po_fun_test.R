utils::globalVariables(c("i", "m", "Method", "value", "time", "var", "sd"))

#' Insert Values at Specified Positions in a Vector
#'
#' This function inserts specified values at given positions in a vector.
#'
#' @param a A vector.
#' @param pos A numeric vector specifying the positions to insert the values.
#' @param ... Values to be inserted at the specified positions.
#' @return A vector with the inserted values.
#' @export
insert.at <- function(a, pos, ...){
  dots <- list(...)
  stopifnot(length(dots)==length(pos))
  result <- vector("list",2*length(pos)+1)
  result[c(TRUE,FALSE)] <- split(a, cumsum(seq_along(a) %in% (pos+1)))
  result[c(FALSE,TRUE)] <- dots
  result <- unlist(result)
  if (pos == 0){
    result <- unlist(c(dots,a))
  }
  if (pos == length(a)){
    result <- unlist(c(a,dots))
  }
  return(result)
}

#' Perform LASSO or Ridge Regression for Portfolio Optimization
#'
#' This function performs LASSO, Ridge, or Elastic Net regression for portfolio optimization.
#'
#' @param y0 A numeric vector of response values.
#' @param x0 A numeric matrix of predictors.
#' @param method The regularization method: "LASSO", "RIDGE", or "EN" (Elastic Net).
#' @return A numeric vector of optimized portfolio weights.
#' @export
po.avg <- function (y0, x0, method = "LASSO"){
  p0 <- dim(x0)[2]
  j <- which(apply(x0, 2, var) <= 0, arr.ind = TRUE)
  x0 <- x0[,which(apply(x0, 2, var) > 0, arr.ind = TRUE)]
  p <- dim(x0)[2]
  w <- foreach (i=1:p, .combine = rbind, .packages='glmnet', .export='insert.at') %dopar% {
    y_tmp <- x0[,i] - y0
    x_tmp <- x0[,i] - x0[,-i]
    switch(method,
           LASSO = {cv <- cv.glmnet(x_tmp,y_tmp,parallel = TRUE)
           w_tmp <- cv$glmnet.fit$beta[,which(cv$lambda == cv$lambda.min)]
           insert.at(w_tmp,i-1,1-sum(as.numeric(w_tmp)))},

           RIDGE = {cv <- cv.glmnet(x_tmp,y_tmp,alpha = 0,parallel = TRUE)
           w_tmp <- cv$glmnet.fit$beta[,which(cv$lambda == cv$lambda.min)]
           insert.at(w_tmp,i-1,1-sum(as.numeric(w_tmp)))},

           EN = {a <- seq(0.1, 0.9, 0.1)
           search <- foreach(m = a, .combine = rbind) %dopar% {
             cv <- cv.glmnet(x_tmp,y_tmp,alpha = m)
             data.frame(cvm = cv$cvm[cv$lambda == cv$lambda.min], lambda.1se = cv$lambda.min, alpha = m)}
           cv <- search[search$cvm == min(search$cvm), ]
           w_tmp <- as.matrix(glmnet(x_tmp,y_tmp, lambda = cv$lambda.min, alpha = cv$alpha)$beta)
           insert.at(w_tmp,i-1,1-sum(as.numeric(w_tmp)))})
  }
  w <- colMeans(w)
  if (length(w) != p0) {
    w <- insert.at(w,j-1,0)
  }
  return(w)}

#' Perform Gross Exposure Portfolio Optimization
#'
#' This function performs gross exposure portfolio optimization using LASSO.
#'
#' @param y0 A numeric vector of response values.
#' @param x0 A numeric matrix of predictors.
#' @param method The regularization method: "NOSHORT" or "EQUAL".
#' @return A numeric vector of optimized portfolio weights.
#' @export
po.grossExp <- function (y0, x0, method = "NOSHORT"){
  p0 <- dim(x0)[2]
  j <- which(apply(x0, 2, var) <= 0, arr.ind = TRUE)
  x0 <- x0[,which(apply(x0, 2, var) > 0, arr.ind = TRUE)]
  p <- dim(x0)[2]
  switch(method,
         NOSHORT = {Dmat <- matrix(nearPD(cov(x0))["mat"][[1]]@x,p,p)   # cov(x0)
         dvec <- as.matrix(rep(0, p))
         A.ge <- diag(p)
         b.ge <- as.matrix(rep(0, p))
         A.eq <- as.matrix(rep(1, p))
         b.eq <- 1
         sol  <- solve.QP(Dmat,dvec,Amat=cbind(A.ge,A.eq),bvec=c(b.ge,b.eq))$solution
         if (any(round(sol,5) == 1)){
           k <- which(round(sol,5)==1)
           y_tmp <- x0[,k] - y0
           x_tmp <- x0[,k] - x0[,-k]
           cv <- cv.glmnet(x_tmp,y_tmp,parallel = TRUE)
           w <- cv$glmnet.fit$beta[,which(cv$lambda == cv$lambda.min)]
           w <- insert.at(w,k-1,1-sum(as.numeric(w)))
         } else {
           y_tmp <- x0%*%sol - y0
           x_tmp <- as.numeric(x0%*%sol) - x0
           cv <- cv.glmnet(x_tmp,y_tmp,parallel = TRUE)
           w <- cv$glmnet.fit$beta[,which(cv$lambda == cv$lambda.min)]
           w <- (1-sum(cv$glmnet.fit$beta[,which(cv$lambda == cv$lambda.min)]))*sol + cv$glmnet.fit$beta[,which(cv$lambda == cv$lambda.min)]
         }},

         EQUAL = {y_tmp <- rowMeans(x0) - y0
         x_tmp <- rowMeans(x0) - x0
         cv <- cv.glmnet(x_tmp,y_tmp)
         w <- (1-sum(cv$glmnet.fit$beta[,which(cv$lambda == cv$lambda.min)]))*rep(1/p,p) + cv$glmnet.fit$beta[,which(cv$lambda == cv$lambda.min)]})
  if (length(w) != p0) {
    w <- insert.at(w,j-1,0)
  }
  return(w)
}

#' Perform Shrinkage-Based Portfolio Optimization
#'
#' This function uses covariance shrinkage techniques for portfolio optimization.
#'
#' @param y0 A numeric vector of response values.
#' @param x0 A numeric matrix of predictors.
#' @return A numeric vector of optimized portfolio weights.
#' @export
po.covShrink <- function (y0, x0){
  p0 <- dim(x0)[2]
  j <- which(apply(x0, 2, var) <= 0, arr.ind = TRUE)
  x0 <- x0[,which(apply(x0, 2, var) > 0, arr.ind = TRUE)]
  p <- dim(x0)[2]
  Dmat <- matrix(nearPD(cov.shrink(cov(x0)))["mat"][[1]]@x,p,p) #cov(x0)
  dvec <- as.matrix(rep(0, p))
  A.eq <- as.matrix(rep(1, p))
  b.eq <- 1
  w  <- solve.QP(Dmat,dvec,Amat=A.eq,bvec=b.eq)$solution
  if (length(w) != p0) {
    w <- insert.at(w,j-1,0)
  }
  return(w)
}

#' Perform Simple Linear Model for Portfolio Optimization
#'
#' This function uses simple linear regression to perform portfolio optimization.
#'
#' @param y0 A numeric vector of response values.
#' @param x0 A numeric matrix of predictors.
#' @return A numeric vector of optimized portfolio weights.
#' @export
po.cols <- function (y0, x0){
  p0 <- dim(x0)[2]
  j <- which(apply(x0, 2, var) <= 0, arr.ind = TRUE)
  x0 <- x0[,which(apply(x0, 2, var) > 0, arr.ind = TRUE)]
  y_tmp <- x0[,1] - y0
  x_tmp <- x0[,1] - x0[,-1]
  w_tmp <- lm(y_tmp~x_tmp -1)$coefficients
  w <- insert.at(w_tmp,0,1-sum(as.numeric(w_tmp)))
  if (length(w) != p0) {
    w <- insert.at(w,j-1,0)
  }
  return(w)
}

#' Perform James-Stein Portfolio Optimization
#'
#' This function uses James-Stein estimation for portfolio optimization.
#'
#' @param x0 A numeric matrix of asset returns.
#' @return A numeric vector of optimized portfolio weights.
#' @export
po.JM <- function (x0){
  p0 <- dim(x0)[2]
  j <- which(apply(x0, 2, var) <= 0, arr.ind = TRUE)
  x0 <- x0[,which(apply(x0, 2, var) > 0, arr.ind = TRUE)]
  p <- dim(x0)[2]
  Dmat <- matrix(nearPD(cov(x0))["mat"][[1]]@x,p,p)  #  cov(x0)
  dvec <- as.matrix(rep(0, p))
  A.ge <- diag(p)
  b.ge <- as.matrix(rep(0, p))
  A.eq <- as.matrix(rep(1, p))
  b.eq <- 1
  w  <- solve.QP(Dmat,dvec,Amat=cbind(A.ge,A.eq),bvec=c(b.ge,b.eq))$solution
  if (length(w) != p0) {
    w <- insert.at(w,j-1,0)
  }
  return(w)
}

#' Perform Hierarchical Clustering on Asset Correlations
#'
#' This function performs hierarchical clustering on asset correlations and returns the clustered groups.
#'
#' @param x A numeric matrix of asset returns.
#' @return A list of asset clusters.
#' @export
buh.clust <- function (x){
  p0 <- dim(x)[2]
  j <- which(apply(x, 2, var) <= 0, arr.ind = TRUE)
  x <- x[,which(apply(x, 2, var) > 0, arr.ind = TRUE)]
  p <- dim(x)[2]
  CorMatix <- cor(x)-diag(p)
  index <- seq(1,p,1)
  group <- list()
  b <- 0
  for (i in (1:p)) {
    group[[i]] <- index[i]
  }
  rho <- rep(NA,p-1)
  tau <- list()
  while (b <= (p-2)) {
    b <- b+1
    rho.max <- max(CorMatix)
    rho[b] <- rho.max
    tau[[b]] <- group
    new.clust <- unique(as.numeric(which(CorMatix==max(CorMatix), arr.ind = T)))
    group <- c(0,group)
    group[[1]] <- c(group[[min(new.clust)+1]],group[[max(new.clust)+1]])
    group[[min(new.clust)+1]] <- NULL
    group[[max(new.clust)]] <- NULL
    tmp <- rep(0,length(group))
    for (i in (2:length(group))) {
      if(length(group)<2) next
      tmp[i] <- max(cancor(x[,group[[1]]],x[,group[[i]]])$`cor`)
    }
    CorMatix <- cbind(tmp, rbind(tmp[-1], CorMatix[-new.clust,-new.clust]))
  }
  group <- tau[[which.min(rho)]]
  return(group)
}

# po.bhu <- function (y0, x0, group, rep){
#   p0 <- dim(x0)[2]
#   j <- which(apply(x0, 2, var) <= 0, arr.ind = TRUE)
#   x0 <- x0[,which(apply(x0, 2, var) > 0, arr.ind = TRUE)]
#   p <- dim(x0)[2]
#   g <- length(group)
#   x.bar <- foreach (i=1:g, .combine = cbind) %dopar% {
#     if (length(group[[i]])>1){
#       rowMeans(x0[,group[[i]]])
#     } else {
#       x0[,group[[i]]]
#     }
#   }
#   Dmat <- matrix(nearPD(cov(x.bar))["mat"][[1]]@x,g,g) #cov(x0)
#   dvec <- as.matrix(rep(0, g))
#   A.ge <- diag(g)
#   b.ge <- as.matrix(rep(0, g))
#   A.eq <- as.matrix(rep(1, g))
#   b.eq <- 1
#   sol  <- round(solve.QP(Dmat,dvec,Amat=cbind(A.ge,A.eq),bvec=c(b.ge,b.eq))$solution,10)
#   y_tmp <- x.bar%*%sol - y0
#   x_tmp <- as.numeric(x.bar%*%sol) - x.bar
#   lambda <- exp(seq(log(0.001), log(10), length.out=100))
#   cv <- cv.glmnet(x_tmp,y_tmp,lambda=lambda)
#   w.bar <- (1-sum(cv$glmnet.fit$beta[,which(cv$lambda == cv$lambda.min)]))*sol + cv$glmnet.fit$beta[,which(cv$lambda == cv$lambda.min)]
#   p.bar <- length(which(w.bar!=0,arr.ind = T))
#   w <- foreach (i=1:rep, .combine = rbind, .packages=c('glmnet','quadprog'), .export='insert.at') %dopar% {
#         index <- which(w.bar!=0,arr.ind = T)
#         for (k in (1:p.bar)){
#             index[k] <- group[[index[k]]][sample(1:length(group[[index[k]]]), 1)]
#         }
#         Dmat <- cov(x0[,index])
#         dvec <- as.matrix(rep(0, p.bar))
#         A.eq <- as.matrix(rep(1, p.bar))
#         b.eq <- 1
#         sol  <- solve.QP(Dmat,dvec,Amat=A.eq,bvec=b.eq)$solution
#         # construct gross exposure portfilio
#         y_tmp <- x0[,index]%*%sol - y0
#         x_tmp <- as.numeric(x0[,index]%*%sol) - x0
#         cv <- cv.glmnet(x_tmp,y_tmp,parallel = TRUE)
#         w <- cv$glmnet.fit$beta[,which(cv$lambda == cv$lambda.min)]
#         w[index] <- w[index] + (1-sum(cv$glmnet.fit$beta[,which(cv$lambda == cv$lambda.min)]))*sol
#         w
#       }
#       w <- colMeans(w)
#   if (length(w) != p0) {
#     w <- insert.at(w,j-1,0)
#   }
#   return(w)
# }

#' Perform Portfolio Optimization Using Clusters and LASSO
#'
#' This function performs portfolio optimization using clustering and LASSO regularization.
#'
#' @param y0 A numeric vector of response values.
#' @param x0 A numeric matrix of predictors.
#' @param group A list of asset clusters.
#' @param rep The number of repetitions for optimization.
#' @return A numeric vector of optimized portfolio weights.
#' @export
po.bhu <- function (y0, x0, group, rep){
  p0 <- dim(x0)[2]
  j <- which(apply(x0, 2, var) <= 0, arr.ind = TRUE)
  x0 <- x0[,which(apply(x0, 2, var) > 0, arr.ind = TRUE)]
  #p <- dim(x0)[2]
  w <- foreach (i=1:rep, .combine = rbind, .packages=c('glmnet','quadprog'), .export='insert.at') %dopar% {
    ind <- rep(NA,length(group))
    for (k in (1:length(group))){
      ind[k] <- group[[k]][sample(1:length(group[[k]]), 1)]
    }
    Dmat <- cov(x0[,ind])
    dvec <- as.matrix(rep(0, length(group)))
    A.ge <- diag(length(group))
    b.ge <- as.matrix(rep(0, length(group)))
    A.eq <- as.matrix(rep(1, length(group)))
    b.eq <- 1
    sol <- solve.QP(Dmat,dvec,Amat=cbind(A.ge,A.eq),bvec=c(b.ge,b.eq))$solution
    if (any(round(sol,5) == 1)){
      k <- which(round(sol,5)==1)
      k <- ind[k]
      y_tmp <- x0[,k] - y0
      x_tmp <- x0[,k] - x0[,-k]
      cv <- cv.glmnet(x_tmp,y_tmp,parallel = TRUE)
      w.tmp <- cv$glmnet.fit$beta[,which(cv$lambda == cv$lambda.min)]
      insert.at(w.tmp,k-1,1-sum(as.numeric(w.tmp)))
    } else {
      y_tmp <- x0[,ind]%*%sol - y0
      x_tmp <- as.numeric(x0[,ind]%*%sol) - x0
      cv <- cv.glmnet(x_tmp,y_tmp,parallel = TRUE)
      w.tmp <- cv$glmnet.fit$beta[,which(cv$lambda == cv$lambda.min)]
      w.tmp[ind] <- w.tmp[ind] + (1-as.numeric(sum(cv$glmnet.fit$beta[,which(cv$lambda == cv$lambda.min)])))*sol
      w.tmp
    }
  }
  w <- colMeans(w)
  if (length(w) != p0) {
    w <- insert.at(w,j-1,0)
  }
  return(w)
}

#' Perform Portfolio Optimization Using TZT Method
#'
#' This function performs portfolio optimization using the TZT method.
#'
#' @param x0 A numeric matrix of asset returns.
#' @param gamma A numeric parameter for the TZT method.
#' @return A numeric vector of optimized portfolio weights.
#' @export
po.TZT <- function (x0, gamma){
  p0 <- dim(x0)[2]
  j <- which(apply(x0, 2, var) <= 0, arr.ind = TRUE)
  x0 <- x0[,which(apply(x0, 2, var) > 0, arr.ind = TRUE)]
  p <- dim(x0)[2]
  t <- dim(x0)[1]
  Dmat <- matrix(nearPD(cov(x0))["mat"][[1]]@x,p,p) #cov(x0)
  dvec <- as.matrix(rep(0, p))
  A.eq <- as.matrix(rep(1, p))
  b.eq <- 1
  w.mw <- solve.QP(Dmat,dvec,Amat=A.eq,bvec=b.eq)$solution
  w.eq <- rep(1/p,p)
  theta.sq <- colMeans(x0)%*%solve(Dmat)%*%colMeans(x0)
  c1 <- ((t-2)*(t-p-2))/((t-p-1)*(t-p-4))
  pi.1 <- w.eq%*%Dmat%*%w.eq - (2/gamma)*w.eq%*%colMeans(x0) + (1/gamma^2)*theta.sq
  pi.2 <- (1/gamma^2)*(c1-1)*theta.sq + (c1/gamma^2)*(p/t)
  delta <- pi.1/(pi.1+pi.2)
  w <- (1-delta)%*%w.eq + delta%*%w.mw
  if (length(w) != p0) {
    w <- insert.at(w,j-1,0)
  }
  return(as.numeric(w))
}

#' Perform Stochastic Weight Portfolio Optimization
#'
#' This function performs stochastic weight portfolio optimization.
#'
#' @param x0 A numeric matrix of asset returns.
#' @param b Number of assets to select in each sample.
#' @param sample Number of random samples to generate.
#' @return A numeric vector of optimized portfolio weights.
#' @export
po.SW <- function (x0, b, sample){
  p <- dim(x0)[2]
  j <- which(apply(x0, 2, var) <= 0, arr.ind = TRUE)
  x0 <- x0[,which(apply(x0, 2, var) > 0, arr.ind = TRUE)]
  p0 <- dim(x0)[2]
  w <- foreach (i=1:sample, .combine = rbind, .packages=c('glmnet','quadprog'), .export='insert.at') %dopar% {
    index <- sample(seq(1,p0,1), b, replace = FALSE)
    Dmat <- cov(x0[,index])
    dvec <- as.matrix(rep(0, b))
    A.eq <- as.matrix(rep(1, b))
    b.eq <- 1
    sol  <- solve.QP(Dmat,dvec,Amat=A.eq,bvec=b.eq)$solution
    w <- rep(0,p0)
    w[index] <- sol
    w
  }
  w <- colMeans(w)
  if (length(w) != p) {
    w <- insert.at(w,j-1,0)
  }
  return(as.numeric(w))
}

#' Perform LASSO Regularization with Stochastic Weight Portfolio Optimization
#'
#' This function performs portfolio optimization using LASSO regularization and stochastic weight selection.
#'
#' @param y0 A numeric vector of response values.
#' @param x0 A numeric matrix of predictors.
#' @param b Number of assets to select in each sample.
#' @param sample Number of random samples to generate.
#' @return A numeric vector of optimized portfolio weights.
#' @export
po.SW.lasso <- function (y0, x0, b, sample){
  p <- dim(x0)[2]
  j <- which(apply(x0, 2, var) <= 0, arr.ind = TRUE)
  x0 <- x0[,which(apply(x0, 2, var) > 0, arr.ind = TRUE)]
  p0 <- dim(x0)[2]
  w <- foreach (i=1:sample, .combine = rbind, .packages=c('glmnet','quadprog'), .export='insert.at') %dopar% {
    ind <- sample(seq(1,p0,1), b, replace = FALSE)
    Dmat <- cov(x0[,ind])
    dvec <- as.matrix(rep(0, b))
    A.ge <- diag(b)
    b.ge <- as.matrix(rep(0, b))
    A.eq <- as.matrix(rep(1, b))
    b.eq <- 1
    sol  <- solve.QP(Dmat,dvec,Amat=cbind(A.ge,A.eq),bvec=c(b.ge,b.eq))$solution
    if (any(round(sol,5) == 1)){
      k <- which(round(sol,5)==1)
      k <- ind[k]
      y_tmp <- x0[,k] - y0
      x_tmp <- x0[,k] - x0[,-k]
      cv <- cv.glmnet(x_tmp,y_tmp,parallel = TRUE)
      w.tmp <- cv$glmnet.fit$beta[,which(cv$lambda == cv$lambda.min)]
      insert.at(w.tmp,k-1,1-sum(as.numeric(w.tmp)))
    } else {
      y_tmp <- x0[,ind]%*%sol - y0
      x_tmp <- as.numeric(x0[,ind]%*%sol) - x0
      cv <- cv.glmnet(x_tmp,y_tmp,parallel = TRUE)
      w.tmp <- cv$glmnet.fit$beta[,which(cv$lambda == cv$lambda.min)]
      w.tmp[ind] <- w.tmp[ind] + (1-as.numeric(sum(cv$glmnet.fit$beta[,which(cv$lambda == cv$lambda.min)])))*sol
      w.tmp
    }
  }
  w <- colMeans(w)
  if (length(w) != p) {
    w <- insert.at(w,j-1,0)
  }
  return(as.numeric(w))
}
