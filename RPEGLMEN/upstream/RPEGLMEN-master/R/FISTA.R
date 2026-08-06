## This is the R implementation of FISTA
FISTA = function(x0, f, g, f_grad, g_prox, ABSTOL = 1e-8, maxiter = 1000){
  lambda = 1
  beta = 0.5
  
  y = x0
  xprev = x0
  h.fast_optval = rep(0, maxiter)
  t = 1
  tnext = t
  
  for(k in 1:maxiter){
    # y = x + (k/(k+3)) * (x - xprev)
    
    # line search
    while(TRUE){
      grad_y = f_grad(y)
      x = g_prox(y - lambda * grad_y, lambda)
      if (f(x) <= f(y) + sum(grad_y * (x - y)) + (1/(2*lambda)) * sum((x-y)^2))
        break
      lambda = beta * lambda
    }
    
    # compute objective function value
    h.fast_optval[k] = f(x) + g(x)
    
    # FISTA iteration
    tnext = (1 + sqrt(1 + 4 * t^2)) / 2
    y = x + (t - 1) / tnext * (x - xprev)
    xprev = x
    
    # stop if converged
    if ( k > 1)
      if (abs(h.fast_optval[k] - h.fast_optval[k-1]) < ABSTOL)
        break
  }
  return(list(x = x, objvals = h.fast_optval[1:k]))
}


# proximal of L1 norm
prox_L1 = function(x, threshold) {
  return(sapply(x, function(x)
    sign(x) * max(abs(x) - threshold, 0)))
}

# Elastic Net regularizer in  the form of
# lambda.EN * ( alpha.EN * ||x||_1 + (1 - alpha.EN) / 2 * ||x||_2^2)

regularizer_EN = function(x, lambda.EN, alpha.EN){
  lambda.EN * ( alpha.EN * sum(abs(x)) + (1 - alpha.EN) / 2 * sum(x^2))
}

# proximal of the Elastic Net regularizer in the form of
# lambda.EN * ( alpha.EN * ||x||_1 + (1 - alpha.EN) / 2 * ||x||_2^2)
prox_EN = function(x, lambda, lambda.EN, alpha.EN){
  1 / (1 + lambda * (1 - alpha.EN) * lambda.EN) * prox_L1(x, lambda * alpha.EN * lambda.EN)
}

# cost function of MSE

mse = function(x, A, b){
  0.5 * sum((A %*% x - b)^2)
}

# grad of mse

mse_grad = function(x, A, b){
  t(A) %*% A %*% x - t(A) %*% b
}

# Negative loglikelihood of glm gamma for fixed shape parameter
nll_gamma_glm = function(x, A, b, shape){
  n.vars = length(x)
  linear_predictor <- A %*% x
  # rate = shape / mean:
  rate <- shape / exp(linear_predictor)
  # sum of negative log likelihoods:
  -sum(dgamma(b, rate = rate, shape = shape,
              log = TRUE)) / length(b)
}

# the gradient of nll_gamma_glm

nll_gamma_glm_grad = function(x, A, b, shape){
  scale_vec = exp(A %*% x) / shape
  t(A) %*% ( shape - b / scale_vec) / length(b)
}

# find the smallest lambda that makes solution zero for Gamma_GLM_EN
# note the normalizing factor N

find_lambda_max_Gamma_GLM_EN = function(A,
                                        b,
                                        alpha.EN,
                                        shape){
  return(max(abs(shape/length(b)*(b-1)%*%A))/alpha.EN)
}

# generate vector of candidate lambdas

generate_lambda_grid = function(A, b, alpha, shape,
                                n.lambda = 100,
                                min.lambda.ratio = 1e-4
){
  lambda_max = find_lambda_max_Gamma_GLM_EN(A, b, alpha, shape)
  lambdas = 10^(seq(log10(lambda_max * min.lambda.ratio), log10(lambda_max), length.out = n.lambda))
}

# NLL of gamma glm to be used for optim
# the last coeff is the shape, rest is linear coeffs
# x input data, y output data

nll_gamma_optim = function(coeffs, x, y){
  n.vars = length(coeffs)
  # a = coeffs[1]
  shape = coeffs[n.vars]
  b = coeffs[1:(n.vars - 1)]
  linear_predictor <- x %*% b
  # rate = shape / mean:
  rate <- shape / exp(linear_predictor)
  # sum of negative log likelihoods:
  -sum(dgamma(y, rate = rate, shape = shape,
              log = TRUE)) / length(y)
}

# MLE of gamma GLM with shape unknown using optim BFGS
# x0 should contain nvar + 1 elements, last one must be positive

fit_glm_gamma = function(x0, x, y){
  res.optim = optim(x0, nll_gamma_optim, x = x, y = y, method = "BFGS")
}

# wrapper function to fit glmGammaNet for fixed shape, lambda.EN and alpha.EN

glmGammaNet = function(x0, A, b, lambda.EN, shape0, alpha.EN = 0.5, ABSTOL = 1e-8, maxiter = 1000){
  # First fit mle without regularization to estimate shape parameter
  # res.glm.gamma = fit_glm_gamma(c(rep(0, ncol(A)), 1), A, b)
  # res.glm.gamma.coeffs = res.glm.gamma$par
  # x0 = head(res.glm.gamma.coeffs, length(res.glm.gamma.coeffs) - 1)
  # shape0 = tail(res.glm.gamma.coeffs, 1)
  f = function(x){
    nll_gamma_glm(x, A, b, shape0)
  }
  
  f_grad = function(x){
    nll_gamma_glm_grad(x, A, b, shape0)
  }
  
  g = function(x){
    regularizer_EN(x, lambda.EN = lambda.EN, alpha.EN = alpha.EN)
  }
  
  g_prox = function(x, lambda){
    prox_EN(x, lambda = lambda, lambda.EN = lambda.EN, alpha.EN = alpha.EN)
  }
  res.FISTA = FISTA(x0, f, g, f_grad, g_prox, ABSTOL = ABSTOL, maxiter = maxiter)
  return(list(x = res.FISTA$x, shape = shape0, objvals = res.FISTA$objvals))
}


# cross validation for glmGammaNet

cv.glmGammaNet = function(A, b, ..., alpha.EN = 0.5, nfolds = 10,
                          nlambda = 100, ABSTOL = 1e-8,
                          maxiter = 1000, min.lambda.ratio = 1e-4){
  # First fit mle without regularization to estimate shape parameter
  res.glm.gamma = fit_glm_gamma(c(rep(0, ncol(A)), 1), A, b)
  res.glm.gamma.coeffs = res.glm.gamma$par
  x0 = head(res.glm.gamma.coeffs, length(res.glm.gamma.coeffs) - 1)
  shape0 = tail(res.glm.gamma.coeffs, 1)
  lambda.EN.vec = generate_lambda_grid(A, b, alpha.EN,
                                       shape0, n.lambda = nlambda,
                                       min.lambda.ratio = min.lambda.ratio)
  NLL_vec = rep(Inf, nlambda)
  NLL_sds = rep(Inf, nlambda)
  for(i in 1:nlambda){
    lambda.EN = lambda.EN.vec[i]
    NLLs = rep(Inf, nfolds)
    for(k in 1:nfolds){
      train.idx = sample(1:length(b), length(b) * (nfolds - 1) / nfolds)
      A.train = A[train.idx,]
      b.train = b[train.idx]
      A.test = A[-train.idx,]
      b.test = b[-train.idx]
      train.fit = glmGammaNet(x0,
                              A.train,
                              b.train,
                              lambda.EN = lambda.EN,
                              shape0 = shape0,
                              alpha.EN = alpha.EN,
                              ABSTOL = ABSTOL,
                              maxiter = maxiter)
      NLLs[k] = nll_gamma_glm(train.fit$x, A.test, b.test, shape0)
    }
    NLL_vec[i] = mean(NLLs)
    NLL_sds[i] = sd(NLLs)
  }
  best.idx = which.min(NLL_vec)
  # best.idx = sort(NLL_vec, index.return = TRUE)$ix[smallest.idx]
  best.lambda.EN = lambda.EN.vec[best.idx]
  best.fit = glmGammaNet(x0, A, b,
                         lambda.EN = best.lambda.EN,
                         shape0 = shape0,
                         alpha.EN = alpha.EN,
                         ABSTOL = ABSTOL,
                         maxiter = maxiter)
  return(best.fit$x)
}

