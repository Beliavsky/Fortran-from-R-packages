sim_ACD <- function(N = 1000, 
                    model = "ACD",
                    dist = "exponential",
                    param = NULL, 
                    order = NULL,
                    Nburn = 50,
                    startX = c(1),
                    startMu = c(1),
                    errors = NULL,
                    sampleErrors = TRUE,
                    roundToSec = FALSE,
                    rm0 = FALSE,
                    bp = NULL){
  

  if("acdFit" %in% class(model)){
    param <- coef(model)
    order <- model$order
    dist <- model$distribution
    bp <- model$breakPoints
    model <- model$model  
  }
  
  #provides the possibility of entering truncated and/or case mismatched arguments:
  model <- match.arg(toupper(model), c("ACD", "LACD1", "LACD2", "EXACD", "AMACD", "ABACD", "AACD", "BCACD", "BACD", "SNIACD", "LSNIACD", "TACD", "TAMACD"))
  dist <- match.arg(tolower(dist), c("exponential", "weibull", "burr", "gengamma", "genf"))
  
  # checks bp (breakpoints):
  if (model %in% c("SNIACD", "LSNIACD", "TACD", "TAMACD")) {
    if (is.null(bp)) {
      bp <- .setBP(1)
    }
    
    if (!is.numeric(bp) || !is.vector(bp)) {
      stop("'bp' must be a numeric vector")
    }
    
    if (any(bp < 0)) {
      stop("'bp' cannot contain negative values")
    }
    
    if (is.unsorted(bp)) {
      bp <- sort(bp)
      warning("'bp' was not sorted; it has been reordered to ascending")
    }
    
    J <- length(bp) + 1
  } else {
    J <- NULL
  }
  
  distCode <- .getDistCode(dist)
  #checks param and order input:
  if(length(param) != 0){
    if(length(order) == 0) order <- .setOrder(model)
    else .checkOrder(order, model)
    .checkPara(order, param, distCode, model)
    paraTemp <- .seperateStartPara(param, model, distCode, order, 0, J)
    distPara <- paraTemp$distStartPara
    startPara <- paraTemp$startPara
  }else{
    if(length(order) != 0){
      .checkOrder(order, model)
    } else{
      order <- .setOrder(model)
    }
    paraTemp <- .setStartPara(model, distCode, 1, order, 0, J)
    distPara <- paraTemp$distStartPara
    startPara <- paraTemp$startPara
    param <- c(distPara, startPara)
  }
  
  #simulates the error terms:
  if(length(errors) == 0){
    if(dist == "exponential") e <- stats::rexp(N + Nburn)
    else if(dist == "weibull"){
      e <- stats::rweibull(N + Nburn, shape = distPara, scale = 1/(gamma(1+1/distPara)))
    }
    else if(dist == "burr"){
      kappa <- distPara[1]
      sig2 <- distPara[2]
      muPara <- burrExpectation(theta = 1, kappa = kappa, sig2 = sig2)^kappa
      e <- rburr(N + Nburn, theta = 1, kappa = 1.2, sig2 = .3)
    } else if(dist == "gengamma"){
      kappa <- distPara[1]
      gammaPara <- distPara[2]
      
      e <- rgengamma(N + Nburn, gamma = gammaPara, kappa = kappa, forceExpectation = T)
    } else if(dist == "genf"){
      stop("Simulations are not available for the generelized F distribution")
    }
  } else{
    if(sampleErrors) e <- sample(errors, size = N + Nburn, replace = TRUE)
    else{
      if(length(errors) != N + Nburn) stop("the 'errors' vector needs to be of length N + Nburn if sampleErrors = FALSE")
      e <- errors
    }
  } 
  
  maxpq = max(order)
  #if the start value vector is smaller than the order, the start values are repeated:
  if(maxpq > min(length(startX), length(startMu))){
    startX <- rep(startX, length.out = maxpq)
    startMu <- rep(startMu, length.out = maxpq)
  } 
  
  
  if(model == "ACD"){
    temp <- .Call("sim_ACDCALL",
                  as.integer(N),
                  startPara,
                  order,
                  startX,
                  startMu,
                  e,
                  as.integer(Nburn), PACKAGE = "ACDm")
    
  } else if(model == "LACD1"){
    temp <- .Call("sim_LACD1",
                  as.integer(N),
                  startPara,
                  order,
                  startX,
                  startMu,
                  e,
                  as.integer(Nburn), PACKAGE = "ACDm")
    
  } else if(model == "LACD2"){
    temp <- .Call("sim_LACD2",
                  as.integer(N),
                  startPara,
                  order,
                  startX,
                  startMu,
                  e,
                  as.integer(Nburn), PACKAGE = "ACDm")
  } else if(model == "EXACD"){
    temp <- .Call("sim_EXACD",
                  as.integer(N),
                  startPara,
                  order,
                  startX,
                  startMu,
                  e,
                  as.integer(Nburn), PACKAGE = "ACDm")
  } else if(model == "AMACD"){
    temp <- .Call("sim_AMACD",
                  as.integer(N),
                  startPara,
                  order,
                  startX,
                  startMu,
                  e,
                  as.integer(Nburn), PACKAGE = "ACDm")
  } else if(model == "ABACD"){
    temp <- .Call("sim_ABACD",
                  as.integer(N),
                  startPara,
                  order,
                  startX,
                  startMu,
                  e,
                  as.integer(Nburn), PACKAGE = "ACDm")
  } else if(model == "AACD"){
    temp <- .Call("sim_AACD",
                  as.integer(N),
                  startPara,
                  order,
                  startX,
                  startMu,
                  e,
                  as.integer(Nburn), PACKAGE = "ACDm")
  } else if(model == "BCACD"){
    temp <- .Call("sim_BCACD",
                  as.integer(N),
                  startPara,
                  order,
                  startX,
                  startMu,
                  e,
                  as.integer(Nburn), PACKAGE = "ACDm")
  } else if(model == "BACD"){
    temp <- .Call("sim_BACD",
                  as.integer(N),
                  startPara,
                  order,
                  startX,
                  startMu,
                  e,
                  as.integer(Nburn), PACKAGE = "ACDm")
  } else if(model == "SNIACD"){
    temp <- .Call("sim_SNIACD",
                  as.integer(N),
                  startPara,
                  order,
                  startX,
                  startMu,
                  e,
                  as.integer(Nburn), 
                  bp, PACKAGE = "ACDm")
  } else if(model == "LSNIACD"){
    temp <- .Call("sim_LSNIACD",
                  as.integer(N),
                  startPara,
                  order,
                  startX,
                  startMu,
                  e,
                  as.integer(Nburn), 
                  bp, PACKAGE = "ACDm")
  } else if(model == "TACD"){
    temp <- .Call("sim_TACD",
                  as.integer(N),
                  startPara,
                  order,
                  startX,
                  startMu,
                  e,
                  as.integer(Nburn), 
                  bp, PACKAGE = "ACDm")
  } else if(model == "TAMACD"){
    temp <- .Call("sim_TAMACD",
                  as.integer(N),
                  startPara,
                  order,
                  startX,
                  startMu,
                  e,
                  as.integer(Nburn), 
                  bp, PACKAGE = "ACDm")
  }
  
  
  if(!roundToSec){ 
    return(temp)
  } else{
    temp <- ceiling(cumsum(temp))
    if(!rm0)  return(c(temp[1], temp[-1]-temp[-N]))
    else{
      temp <- c(temp[1], temp[-1]-temp[-N])
      return(temp[temp != 0])
    }
  }
}