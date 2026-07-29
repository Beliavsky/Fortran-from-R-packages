.computeLLcpp <- function(param, getLL.args = list()){  
  
  if(length(getLL.args$mean) == 0 || is.na(getLL.args$mean)) getLL.args$mean <- mean(getLL.args$dur, na.rm = TRUE)
  
  #combines the param and fixedParam into the full param vector if there are fixed parameters:
  if(length(getLL.args$fixedParamPos) != 0)    
    param <- .returnfixedPara(freePara = param, fixedParam = getLL.args$fixedParam, fixedParamPos = getLL.args$fixedParamPos)
  
  
  distPara <- .seperateStartPara(param, getLL.args$model, getLL.args$distCode, getLL.args$order, J = getLL.args$J)$distStartPara
  
  if(is.null(distPara)) distPara <- 0
  
  CppFunction <- switch(getLL.args$model,
                        ACD = .computeScoreACD,
                        LACD1 = .computeScoreLACD1,
                        LACD2 = .computeScoreLACD2,
                        EXACD = .computeScoreEXACD,
                        BACD = .computeScoreBACD,
                        BCACD = .computeScoreBCACD,
                        AMACD = .computeScoreAMACD,
                        ABACD = .computeScoreABACD,
                        SNIACD = .computeScoreSNIACD,
                        LSNIACD = .computeScoreLSNIACD,
                        TACD = .computeScoreTACD,
                        TAMACD = .computeScoreTAMACD,
                        AACD = .computeScoreAACD,
                        stop("model not supported"))
  
  if(is.null(getLL.args$returnIndex)) getLL.args$returnIndex <- 1
  
  if(getLL.args$model %in% c("SNIACD", "LSNIACD", "TACD", "TAMACD")){
    temp <- CppFunction(x = getLL.args$dur,
                        param = param,
                        order = getLL.args$order,
                        mean = getLL.args$mean,
                        dist = getLL.args$distCode,
                        distPara = distPara,
                        newDay = getLL.args$newDay,
                        forceErrExpec = getLL.args$forceErrExpec,
                        returnIndex = ifelse(getLL.args$returnIndex == 5, 3, getLL.args$returnIndex),
                        startType = getLL.args$startType,
                        bp = getLL.args$breakPoints)
  } else{
    temp <- CppFunction(x = getLL.args$dur,
                        param = param,
                        order = getLL.args$order,
                        mean = getLL.args$mean,
                        dist = getLL.args$distCode,
                        distPara = distPara,
                        newDay = getLL.args$newDay,
                        forceErrExpec = getLL.args$forceErrExpec,
                        returnIndex = ifelse(getLL.args$returnIndex == 5, 3, getLL.args$returnIndex),
                        startType = getLL.args$startType)
  }
  
  if(getLL.args$returnIndex == 5){
    
    unfixedParams <- rep(TRUE, length(param))
    if(length(getLL.args$fixedParamPos) != 0) if(sum(getLL.args$fixedParamPos) > 0) unfixedParams <- !getLL.args$fixedParamPos
    
    assign("score", -temp$score[unfixedParams], envir = ACDmGlobalEnv)
    assign("scoreParam", param, envir = ACDmGlobalEnv)
  } 
  
  if(getLL.args$trace != 0) 
    assign("ACDmOptimTrace", c(get("ACDmOptimTrace", envir = ACDmGlobalEnv), param, temp$LL), envir = ACDmGlobalEnv)
  
  if(getLL.args$returnIndex == 1 || getLL.args$returnIndex == 5) return(-temp$LL)
  else return(temp)
}

.getScoreCpp <- function(param, getLL.args = list()){  
  ## This function checks if the score was already computed, and if so returns it.
  ##  If not, the score is computed and then returned
  
  score <- get("score", envir = ACDmGlobalEnv)
  scoreParam <- get("scoreParam", envir = ACDmGlobalEnv)

  unfixedParams <- rep(TRUE, length(param))
  if(length(getLL.args$fixedParamPos) != 0) if(sum(getLL.args$fixedParamPos) > 0) unfixedParams <- !getLL.args$fixedParamPos
  
  getLL.args$returnIndex <- 3
  if(is.null(score)) return(-.computeLLcpp(param, getLL.args)$score[unfixedParams])
  if(is.null(scoreParam)) return(-.computeLLcpp(param, getLL.args)$score[unfixedParams])
  if(length(scoreParam) != length(param)) return(-.computeLLcpp(param, getLL.args)$score[unfixedParams])
  if(!isTRUE(all.equal(scoreParam, param))) return(-.computeLLcpp(param, getLL.args)$score[unfixedParams])
    
  return(score)
}

