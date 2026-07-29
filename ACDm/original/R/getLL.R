.getLLcall <- function(param, getLL.args = list()){  
  
  if(length(getLL.args$mean) == 0 || is.na(getLL.args$mean)) getLL.args$mean <- mean(getLL.args$dur, na.rm = TRUE)
  
  #combines the param and fixedParam into the full param vector if there are fixed parameters:
  if(length(getLL.args$fixedParamPos) != 0)    
    param <- .returnfixedPara(freePara = param, fixedParam = getLL.args$fixedParam, fixedParamPos = getLL.args$fixedParamPos)
  
  distPara <- .seperateStartPara(param, getLL.args$model, getLL.args$distCode, getLL.args$order, J = getLL.args$J)$distStartPara

  if(length(getLL.args$exogenousVar) == 0){
    if(getLL.args$model == "ACD"){
      temp <- .Call("getLL_ACDcall",
                    as.double(getLL.args$dur),
                    as.double(param),                     
                    as.integer(getLL.args$order),
                    as.double(getLL.args$mean),
                    as.integer(getLL.args$distCode),
                    as.double(distPara),
                    as.integer(getLL.args$newDay),
                    as.integer(getLL.args$forceErrExpec), PACKAGE = "ACDm")
    } else if(getLL.args$model == "LACD1"){
      temp <- .Call("getLL_LACD1call",
                    as.double(getLL.args$dur),
                    as.double(param),                     
                    as.integer(getLL.args$order),
                    as.double(getLL.args$mean),
                    as.integer(getLL.args$distCode),
                    as.double(distPara),
                    as.integer(getLL.args$newDay),
                    as.integer(getLL.args$forceErrExpec), PACKAGE = "ACDm")
    } else if(getLL.args$model == "LACD2"){
      temp <- .Call("getLL_LACD2call",
                    as.double(getLL.args$dur),
                    as.double(param),                     
                    as.integer(getLL.args$order),
                    as.double(getLL.args$mean),
                    as.integer(getLL.args$distCode),
                    as.double(distPara),
                    as.integer(getLL.args$newDay),
                    as.integer(getLL.args$forceErrExpec), PACKAGE = "ACDm")
    } else if(getLL.args$model == "EXACD"){
      temp <- .Call("getLL_EXACDcall",
                    as.double(getLL.args$dur),
                    as.double(param),                     
                    as.integer(getLL.args$order),
                    as.double(getLL.args$mean),
                    as.integer(getLL.args$distCode),
                    as.double(distPara),
                    as.integer(getLL.args$newDay),
                    as.integer(getLL.args$forceErrExpec), PACKAGE = "ACDm")
    } else if(getLL.args$model == "AMACD"){
      temp <- .Call("getLL_AMACDcall",
                    as.double(getLL.args$dur),
                    as.double(param),                     
                    as.integer(getLL.args$order),
                    as.double(getLL.args$mean),
                    as.integer(getLL.args$distCode),
                    as.double(distPara),
                    as.integer(getLL.args$newDay),
                    as.integer(getLL.args$forceErrExpec), PACKAGE = "ACDm")
    } else if(getLL.args$model == "AACD"){
      temp <- .Call("getLL_AACDcall",
                    as.double(getLL.args$dur),
                    as.double(param),                     
                    as.integer(getLL.args$order),
                    as.double(getLL.args$mean),
                    as.integer(getLL.args$distCode),
                    as.double(distPara),
                    as.integer(getLL.args$newDay),
                    as.integer(getLL.args$forceErrExpec), PACKAGE = "ACDm")
    } else if(getLL.args$model == "ABACD"){
      temp <- .Call("getLL_ABACDcall",
                    as.double(getLL.args$dur),
                    as.double(param),                     
                    as.integer(getLL.args$order),
                    as.double(getLL.args$mean),
                    as.integer(getLL.args$distCode),
                    as.double(distPara),
                    as.integer(getLL.args$newDay),
                    as.integer(getLL.args$forceErrExpec), PACKAGE = "ACDm")
    } else if(getLL.args$model == "BACD"){
      temp <- .Call("getLL_BACDcall",
                    as.double(getLL.args$dur),
                    as.double(param),                     
                    as.integer(getLL.args$order),
                    as.double(getLL.args$mean),
                    as.integer(getLL.args$distCode),
                    as.double(distPara),
                    as.integer(getLL.args$newDay),
                    as.integer(getLL.args$forceErrExpec), PACKAGE = "ACDm")
    } else if(getLL.args$model == "BCACD"){
      temp <- .Call("getLL_BCACDcall",
                    as.double(getLL.args$dur),
                    as.double(param),                     
                    as.integer(getLL.args$order),
                    as.double(getLL.args$mean),
                    as.integer(getLL.args$distCode),
                    as.double(distPara),
                    as.integer(getLL.args$newDay),
                    as.integer(getLL.args$forceErrExpec), PACKAGE = "ACDm")
    } else if(getLL.args$model == "SNIACD"){
      temp <- .Call("getLL_SNIACDcall",
                    as.double(getLL.args$dur),
                    as.double(param),                     
                    as.integer(getLL.args$order),
                    as.double(getLL.args$mean),
                    as.integer(getLL.args$distCode),
                    as.double(distPara),
                    as.integer(getLL.args$newDay),
                    as.double(getLL.args$breakPoints),
                    as.integer(getLL.args$forceErrExpec), PACKAGE = "ACDm")
    } else if(getLL.args$model == "LSNIACD"){
      temp <- .Call("getLL_logSNIACDcall",
                    as.double(getLL.args$dur),
                    as.double(param),                     
                    as.integer(getLL.args$order),
                    as.double(getLL.args$mean),
                    as.integer(getLL.args$distCode),
                    as.double(distPara),
                    as.integer(getLL.args$newDay),
                    as.double(getLL.args$breakPoints),
                    as.integer(getLL.args$forceErrExpec), PACKAGE = "ACDm")
    } else if(getLL.args$model == "TACD"){
      temp <- .Call("getLL_TACDcall",
                    as.double(getLL.args$dur),
                    as.double(getLL.args$dur),
                    as.integer(0),
                    as.double(getLL.args$breakPoints),
                    as.double(param),                     
                    as.integer(getLL.args$order),
                    as.double(getLL.args$mean),
                    as.integer(getLL.args$distCode),
                    as.double(distPara),
                    as.integer(getLL.args$newDay),
                    as.integer(getLL.args$forceErrExpec), PACKAGE = "ACDm")
    } else if(getLL.args$model == "TAMACD"){
      temp <- .Call("getLL_TAMACDcall",
                    as.double(getLL.args$dur),
                    as.double(getLL.args$dur),
                    as.integer(0),
                    as.double(getLL.args$breakPoints),
                    as.double(param),                     
                    as.integer(getLL.args$order),
                    as.double(getLL.args$mean),
                    as.integer(getLL.args$distCode),
                    as.double(distPara),
                    as.integer(getLL.args$newDay),
                    as.integer(getLL.args$forceErrExpec), PACKAGE = "ACDm")
    }
    
  } else { #if there are exogenous variables:
    
    if(getLL.args$model == "ACD"){
      temp <- .Call("getLL_ACDcallEx",
                    as.double(getLL.args$dur),
                    as.double(getLL.args$exogenousVar),
                    as.double(param),                     
                    as.integer(getLL.args$order),
                    as.double(getLL.args$mean),
                    as.integer(getLL.args$distCode),
                    as.double(distPara),
                    as.integer(getLL.args$newDay),
                    as.integer(getLL.args$forceErrExpec), PACKAGE = "ACDm")
    } else if(getLL.args$model == "LACD1"){
      temp <- .Call("getLL_LACD1callEx",
                    as.double(getLL.args$dur),
                    as.double(getLL.args$exogenousVar),
                    as.double(param),                     
                    as.integer(getLL.args$order),
                    as.double(getLL.args$mean),
                    as.integer(getLL.args$distCode),
                    as.double(distPara),
                    as.integer(getLL.args$newDay),
                    as.integer(getLL.args$forceErrExpec), PACKAGE = "ACDm")
    } else if(getLL.args$model == "LACD2"){
      temp <- .Call("getLL_LACD2callEx",
                    as.double(getLL.args$dur),
                    as.double(getLL.args$exogenousVar),
                    as.double(param),                     
                    as.integer(getLL.args$order),
                    as.double(getLL.args$mean),
                    as.integer(getLL.args$distCode),
                    as.double(distPara),
                    as.integer(getLL.args$newDay),
                    as.integer(getLL.args$forceErrExpec), PACKAGE = "ACDm")
    } else if(getLL.args$model == "EXACD"){
      temp <- .Call("getLL_EXACDcallEx",
                    as.double(getLL.args$dur),
                    as.double(getLL.args$exogenousVar),
                    as.double(param),                     
                    as.integer(getLL.args$order),
                    as.double(getLL.args$mean),
                    as.integer(getLL.args$distCode),
                    as.double(distPara),
                    as.integer(getLL.args$newDay),
                    as.integer(getLL.args$forceErrExpec), PACKAGE = "ACDm")
    } else if(getLL.args$model == "AMACD"){
      temp <- .Call("getLL_AMACDcallEx",
                    as.double(getLL.args$dur),
                    as.double(getLL.args$exogenousVar),
                    as.double(param),                     
                    as.integer(getLL.args$order),
                    as.double(getLL.args$mean),
                    as.integer(getLL.args$distCode),
                    as.double(distPara),
                    as.integer(getLL.args$newDay),
                    as.integer(getLL.args$forceErrExpec), PACKAGE = "ACDm")
    } else if(getLL.args$model == "AACD"){
      temp <- .Call("getLL_AACDcallEx",
                    as.double(getLL.args$dur),
                    as.double(getLL.args$exogenousVar),
                    as.double(param),                     
                    as.integer(getLL.args$order),
                    as.double(getLL.args$mean),
                    as.integer(getLL.args$distCode),
                    as.double(distPara),
                    as.integer(getLL.args$newDay),
                    as.integer(getLL.args$forceErrExpec), PACKAGE = "ACDm")
    } else if(getLL.args$model == "ABACD"){
      temp <- .Call("getLL_ABACDcallEx",
                    as.double(getLL.args$dur),
                    as.double(getLL.args$exogenousVar),
                    as.double(param),                     
                    as.integer(getLL.args$order),
                    as.double(getLL.args$mean),
                    as.integer(getLL.args$distCode),
                    as.double(distPara),
                    as.integer(getLL.args$newDay),
                    as.integer(getLL.args$forceErrExpec), PACKAGE = "ACDm")
    } else if(getLL.args$model == "BACD"){
      temp <- .Call("getLL_BACDcallEx",
                    as.double(getLL.args$dur),
                    as.double(getLL.args$exogenousVar),
                    as.double(param),                     
                    as.integer(getLL.args$order),
                    as.double(getLL.args$mean),
                    as.integer(getLL.args$distCode),
                    as.double(distPara),
                    as.integer(getLL.args$newDay),
                    as.integer(getLL.args$forceErrExpec), PACKAGE = "ACDm")
    }  else if(getLL.args$model == "BCACD"){
      temp <- .Call("getLL_BCACDcallEx",
                    as.double(getLL.args$dur),
                    as.double(getLL.args$exogenousVar),
                    as.double(param),                     
                    as.integer(getLL.args$order),
                    as.double(getLL.args$mean),
                    as.integer(getLL.args$distCode),
                    as.double(distPara),
                    as.integer(getLL.args$newDay),
                    as.integer(getLL.args$forceErrExpec), PACKAGE = "ACDm")
    } else if(getLL.args$model == "SNIACD"){
      temp <- .Call("getLL_SNIACDcallEx",
                    as.double(getLL.args$dur),
                    as.double(getLL.args$exogenousVar),
                    as.double(param),                     
                    as.integer(getLL.args$order),
                    as.double(getLL.args$mean),
                    as.integer(getLL.args$distCode),
                    as.double(distPara),
                    as.integer(getLL.args$newDay),
                    as.double(getLL.args$breakPoints),
                    as.integer(getLL.args$forceErrExpec), PACKAGE = "ACDm")
    } else if(getLL.args$model == "LSNIACD"){
      temp <- .Call("getLL_logSNIACDcallEx",
                    as.double(getLL.args$dur),
                    as.double(getLL.args$exogenousVar),
                    as.double(param),                     
                    as.integer(getLL.args$order),
                    as.double(getLL.args$mean),
                    as.integer(getLL.args$distCode),
                    as.double(distPara),
                    as.integer(getLL.args$newDay),
                    as.double(getLL.args$breakPoints),
                    as.integer(getLL.args$forceErrExpec), PACKAGE = "ACDm")
    } else if(getLL.args$model == "TACD"){
      temp <- .Call("getLL_TACDcallEx",
                    as.double(getLL.args$dur),
                    as.double(getLL.args$exogenousVar),
                    as.double(getLL.args$dur),
                    as.integer(0),
                    as.double(getLL.args$breakPoints),
                    as.double(param),                     
                    as.integer(getLL.args$order),
                    as.double(getLL.args$mean),
                    as.integer(getLL.args$distCode),
                    as.double(distPara),
                    as.integer(getLL.args$newDay),
                    as.integer(getLL.args$forceErrExpec), PACKAGE = "ACDm")
    } else if(getLL.args$model == "TAMACD"){
      temp <- .Call("getLL_TAMACDcallEx",
                    as.double(getLL.args$dur),
                    as.double(getLL.args$exogenousVar),
                    as.double(getLL.args$dur),
                    as.integer(0),
                    as.double(getLL.args$breakPoints),
                    as.double(param),                     
                    as.integer(getLL.args$order),
                    as.double(getLL.args$mean),
                    as.integer(getLL.args$distCode),
                    as.double(distPara),
                    as.integer(getLL.args$newDay),
                    as.integer(getLL.args$forceErrExpec), PACKAGE = "ACDm")
    }
  }
  
  if(getLL.args$returnMu){
    return(list(LL = temp[[3]], mu = temp[[1]], resi = temp[[2]]))
  }
  if(!getLL.args$returnMu){
    if(getLL.args$trace != 0) assign("ACDmOptimTrace", c(get("ACDmOptimTrace", envir = ACDmGlobalEnv), param, temp[[3]]), envir = ACDmGlobalEnv)
     return(-temp[[3]])
  }
  
}