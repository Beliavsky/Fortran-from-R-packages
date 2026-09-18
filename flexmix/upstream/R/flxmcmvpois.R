#
#  Copyright (C) 2004-2016 Friedrich Leisch and Bettina Gruen
#  $Id: flxmcmvpois.R 5079 2016-01-31 12:21:12Z gruen $
#

FLXMCmvpois <- function(formula=.~.)
{
  z <- new("FLXMC", weighted=TRUE, formula=formula,
                dist="mvpois", name="model-based Poisson clustering")
  
  z@preproc.y <- function(x){
    storage.mode(x) <- "integer"
    x
  }

  z@defineComponent <- function(para) {
    logLik <- function(x, y){
      colSums(dpois(t(y), para$lambda, log=TRUE))
    }

    predict <- function(x, ...){
      matrix(para$lambda, nrow = nrow(x), ncol=length(para$lambda),
             byrow=TRUE)
    }
    
    new("FLXcomponent", parameters=list(lambda=para$lambda), df=para$df,
        logLik=logLik, predict=predict)
  }
  z@fit <- function(x, y, w, ...){
      z@defineComponent(list(lambda = colSums(w*y)/sum(w), df = ncol(y)))
  }
  z 
}

