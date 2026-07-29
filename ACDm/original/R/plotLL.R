plotLL <- function(fitModel, parameter1 = 1, parameter2 = NULL, param1sequence = NULL, param2sequence = NULL, startpoint = NULL,
                   length.out = NULL, returnOutput = FALSE){
  
  logLikelihood <- NULL
  
  if(!is.null(length.out) && (!is.null(param1sequence) || !is.null(param1sequence))) warning("'length.out' is not used when paramXsequence is supplied")
   
  
  parameters <- parameter1
  if(length(parameter2) != 0){
    
    #checks if the required graphical package 'rgl' is installed
    if (!requireNamespace("rgl", quietly = TRUE)) {
      stop("The package 'rgl' needs to be installed to plot with thwo parameters. 
           Install it by entering 'install.packages(\"rgl\")' in the console, and then try again.",
           call. = FALSE)
    }
    
    parameters <- c(parameters, parameter2)
  } 
  #checks if the parameters were correctly entered, and transforms parameter string to int
  for(i in 1:length(parameters)){
    if(!is.numeric(parameters[i])){
      if(!is.character(parameters[i])){
        msg <- paste0("The argument 'parameter", i, "' must be an integer or a string of characters")
        stop(msg)
      } else if (parameters[i] %in% names(stats::coef(fitModel))){
        parameters[i] <-which(names(stats::coef(fitModel)) == parameters[i])
      } else {
        msg <- paste0("Couldn't find parameter", i, ": \"", parameters[i], "\"")
        stop("Couldn't find parameter1 \"", parameters[i], "\"")
      }
    }
    if(!(parameters[i] %in% 1:length(stats::coef(fitModel)))){
      msg <- paste0("There are only ", length(stats::coef(fitModel)), " parameters in the fitted model")
      stop(msg)
    }
  }
  
  if(!is.numeric(parameters)) parameters <- as.numeric(parameters)
  
  #computes the mean of the durations:
  if(length(fitModel$durations$adjDur) != 0){ 
    mean <- mean(fitModel$durations$adjDur) 
    dur <- fitModel$durations$adjDur
  } else {
    mean <- mean(fitModel$durations$durations) 
    dur <- fitModel$durations$durations
  }
  
  #fills the matrix 'exVar' with the exogenous variables if they are present in the model:
  exVar <- NULL 
  if(length(fitModel$exogenousVariables) != 0)
    exVar <- as.matrix(fitModel$durations[ , fitModel$exogenousVariables])
  
  if(length(parameter2) == 0){
    if(missing(param1sequence)){
      middle <- stats::coef(fitModel)[parameters[1]]
      se <- fitModel$parameterInference$SE[parameters[1]]
      if(is.na(se)){
        if("robustSE" %in% names(fitModel$parameterInference)){
          se <- fitModel$parameterInference$robustSE[parameters[1]]
        } else{
          msg <- paste0("Can't use default values for the range of '",
                        names(stats::coef(fitModel))[parameters[1]], 
                        "' since it's standard error is 'NA', 'param1sequence' must be given instead")
          stop(msg)
        }
      }
      param1sequence <- seq(from = middle - se*4, to = middle + se*4, length.out = ifelse(is.null(length.out), 21, length.out))
    }
  } else{
    if(missing(param1sequence)){
      middle <- stats::coef(fitModel)[parameters[1]]
      se <- fitModel$parameterInference$SE[parameters[1]]
      if(is.na(se)){
        if("robustSE" %in% names(fitModel$parameterInference)){
          se <- fitModel$parameterInference$robustSE[parameters[1]]
        } else{
          msg <- paste0("Can't use default values for the range of '",
                        names(stats::coef(fitModel))[parameters[1]], 
                        "' since it's standard error is 'NA', 'param1sequence' must be given instead")
          stop(msg)
        }
      }
      param1sequence <- seq(from = middle - se*4, to = middle + se*4, length.out = ifelse(is.null(length.out), 21, length.out))
    }
    if(missing(param2sequence)){
      middle <- stats::coef(fitModel)[parameters[2]]
      se <- fitModel$parameterInference$SE[parameters[2]]
      if(is.na(se)){
        if("robustSE" %in% names(fitModel$parameterInference)){
          se <- fitModel$parameterInference$robustSE[parameters[2]]
        } else{
          msg <- paste0("Can't use default values for the range of '",
                        names(stats::coef(fitModel))[parameters[2]], 
                        "' since it's standard error is 'NA', 'param2sequence' must be given instead")
          stop(msg)
        }
      }
      param2sequence <- seq(from = middle - se*4, to = middle + se*4, length.out = ifelse(is.null(length.out), 21, length.out))
    }
  }

  
  getLL.args <- list(dur = dur, 
                     exogenousVar = exVar , 
                     model = fitModel$model, 
                     order = fitModel$order,
                     mean = mean, 
                     distCode = fitModel$distCode, 
                     newDay = fitModel$newDayVector, 
                     returnMu = FALSE,
                     breakPoints = fitModel$breakPoints, 
                     forceErrExpec = fitModel$forceErrExpec, 
                     fixedParam = NULL, 
                     fixedParamPos = NULL, 
                     trace = 0,
                     returnIndex = 1,
                     startType = 1)
  

  if(fitModel$model %in% c("TACD", "SNIACD", "LSNIACD")){
    getLL.args <- c(getLL.args, list(J = length(fitModel$breakPoints) + 1))
  } 
  
  if(length(parameters) == 1){ #if only one parameter should be plotted:
    #function to compute the loglikelihood for diffrent values of one parameter, holding the others fixed:
    f1 <- function(x){
      
      #sets the fixed parameters to either the fitted values or the 'startpoint' values
      if(length(startpoint) == 0){
        internalCoef <- stats::coef(fitModel)
      } else {
        internalCoef <- startpoint
      }
      
      internalCoef[parameters[1]] <- x
      #-.getLLcall(param = internalCoef, getLL.args = getLL.args)
      -.computeLLcpp(param = internalCoef, getLL.args = getLL.args)
    } 
    
    #creates a data.frame of the log likelihood and the parameter values:
    df <- data.frame(param1sequence = param1sequence, logLikelihood = sapply(param1sequence, FUN = f1))
    g <- ggplot2::ggplot(df, aes(x = param1sequence, y = logLikelihood)) 
    if(min(param1sequence) <= stats::coef(fitModel)[parameters[1]] && 
       max(param1sequence) >= stats::coef(fitModel)[parameters[1]]){
      g <- g + ggplot2::geom_point(x = stats::coef(fitModel)[parameters[1]], y = fitModel$goodnessOfFit$value[1], color = "red", size = 3)
      maxll <- max(c(df$logLikelihood, fitModel$goodnessOfFit$value[1]), na.rm = TRUE)
      minll <- min(df$logLikelihood, na.rm = TRUE)
      g <- g + ggplot2::annotate("text", x = stats::coef(fitModel)[parameters[1]], y = maxll - (maxll-minll)/20, label = "MLE", color = "red")
    }
    g <- g + ggplot2::geom_point() + ggplot2::geom_line()
    g <- g + ggplot2::ylab("log likelihood") + ggplot2::xlab(names(stats::coef(fitModel))[parameters[1]])
    
    print(g)
    if(returnOutput) return(df)
    
  } else { #if two parameters should be plotted together:
    
    #function to compute the loglikelihood for diffrent values of the two parameters, holding the others fixed:
    f2 <- function(x, y){
      
      #sets the fixed parameters to either the fitted values or the 'startpoint' values
      if(length(startpoint) == 0){
        internalCoef <- stats::coef(fitModel)
      } else {
        internalCoef <- startpoint
      }
      
      internalCoef[parameter1] <- x
      internalCoef[parameter2] <- y
      -.getLLcall(param = internalCoef, getLL.args = getLL.args)
    }  
    
    
    #creates a matrix of the log likelihood values for diffrent values of the two parameters:
    z <- matrix(nrow = length(param1sequence), ncol = length(param2sequence))
    fillz <- function(){
      for(i in 1:length(param1sequence))  {
        for(j in 1:length(param2sequence))  {
          z[i, j] <<- f2(param1sequence[i], param2sequence[j])
        }
      }
    } 
    fillz()
    
    nbcol = 100
    color = grDevices::heat.colors(nbcol)
    z1 <- z
    z1[!is.finite(z)] <- NaN
    zcol  = cut(z1, nbcol)
    rgl::persp3d(x = param1sequence, y = param2sequence,  z = z1, col=color[zcol],
            ticktype="detailed", xlab = names(stats::coef(fitModel))[parameters[1]],
            ylab = names(stats::coef(fitModel))[parameters[2]], zlab = "log likelihood",axes=TRUE,
            phi = 30, theta = 40)
    
    xtemp <- rep(param1sequence, length(param2sequence))
    ytemp <- rep(param2sequence, each = length(param1sequence))
    ztemp <- as.vector(z1)
    
    rgl::points3d(x = xtemp, y = ytemp,  z = ztemp, color = "blue")
    
    min1 <- min(param1sequence); max1 <- max(param1sequence)
    min2 <- min(param2sequence); max2 <- max(param2sequence)
    mlepoint1 <- stats::coef(fitModel)[parameters[1]]
    mlepoint2 <- stats::coef(fitModel)[parameters[2]]
    
    if(min1 <= mlepoint1 && max1 >= mlepoint1 &&
       min2 <= mlepoint2 && max2 >= mlepoint2){
      
      
      rgl::points3d(x = mlepoint1, y = mlepoint2,  z = fitModel$goodnessOfFit$value[1], color = "red", size = 6)
    }
    
    if(returnOutput) return(list(para1 = param1sequence, para2 = param2sequence, z = z))
  }
}