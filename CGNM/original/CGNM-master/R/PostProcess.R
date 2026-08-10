## TODO make autodiagnosis function, is Lambda all big?  Is the minimum SSR parameter set within the initial range?
## TODO implement best parameter
## TODO implement worst accepted parameter


makeParaDistributionPlotDataFrame=function(CGNM_result, indicesToInclude=NA, cutoff_pvalue=0.05, numParametersIncluded=NA, ParameterNames=NA, ReparameterizationDef=NA, useAcceptedApproximateMinimizers=TRUE){

  ReReparameterise=TRUE

  Re_ReparameterizationDef=ReparameterizationDef
  Re_ParameterNames=ParameterNames

  if(is.null(CGNM_result$runSetting$ReparameterizationDef)){
    ReparameterizationDef=paste0("x", seq(1,length(CGNM_result$runSetting$initial_lowerRange)))
  }else{
    ReparameterizationDef=CGNM_result$runSetting$ReparameterizationDef
  }

  if(length(CGNM_result$runSetting$ParameterNames)!=length(ReparameterizationDef)){
    ParameterNames=ReparameterizationDef
  }else{
    ParameterNames=CGNM_result$runSetting$ParameterNames
  }

  if(is.na(Re_ReparameterizationDef)[1]){
    ReReparameterise=FALSE
  }




  if(length(Re_ParameterNames)!=length(Re_ReparameterizationDef)|is.na(Re_ParameterNames)[1]){
    Re_ParameterNames= Re_ReparameterizationDef
  }


  Kind_iter=NULL
  X_value=NULL
  cluster=NULL
  freeParaValues=data.frame()
  SSR_vec=CGNM_result$residual_history[,dim(CGNM_result$residual_history)[2]]


  # if(is.na(ParameterNames[1])|length(ParameterNames)!=length(ReparameterizationDef)){
  #   ParameterNames=ReparameterizationDef
  # }
  #
  # if(is.na(ParameterNames[1])){
  #   ParameterNames=paste0("x",seq(1,dim(CGNM_result$initialX)[2]))
  #   ReparameterizationDef=ParameterNames
  # }

  if(!is.na(indicesToInclude[1])){
    useIndecies_b=seq(1,dim(CGNM_result$X)[1])%in%indicesToInclude
  }else if(useAcceptedApproximateMinimizers){
    useIndecies_b=acceptedIndices_binary(CGNM_result,cutoff_pvalue, numParametersIncluded, useAcceptedApproximateMinimizers)
  }else{
    if(is.na(numParametersIncluded)|numParametersIncluded>dim(CGNM_result$X)[1]){
      useIndecies_b=rep(TRUE, dim(CGNM_result$X)[1])
    }else{
      useIndecies_b=(SSR_vec<=sort(SSR_vec)[numParametersIncluded])
    }
  }

  initialX_df=CGNM_result$initialX
  colnames(initialX_df)=paste0("x",seq(1,dim(CGNM_result$initialX)[2]))
  initialX_df=data.frame(initialX_df)

  repara_initialX_df=data.frame(row.names = seq(1,dim(initialX_df)[1]))
  temp_repara_initialX_df=data.frame(row.names = seq(1,dim(initialX_df)[1]))

  for(i in seq(1,length(ParameterNames))){
    temp_repara_initialX_df[,ParameterNames[i]]=with(initialX_df, eval(parse(text=ReparameterizationDef[i])))
  }
  if(ReReparameterise){
    for(i in seq(1,length(Re_ParameterNames))){
      repara_initialX_df[,Re_ParameterNames[i]]=with(temp_repara_initialX_df, eval(parse(text=Re_ReparameterizationDef[i])))
    }
  }else{
    repara_initialX_df=temp_repara_initialX_df

  }



  finalX_df=CGNM_result$X
  colnames(finalX_df)=paste0("x",seq(1,dim(CGNM_result$X)[2]))
  finalX_df=data.frame(finalX_df)

  repara_finalX_df=data.frame(row.names = seq(1,dim(finalX_df)[1]))
  temp_repara_finalX_df=data.frame(row.names = seq(1,dim(finalX_df)[1]))

  for(i in seq(1,length(ParameterNames))){
    temp_repara_finalX_df[,ParameterNames[i]]=with(finalX_df, eval(parse(text=ReparameterizationDef[i])))
  }
  if(ReReparameterise){
    for(i in seq(1,length(Re_ParameterNames))){
      repara_finalX_df[,Re_ParameterNames[i]]=with(temp_repara_finalX_df, eval(parse(text=Re_ReparameterizationDef[i])))
    }
  }else{
    repara_finalX_df=temp_repara_finalX_df

  }


  for(i in seq(1,dim(repara_initialX_df)[2])){
    freeParaValues=rbind(freeParaValues, data.frame(Name=names(repara_initialX_df)[i],X_value=repara_initialX_df[,i], Kind_iter="Initial", SSR=NA))
  }


  for(i in seq(1,dim(repara_finalX_df)[2])){
    freeParaValues=rbind(freeParaValues, data.frame(Name=names(repara_finalX_df)[i],X_value=repara_finalX_df[useIndecies_b,i], Kind_iter="Final Accepted", SSR=SSR_vec[useIndecies_b]))
  }

  if(!is.null(CGNM_result$bootstrapX)){

    bootstrapX_df=CGNM_result$bootstrapX
    colnames(bootstrapX_df)=paste0("x",seq(1,dim(CGNM_result$bootstrapX)[2]))
    bootstrapX_df=data.frame(bootstrapX_df)

    repara_bootstrapX_df=data.frame(row.names = seq(1,dim(bootstrapX_df)[1]))
    temp_repara_bootstrapX_df=data.frame(row.names = seq(1,dim(bootstrapX_df)[1]))
    for(i in seq(1,length(ParameterNames))){
      temp_repara_bootstrapX_df[,ParameterNames[i]]=with(bootstrapX_df, eval(parse(text=ReparameterizationDef[i])))
    }
    if(ReReparameterise){
      for(i in seq(1,length(Re_ParameterNames))){
        repara_bootstrapX_df[,Re_ParameterNames[i]]=with(temp_repara_bootstrapX_df, eval(parse(text=Re_ReparameterizationDef[i])))
      }

    }else{
      repara_bootstrapX_df=temp_repara_bootstrapX_df
    }

    for(i in seq(1,dim(repara_bootstrapX_df)[2])){
      freeParaValues=rbind(freeParaValues, data.frame(Name=names(repara_bootstrapX_df)[i],X_value=repara_bootstrapX_df[,i], Kind_iter="Bootstrap",SSR=NA))
    }

    freeParaValues$Kind_iter=factor(freeParaValues$Kind_iter, levels = c("Initial", "Final Accepted", "Bootstrap"))

  }else{
    freeParaValues$Kind_iter=factor(freeParaValues$Kind_iter, levels = c("Initial", "Final Accepted"))

  }

  freeParaValues$Name=factor(freeParaValues$Name, levels = names(repara_finalX_df))

  return(freeParaValues)
}

#' @title acceptedMaxSSR
#' @description
#' CGNM find multiple sets of minimizers of the nonlinear least squares (nls) problem by solving nls from various initial iterates.  Although CGNM is shown to be robust compared to other conventional multi-start algorithms, not all initial iterates minimizes successfully.  By assuming sum of squares residual (SSR) follows the chai-square distribution we first reject the approximated minimiser who SSR is statistically significantly worse than the minimum SSR found by the CGNM.  Then use elbow-method (a heuristic often used in mathematical optimisation to balance the quality and the quantity of the solution found) to find the "acceptable" maximum SSR.
#' @param CGNM_result (required input) \emph{A list} stores the computational result from Cluster_Gauss_Newton_method() function in CGNM package.
#' @param cutoff_pvalue (default: 0.05) \emph{A number} defines the rejection p-value for the first stage of acceptable computational result screening.
#' @param numParametersIncluded (default: NA) \emph{A natural number} defines the number of parameter sets to be included in the assessment of the acceptable parameters.  If set NA then use all the parameters found by the CGNM.
#' @param useAcceptedApproximateMinimizers (default: TRUE) \emph{TRUE or FALSE} If true then use chai-square and elbow method to choose maximum accepted SSR.  If false returnsnumParametersIncluded-th smallest SSR (or if numParametersIncluded=NA then returns the largest SSR).
#' @param algorithm (default: 2) \emph{1 or 2} specify the algorithm used for obtain accepted approximate minimizers. (Algorithm 1 uses elbow method, Algorithm 2 uses Grubbs' Test for Outliers.)
#' @return \emph{A positive real number} that is the maximum sum of squares residual (SSR) the algorithm has selected to accept.
#' @examples
#'
#'model_analytic_function=function(x){
#'
#'  observation_time=c(0.1,0.2,0.4,0.6,1,2,3,6,12)
#'  Dose=1000
#'  F=1
#'
#'  ka=x[1]
#'  V1=x[2]
#'  CL_2=x[3]
#'  t=observation_time
#'
#'  Cp=ka*F*Dose/(V1*(ka-CL_2/V1))*(exp(-CL_2/V1*t)-exp(-ka*t))
#'
#'  log10(Cp)
#'}
#'
#' observation=log10(c(4.91, 8.65, 12.4, 18.7, 24.3, 24.5, 18.4, 4.66, 0.238))
#'
#' CGNM_result=Cluster_Gauss_Newton_method(
#' nonlinearFunction=model_analytic_function,
#' targetVector = observation,
#' initial_lowerRange = c(0.1,0.1,0.1), initial_upperRange =  c(10,10,10),
#' num_iter = 10, num_minimizersToFind = 100, saveLog = FALSE)
#'
#' acceptedMaxSSR(CGNM_result)
#' @export
acceptedMaxSSR=function(CGNM_result, cutoff_pvalue=0.05, numParametersIncluded=NA, useAcceptedApproximateMinimizers=TRUE, algorithm=2){

  SSR_vec=CGNM_result$residual_history[,dim(CGNM_result$residual_history)[2]]

  if(useAcceptedApproximateMinimizers){

    if(algorithm==2){

      numInInitialSet=dim(CGNM_result$residual_history)[1]

      orderedIndex=order(SSR_vec)

      bestIndex=topIndices(CGNM_result,1)

      ptest_vec=c()

      cumResidualFromBestFit=c()

      pvalue_vec=c()
      for(i in seq(2,numInInitialSet)){

        numInSample=i-1
        cumResidualFromBestFit=c(cumResidualFromBestFit, sum((CGNM_result$Y[orderedIndex[i-1],]-CGNM_result$Y[orderedIndex[i],])))

        if(numInSample>3){

          tStat=qt(1-cutoff_pvalue/(2*numInSample), df=(numInSample-2))

          #test statistics from Grubbs' Test for Outliers
          testStatistics=(numInSample-1)/sqrt(numInSample)*sqrt(tStat^2/(numInSample-2+tStat^2))

          #test hypothesis where the null hypothesis is that the last added sum of residual is outlier compared to the ones that are already accepted
          ptest_vec=c(ptest_vec,abs(mean(cumResidualFromBestFit)-cumResidualFromBestFit[length(cumResidualFromBestFit)])/sd(cumResidualFromBestFit)>testStatistics)

        }else{
          ptest_vec=c(ptest_vec,FALSE)
        }
      }

      acceptMaxSSR=sort(SSR_vec)[which(ptest_vec)[1]+1]
    }else{
      min_R=min(SSR_vec)
      minIndex=which(SSR_vec==min_R)[1]

      targetVector=as.numeric(CGNM_result$runSetting$targetVector)
      targetVector[is.na(targetVector)]=0

      residual_vec=CGNM_result$Y[minIndex,]-targetVector

      acceptMaxSSR=max(qchisq(1-cutoff_pvalue, df=length(residual_vec)) *(sd(residual_vec))^2+min_R, sqrt(sqrt(.Machine$double.eps)))

      accept_index=which(SSR_vec<acceptMaxSSR)

      accept_vec=as.vector(SSR_vec<acceptMaxSSR)

      numAccept=sum(accept_vec, na.rm = TRUE)



      if(!is.na(numParametersIncluded)&(sum(accept_vec)>numParametersIncluded)){
        sortedAcceptedSSR=sort(SSR_vec[accept_vec])[seq(1,numParametersIncluded)]
      }else{
        sortedAcceptedSSR=sort(SSR_vec[accept_vec])
      }


      trapizoido_area=c()
      for(i in seq(1,length(sortedAcceptedSSR)-1))
        trapizoido_area=c(trapizoido_area, (sortedAcceptedSSR[1]+sortedAcceptedSSR[i])*i+(sortedAcceptedSSR[i+1]+sortedAcceptedSSR[length(sortedAcceptedSSR)])*(length(sortedAcceptedSSR)-i))

      strictMaxAcceptSSR=sortedAcceptedSSR[which(trapizoido_area==min(trapizoido_area))]

      accept_index=which(SSR_vec<strictMaxAcceptSSR)
      accept_vec=as.vector(SSR_vec<strictMaxAcceptSSR)
      numAccept=sum(accept_vec, na.rm = TRUE)
      acceptMaxSSR=strictMaxAcceptSSR
    }
## numPara*log(R_nu)+qchisq(1-alpha,1)

  }else if(is.na(numParametersIncluded)|numParametersIncluded>length(SSR_vec)){
    acceptMaxSSR=max(SSR_vec, na.rm = TRUE)
  }else{
    acceptMaxSSR=sort(SSR_vec)[numParametersIncluded]
  }

  acceptMaxSSR=max(acceptMaxSSR, min(SSR_vec)+sqrt(.Machine$double.eps))
  return(acceptMaxSSR)
}

#' @title acceptedApproximateMinimizers
#' @description
#' CGNM find multiple sets of minimizers of the nonlinear least squares (nls) problem by solving nls from various initial iterates.  Although CGNM is shown to be robust compared to other conventional multi-start algorithms, not all initial iterates minimizes successfully.  By assuming sum of squares residual (SSR) follows the chai-square distribution we first reject the approximated minimiser who SSR is statistically significantly worse than the minimum SSR found by the CGNM.  Then use elbow-method (a heuristic often used in mathematical optimisation to balance the quality and the quantity of the solution found) to find the "acceptable" maximum SSR. This function outputs the acceptable approximate minimizers of the nonlinear least squares problem found by the CGNM.
#' @param CGNM_result (required input) \emph{A list} stores the computational result from Cluster_Gauss_Newton_method() function in CGNM package.
#' @param cutoff_pvalue (default: 0.05) \emph{A number} defines the rejection p-value for the first stage of acceptable computational result screening.
#' @param numParametersIncluded (default: NA) \emph{A natural number} defines the number of parameter sets to be included in the assessment of the acceptable parameters.  If set NA then use all the parameters found by the CGNM.
#' @param useAcceptedApproximateMinimizers (default: TRUE) \emph{TRUE or FALSE} If true then use chai-square and elbow method to choose maximum accepted SSR.  If false returns the parameters upto numParametersIncluded-th smallest SSR (or if numParametersIncluded=NA then use all the parameters found by the CGNM).
#' @param algorithm (default: 2) \emph{1 or 2} specify the algorithm used for obtain accepted approximate minimizers. (Algorithm 1 uses elbow method, Algorithm 2 uses Grubbs' Test for Outliers.)
#' @param ParameterNames (default: NA) \emph{A vector of strings} the user can supply so that these names are used when making the plot. (Note if it set as NA or vector of incorrect length then the parameters are named as theta1, theta2, ... or as in ReparameterizationDef)
#' @param ReparameterizationDef (default: NA) \emph{A vector of strings} the user can supply definition of reparameterization where each string follows R syntax
#' @return \emph{A dataframe} that each row stores the accepted approximate minimizers found by CGNM.
#' @examples
#'
#'model_analytic_function=function(x){
#'
#'  observation_time=c(0.1,0.2,0.4,0.6,1,2,3,6,12)
#'  Dose=1000
#'  F=1
#'
#'  ka=x[1]
#'  V1=x[2]
#'  CL_2=x[3]
#'  t=observation_time
#'
#'  Cp=ka*F*Dose/(V1*(ka-CL_2/V1))*(exp(-CL_2/V1*t)-exp(-ka*t))
#'
#'  log10(Cp)
#'}
#'
#' observation=log10(c(4.91, 8.65, 12.4, 18.7, 24.3, 24.5, 18.4, 4.66, 0.238))
#'
#' CGNM_result=Cluster_Gauss_Newton_method(
#' nonlinearFunction=model_analytic_function,
#' targetVector = observation,
#' initial_lowerRange = c(0.1,0.1,0.1), initial_upperRange =  c(10,10,10),
#' num_iter = 10, num_minimizersToFind = 100, saveLog = FALSE)
#'
#' acceptedApproximateMinimizers(CGNM_result)
#' @export
acceptedApproximateMinimizers=function(CGNM_result, cutoff_pvalue=0.05, numParametersIncluded=NA, useAcceptedApproximateMinimizers=TRUE, algorithm=2, ParameterNames=NA, ReparameterizationDef=NA){
  out=CGNM_result$X[acceptedIndices(CGNM_result, cutoff_pvalue, numParametersIncluded, useAcceptedApproximateMinimizers, algorithm=algorithm),]

  out=data.frame(out)
  colnames(out)=paste0("x",seq(1,dim(CGNM_result$initialX)[2]))


  if(is.na(ParameterNames)[1]&!is.null(CGNM_result$runSetting$ParameterNames)){
    ParameterNames=CGNM_result$runSetting$ParameterNames
  }


  if(is.na(ReparameterizationDef)[1]&!is.null(CGNM_result$runSetting$ReparameterizationDef)){
    ReparameterizationDef=CGNM_result$runSetting$ReparameterizationDef
  }



  if(is.na(ParameterNames[1])|length(ParameterNames)!=length(ReparameterizationDef)){
    ParameterNames=ReparameterizationDef
  }

  if(is.na(ParameterNames[1])){
    ParameterNames=paste0("x",seq(1,dim(CGNM_result$initialX)[2]))
    ReparameterizationDef=ParameterNames
  }

  repara_finalX_df=data.frame(row.names = seq(1,dim(out)[1]))

  for(i in seq(1,length(ParameterNames))){
    repara_finalX_df[,ParameterNames[i]]=with(out, eval(parse(text=ReparameterizationDef[i])))
  }

  return(repara_finalX_df)
}

#' @title acceptedIndices
#' @description
#' CGNM find multiple sets of minimizers of the nonlinear least squares (nls) problem by solving nls from various initial iterates.  Although CGNM is shown to be robust compared to other conventional multi-start algorithms, not all initial iterates minimizes successfully.  By assuming sum of squares residual (SSR) follows the chai-square distribution we first reject the approximated minimiser who SSR is statistically significantly worse than the minimum SSR found by the CGNM.  Then use elbow-method (a heuristic often used in mathematical optimisation to balance the quality and the quantity of the solution found) to find the "acceptable" maximum SSR. This function outputs the indices of acceptable approximate minimizers of the nonlinear least squares problem found by the CGNM.
#' @param CGNM_result (required input) \emph{A list} stores the computational result from Cluster_Gauss_Newton_method() function in CGNM package.
#' @param cutoff_pvalue (default: 0.05) \emph{A number} defines the rejection p-value for the first stage of acceptable computational result screening.
#' @param numParametersIncluded (default: NA) \emph{A natural number} defines the number of parameter sets to be included in the assessment of the acceptable parameters.  If set NA then use all the parameters found by the CGNM.
#' @param useAcceptedApproximateMinimizers (default: TRUE) \emph{TRUE or FALSE} If true then use chai-square and elbow method to choose maximum accepted SSR.  If false returns the parameters upto numParametersIncluded-th smallest SSR (or if numParametersIncluded=NA then use all the parameters found by the CGNM).
#' @param algorithm (default: 2) \emph{1 or 2} specify the algorithm used for obtain accepted approximate minimizers. (Algorithm 1 uses elbow method, Algorithm 2 uses Grubbs' Test for Outliers.)
#' @return \emph{A vector of natural number} that contains the indices of accepted approximate minimizers found by CGNM.
#' @examples
#'
#'model_analytic_function=function(x){
#'
#'  observation_time=c(0.1,0.2,0.4,0.6,1,2,3,6,12)
#'  Dose=1000
#'  F=1
#'
#'  ka=x[1]
#'  V1=x[2]
#'  CL_2=x[3]
#'  t=observation_time
#'
#'  Cp=ka*F*Dose/(V1*(ka-CL_2/V1))*(exp(-CL_2/V1*t)-exp(-ka*t))
#'
#'  log10(Cp)
#'}
#'
#' observation=log10(c(4.91, 8.65, 12.4, 18.7, 24.3, 24.5, 18.4, 4.66, 0.238))
#'
#' CGNM_result=Cluster_Gauss_Newton_method(
#' nonlinearFunction=model_analytic_function,
#' targetVector = observation,
#' initial_lowerRange = c(0.1,0.1,0.1), initial_upperRange =  c(10,10,10),
#' num_iter = 10, num_minimizersToFind = 100, saveLog = FALSE)
#'
#' acceptedIndices(CGNM_result)
#' @export
acceptedIndices=function(CGNM_result, cutoff_pvalue=0.05, numParametersIncluded=NA, useAcceptedApproximateMinimizers=TRUE, algorithm=2){
  acceptMaxSSR_value=acceptedMaxSSR(CGNM_result, cutoff_pvalue, numParametersIncluded, useAcceptedApproximateMinimizers, algorithm=algorithm)

  SSR_vec=CGNM_result$residual_history[,dim(CGNM_result$residual_history)[2]]

  return(which(SSR_vec<=acceptMaxSSR_value))
}


#' @title topIndices
#' @description
#' CGNM find multiple sets of minimizers of the nonlinear least squares (nls) problem by solving nls from various initial iterates.  Although CGNM is shown to be robust compared to other conventional multi-start algorithms, not all initial iterates minimizes successfully.  One can visually inspect rank v.s. SSR plot and manually choose number of best fit acceptable parameters.  By using this function "topIndices", we can obtain the indices of the "numTopIndices" best fit parameter combinations.
#' @param CGNM_result (required input) \emph{A list} stores the computational result from Cluster_Gauss_Newton_method() function in CGNM package.
#' @param numTopIndices (required input) \emph{An integer} .
#' @return \emph{A vector of natural number} that contains the indices of accepted approximate minimizers found by CGNM.
#' @examples
#'
#'model_analytic_function=function(x){
#'
#'  observation_time=c(0.1,0.2,0.4,0.6,1,2,3,6,12)
#'  Dose=1000
#'  F=1
#'
#'  ka=x[1]
#'  V1=x[2]
#'  CL_2=x[3]
#'  t=observation_time
#'
#'  Cp=ka*F*Dose/(V1*(ka-CL_2/V1))*(exp(-CL_2/V1*t)-exp(-ka*t))
#'
#'  log10(Cp)
#'}
#'
#' observation=log10(c(4.91, 8.65, 12.4, 18.7, 24.3, 24.5, 18.4, 4.66, 0.238))
#'
#' CGNM_result=Cluster_Gauss_Newton_method(
#' nonlinearFunction=model_analytic_function,
#' targetVector = observation,
#' initial_lowerRange = c(0.1,0.1,0.1), initial_upperRange =  c(10,10,10),
#' num_iter = 10, num_minimizersToFind = 100, saveLog = FALSE)
#'
#' topInd=topIndices(CGNM_result, 10)
#'
#' ## This gives top 10 approximate minimizers
#' CGNM_result$X[topInd,]
#' @export
topIndices=function(CGNM_result, numTopIndices){

  SSR_vec=CGNM_result$residual_history[,dim(CGNM_result$residual_history)[2]]
  sortedSSR=sort(SSR_vec)
  outIndices=which(SSR_vec<=sortedSSR[numTopIndices])

  return(outIndices)
}



#' @title acceptedIndices_binary
#' @description
#' CGNM find multiple sets of minimizers of the nonlinear least squares (nls) problem by solving nls from various initial iterates.  Although CGNM is shown to be robust compared to other conventional multi-start algorithms, not all initial iterates minimizes successfully.  By assuming sum of squares residual (SSR) follows the chai-square distribution we first reject the approximated minimiser who SSR is statistically significantly worse than the minimum SSR found by the CGNM.  Then use elbow-method (a heuristic often used in mathematical optimisation to balance the quality and the quantity of the solution found) to find the "acceptable" maximum SSR. This function outputs the indices of acceptable approximate minimizers of the nonlinear least squares problem found by the CGNM. (note that acceptedIndices(CGNM_result) is equal to seq(1,length(acceptedIndices_binary(CGNM_result)))[acceptedIndices_binary(CGNM_result)])
#' @param CGNM_result (required input) \emph{A list} stores the computational result from Cluster_Gauss_Newton_method() function in CGNM package.
#' @param cutoff_pvalue (default: 0.05) \emph{A number} defines the rejection p-value for the first stage of acceptable computational result screening.
#' @param numParametersIncluded (default: NA) \emph{A natural number} defines the number of parameter sets to be included in the assessment of the acceptable parameters.  If set NA then use all the parameters found by the CGNM.
#' @param useAcceptedApproximateMinimizers (default: TRUE) \emph{TRUE or FALSE} If true then use chai-square and elbow method to choose maximum accepted SSR.  If false returns the indicies upto numParametersIncluded-th smallest SSR (or if numParametersIncluded=NA then use all the parameters found by the CGNM).
#' @param algorithm (default: 2) \emph{1 or 2} specify the algorithm used for obtain accepted approximate minimizers. (Algorithm 1 uses elbow method, Algorithm 2 uses Grubbs' Test for Outliers.)
#' @return \emph{A vector of TRUE and FALSE} that indicate if the each of the approximate minimizer found by CGNM is acceptable or not.
#' @examples
#'
#'model_analytic_function=function(x){
#'
#'  observation_time=c(0.1,0.2,0.4,0.6,1,2,3,6,12)
#'  Dose=1000
#'  F=1
#'
#'  ka=x[1]
#'  V1=x[2]
#'  CL_2=x[3]
#'  t=observation_time
#'
#'  Cp=ka*F*Dose/(V1*(ka-CL_2/V1))*(exp(-CL_2/V1*t)-exp(-ka*t))
#'
#'  log10(Cp)
#'}
#'
#' observation=log10(c(4.91, 8.65, 12.4, 18.7, 24.3, 24.5, 18.4, 4.66, 0.238))
#'
#' CGNM_result=Cluster_Gauss_Newton_method(
#' nonlinearFunction=model_analytic_function,
#' targetVector = observation,
#' initial_lowerRange = c(0.1,0.1,0.1), initial_upperRange =  c(10,10,10),
#' num_iter = 10, num_minimizersToFind = 100, saveLog = FALSE)
#'
#' acceptedIndices_binary(CGNM_result)
#' @export
acceptedIndices_binary=function(CGNM_result, cutoff_pvalue=0.05, numParametersIncluded=NA, useAcceptedApproximateMinimizers=TRUE, algorithm=2){
  acceptMaxSSR_value=acceptedMaxSSR(CGNM_result, cutoff_pvalue, numParametersIncluded, useAcceptedApproximateMinimizers, algorithm= algorithm)

  SSR_vec=CGNM_result$residual_history[,dim(CGNM_result$residual_history)[2]]

  return(as.vector(SSR_vec<=acceptMaxSSR_value))
}


#' @title plot_Rank_SSR
#' @description
#' Make SSR v.s. rank plot. This plot is often used to visualize the maximum accepted SSR.
#' @param CGNM_result (required input) \emph{A list} stores the computational result from Cluster_Gauss_Newton_method() function in CGNM package.
#' @param indicesToInclude (default: NA) \emph{A vector of integers} indices to include in the plot (if NA, use indices chosen by the acceptedIndices() function with default setting).
#' @return \emph{A ggplot object} of SSR v.s. rank.
#' @examples
#'
#'model_analytic_function=function(x){
#'
#'  observation_time=c(0.1,0.2,0.4,0.6,1,2,3,6,12)
#'  Dose=1000
#'  F=1
#'
#'  ka=x[1]
#'  V1=x[2]
#'  CL_2=x[3]
#'  t=observation_time
#'
#'  Cp=ka*F*Dose/(V1*(ka-CL_2/V1))*(exp(-CL_2/V1*t)-exp(-ka*t))
#'
#'  log10(Cp)
#'}
#'
#' observation=log10(c(4.91, 8.65, 12.4, 18.7, 24.3, 24.5, 18.4, 4.66, 0.238))
#'
#' CGNM_result=Cluster_Gauss_Newton_method(
#' nonlinearFunction=model_analytic_function,
#' targetVector = observation,
#' initial_lowerRange = c(0.1,0.1,0.1), initial_upperRange =  c(10,10,10),
#' num_iter = 10, num_minimizersToFind = 100, saveLog = FALSE)
#'
#' plot_Rank_SSR(CGNM_result)
#' @export
#' @import ggplot2
plot_Rank_SSR=function(CGNM_result, indicesToInclude=NA){
  SSR_value=NULL
  cutoff_pvalue=0.05
  numParametersIncluded=NA
  useAcceptedApproximateMinimizers=TRUE

  is_accepted=NULL
  SSR_vec=CGNM_result$residual_history[,dim(CGNM_result$residual_history)[2]]

  if(!is.na(indicesToInclude[1])){
    if(sum(indicesToInclude)<=dim(CGNM_result$X)[1]&length(indicesToInclude)==dim(CGNM_result$X)[1]){
      acceptedIndices_b=indicesToInclude
    }else{
      acceptedIndices_b=seq(1,dim(CGNM_result$X)[1]) %in% indicesToInclude
    }
  }else{
    acceptedIndices_b=acceptedIndices_binary(CGNM_result,cutoff_pvalue, numParametersIncluded, useAcceptedApproximateMinimizers)
  }

  numAccept=sum(acceptedIndices_b)
  acceptMaxSSR=max(SSR_vec[acceptedIndices_b])#acceptedMaxSSR(CGNM_result,cutoff_pvalue, numParametersIncluded, useAcceptedApproximateMinimizers)
  min_R=min(SSR_vec)

  plot_df=data.frame(SSR_value=SSR_vec, rank=rank(SSR_vec, na.last=TRUE, ties.method = "last"), is_accepted=acceptedIndices_b)

  ggplot2::ggplot(plot_df, ggplot2::aes(x=rank,y=SSR_value, colour=is_accepted))+ggplot2::geom_point()+ggplot2::coord_cartesian(ylim=c(0,acceptMaxSSR*2))+ggplot2::geom_vline(xintercept = numAccept, color="grey")+ ggplot2::annotate(geom="text", x=numAccept, y=acceptMaxSSR*0.5, label=paste("Accepted: ",numAccept,"\n Accepted max SSR: ",formatC(acceptMaxSSR, format = "g", digits = 3)),angle = 90,
                                                                                                                                                                                                                                       color="black")+ ggplot2::annotate(geom="text", x=length(SSR_vec)*0.1, y=min_R*1.1, label=paste("min SSR: ",formatC(min_R, format = "g", digits = 3)),
                                                                                                                                                                                                                                                                         color="black")+ylab("SSR")



}


#' @title plot_paraDistribution_byViolinPlots
#' @description
#' Make violin plot to compare the initial distribution and distribition of the accepted approximate minimizers found by the CGNM. Bars in the violin plots indicates the interquartile range. The solid line connects the interquartile ranges of the initial distribution and the distribution of the accepted approximate minimizer at the final iterate.  The blacklines connets the minimums and maximums of the initial distribution and the distribution of the accepted approximate minimizer at the final iterate. The black dots indicate the median.
#' @param CGNM_result (required input) \emph{A list} stores the computational result from Cluster_Gauss_Newton_method() function in CGNM package.
#' @param indicesToInclude (default: NA) \emph{A vector of integers} indices to include in the plot (if NA, use indices chosen by the acceptedIndices() function with default setting).
#' @param ParameterNames (default: NA) \emph{A vector of strings} the user can supply so that these names are used when making the plot. (Note if it set as NA or vector of incorrect length then the parameters are named as theta1, theta2, ... or as in ReparameterizationDef)
#' @param ReparameterizationDef (default: NA) \emph{A vector of strings} the user can supply definition of reparameterization where each string follows R syntax
#' @return \emph{A ggplot object} including the violin plot, interquartile range and median, minimum and maximum.
#' @examples
#'
#'model_analytic_function=function(x){
#'
#'  observation_time=c(0.1,0.2,0.4,0.6,1,2,3,6,12)
#'  Dose=1000
#'  F=1
#'
#'  ka=x[1]
#'  V1=x[2]
#'  CL_2=x[3]
#'  t=observation_time
#'
#'  Cp=ka*F*Dose/(V1*(ka-CL_2/V1))*(exp(-CL_2/V1*t)-exp(-ka*t))
#'
#'  log10(Cp)
#'}
#'
#' observation=log10(c(4.91, 8.65, 12.4, 18.7, 24.3, 24.5, 18.4, 4.66, 0.238))
#'
#'
#'
#' CGNM_result=Cluster_Gauss_Newton_method(
#' nonlinearFunction=model_analytic_function,
#' targetVector = observation,
#' initial_lowerRange = rep(0.01,3), initial_upperRange =  rep(100,3),
#' lowerBound=rep(0,3), ParameterNames = c("Ka","V1","CL"),
#' num_iter = 10, num_minimizersToFind = 100, saveLog = FALSE)
#'
#' plot_paraDistribution_byViolinPlots(CGNM_result)
#' plot_paraDistribution_byViolinPlots(CGNM_result,
#'      ReparameterizationDef=c("log10(Ka)","log10(V1)","log10(CL)"))
#'
#'
#' @export
#' @import ggplot2

plot_paraDistribution_byViolinPlots=function(CGNM_result, indicesToInclude=NA, ParameterNames=NA, ReparameterizationDef=NA){

  # if(is.na(ParameterNames)[1]&!is.null(CGNM_result$runSetting$ParameterNames)){
  #   ParameterNames=CGNM_result$runSetting$ParameterNames
  # }
  #
  #
  # if(is.na(ReparameterizationDef)[1]&!is.null(CGNM_result$runSetting$ReparameterizationDef)){
  #   ReparameterizationDef=CGNM_result$runSetting$ReparameterizationDef
  # }


  Kind_iter=NULL
  X_value=NULL

  cutoff_pvalue=0.05
  numParametersIncluded=NA
  useAcceptedApproximateMinimizers=TRUE

  freeParaValues=makeParaDistributionPlotDataFrame(CGNM_result, indicesToInclude, cutoff_pvalue, numParametersIncluded, ParameterNames, ReparameterizationDef, useAcceptedApproximateMinimizers)

  p<-ggplot2::ggplot(freeParaValues,ggplot2::aes(x=Kind_iter,y=X_value))+ggplot2::facet_wrap(Name~., scales = "free")

  p+ggplot2::geom_violin(trim=T,fill="#999999",linetype="blank",alpha=I(1/2))+
    ggplot2::stat_summary(geom="pointrange",fun = median, fun.min = function(x) quantile(x,probs=0.25), fun.max = function(x) quantile(x,probs=0.75), size=0.5,alpha=.5)+
    ggplot2::stat_summary(geom="line",fun =  function(x) quantile(x,probs=0), ggplot2::aes(group=1),size=0.5,alpha=.3,linetype=2)+
    ggplot2::stat_summary(geom="line",fun =  function(x) quantile(x,probs=1), ggplot2::aes(group=1),size=0.5,alpha=.3,linetype=2)+
    ggplot2::stat_summary(geom="line",fun =  function(x) quantile(x,probs=0.25), ggplot2::aes(group=1),size=0.5,alpha=.3)+
    ggplot2::stat_summary(geom="line",fun =  function(x) quantile(x,probs=0.75), ggplot2::aes(group=1),size=0.5,alpha=.3)+
    ggplot2::theme(legend.position="none",axis.text.x = element_text(angle = 25, hjust = 1))+xlab("")+ylab("Value")


}


#' @title plot_paraDistribution_byHistogram
#' @description
#' Make histograms to visualize the initial distribution and distribition of the accepted approximate minimizers found by the CGNM.
#' @param CGNM_result (required input) \emph{A list} stores the computational result from Cluster_Gauss_Newton_method() function in CGNM package.
#' @param indicesToInclude (default: NA) \emph{A vector of integers} indices to include in the plot (if NA, use indices chosen by the acceptedIndices() function with default setting).
#' @param ParameterNames (default: NA) \emph{A vector of strings} the user can supply so that these names are used when making the plot. (Note if it set as NA or vector of incorrect length then the parameters are named as theta1, theta2, ... or as in ReparameterizationDef)
#' @param ReparameterizationDef (default: NA) \emph{A vector of strings} the user can supply definition of reparameterization where each string follows R syntax
#' @param bins (default: 30) \emph{A natural number} Number of bins used for plotting histogram.
#' @return \emph{A ggplot object} including the violin plot, interquartile range and median, minimum and maximum.
#' @examples
#'
#'
#'model_analytic_function=function(x){
#'
#'  observation_time=c(0.1,0.2,0.4,0.6,1,2,3,6,12)
#'  Dose=1000
#'  F=1
#'
#'  ka=x[1]
#'  V1=x[2]
#'  CL_2=x[3]
#'  t=observation_time
#'
#'  Cp=ka*F*Dose/(V1*(ka-CL_2/V1))*(exp(-CL_2/V1*t)-exp(-ka*t))
#'
#'  log10(Cp)
#'}
#'
#' observation=log10(c(4.91, 8.65, 12.4, 18.7, 24.3, 24.5, 18.4, 4.66, 0.238))
#'
#'
#'
#' CGNM_result=Cluster_Gauss_Newton_method(
#' nonlinearFunction=model_analytic_function,
#' targetVector = observation,
#' initial_lowerRange = rep(0.01,3), initial_upperRange =  rep(100,3),
#' lowerBound=rep(0,3), ParameterNames = c("Ka","V1","CL"),
#' num_iter = 10, num_minimizersToFind = 100, saveLog = FALSE)
#'
#' plot_paraDistribution_byHistogram(CGNM_result)
#' plot_paraDistribution_byHistogram(CGNM_result,
#'      ReparameterizationDef=c("log10(Ka)","log10(V1)","log10(CL)"))
#'
#' @export
#' @import ggplot2

plot_paraDistribution_byHistogram=function(CGNM_result, indicesToInclude=NA, ParameterNames=NA, ReparameterizationDef=NA,  bins=30){

  # if(is.na(ParameterNames)[1]&!is.null(CGNM_result$runSetting$ParameterNames)){
  #
  #   ParameterNames=CGNM_result$runSetting$ParameterNames
  # }
  #
  #
  # if(is.na(ReparameterizationDef)[1]&!is.null(CGNM_result$runSetting$ReparameterizationDef)){
  #   ReparameterizationDef=CGNM_result$runSetting$ReparameterizationDef
  # }


  X_value=NULL
  Kind_iter=NULL

  cutoff_pvalue=0.05
  numParametersIncluded=NA
  useAcceptedApproximateMinimizers=TRUE
  freeParaValues=makeParaDistributionPlotDataFrame(CGNM_result, indicesToInclude, cutoff_pvalue, numParametersIncluded, ParameterNames, ReparameterizationDef, useAcceptedApproximateMinimizers)

  p<-ggplot2::ggplot(freeParaValues,ggplot2::aes(X_value))+ggplot2::geom_histogram(bins = bins)+ggplot2::facet_grid(Kind_iter~Name, scales = "free")+xlab("")

  p
}



#' @title plot_profileLikelihood
#' @description
#' Draw profile likelihood surface using the function evaluations conducted during CGNM computation. Note plot_SSRsurface can only be used when log is saved by setting saveLog=TRUE option when running Cluster_Gauss_Newton_method().  The grey horizontal line is the threshold for 95% pointwise confidence interval.
#' @param logLocation (required input) \emph{A string} of folder directory where CGNM computation log files exist.
#' @param alpha (default: 0.25) \emph{a number between 0 and 1} level of significance (used to draw horizontal line on the profile likelihood).
#' @param numBins (default: NA) \emph{A positive integer} SSR surface is plotted by finding the minimum SSR given one of the parameters is fixed and then repeat this for various values.  numBins specifies the number of different parameter values to fix for each parameter. (if set NA the number of bins are set as num_minimizersToFind/10)
#' @param ParameterNames (default: NA) \emph{A vector of strings} the user can supply so that these names are used when making the plot. (Note if it set as NA or vector of incorrect length then the parameters are named as theta1, theta2, ... or as in ReparameterizationDef)
#' @param ReparameterizationDef (default: NA) \emph{A vector of strings} the user can supply definition of reparameterization where each string follows R syntax
#' @param showInitialRange (default: TRUE) \emph{TRUE or FALSE} if TRUE then the initial range appears in the plot.
#' @param Likelihood_function (default: Residual_function_def) \emph{a function} that takes CGNM_result and initial then to calculate a quantity to be sketched in logscale e.g. SSR) this was implemented to conduct pos hoc drawing of the profile likelihood by providing the new definition of likelihood after all CGNM calculations are done.
#' @return \emph{A ggplot object} including the violin plot, interquartile range and median, minimum and maximum.
#' @examples
#'\dontrun{
#'model_analytic_function=function(x){
#'
#'  observation_time=c(0.1,0.2,0.4,0.6,1,2,3,6,12)
#'  Dose=1000
#'  F=1
#'
#'  ka=x[1]
#'  V1=x[2]
#'  CL_2=x[3]
#'  t=observation_time
#'
#'  Cp=ka*F*Dose/(V1*(ka-CL_2/V1))*(exp(-CL_2/V1*t)-exp(-ka*t))
#'
#'  log10(Cp)
#'}
#'
#' observation=log10(c(4.91, 8.65, 12.4, 18.7, 24.3, 24.5, 18.4, 4.66, 0.238))
#'
#' CGNM_result=Cluster_Gauss_Newton_method(
#' nonlinearFunction=model_analytic_function,
#' targetVector = observation,
#' initial_lowerRange = c(0.1,0.1,0.1), initial_upperRange =  c(10,10,10),
#' num_iter = 10, num_minimizersToFind = 100, saveLog=TRUE)
#'
#' plot_profileLikelihood("CGNM_log")
#' }
#' @export
#' @import ggplot2
plot_profileLikelihood=function(logLocation, alpha=0.25, numBins=NA,  ParameterNames=NA, ReparameterizationDef=NA, showInitialRange=TRUE, Likelihood_function=Residual_function_def){
  plot_SSRsurface(logLocation,  alpha= alpha, profile_likelihood=TRUE, numBins=numBins, ParameterNames=ParameterNames, ReparameterizationDef=ReparameterizationDef, showInitialRange=showInitialRange, Residual_function=Likelihood_function)
}



#' @title compare_profileLikelihood
#' @description
#' Draw profile likelihood surface using the function evaluations conducted during CGNM computation. Note plot_SSRsurface can only be used when log is saved by setting saveLog=TRUE option when running Cluster_Gauss_Newton_method().  The grey horizontal line is the threshold for 95% pointwise confidence interval.
#' @param logLocation (required input) \emph{List of strings} of folder directory where CGNM computation log files exist. (also can be list of CGNM_result objects)
#' @param alpha (default: 0.25) \emph{a number between 0 and 1} level of significance (used to draw horizontal line on the profile likelihood).
#' @param numBins (default: NA) \emph{A positive integer} SSR surface is plotted by finding the minimum SSR given one of the parameters is fixed and then repeat this for various values.  numBins specifies the number of different parameter values to fix for each parameter. (if set NA the number of bins are set as num_minimizersToFind/10)
#' @param ParameterNames (default: NA) \emph{A vector of strings} the user can supply so that these names are used when making the plot. (Note if it set as NA or vector of incorrect length then the parameters are named as theta1, theta2, ... or as in ReparameterizationDef)
#' @param ReparameterizationDef (default: NA) \emph{A vector of strings} the user can supply definition of reparameterization where each string follows R syntax
#' @param showInitialRange (default: TRUE) \emph{TRUE or FALSE} if TRUE then the initial range appears in the plot.
#' @param Likelihood_function (default: Residual_function_def) \emph{a function} that takes CGNM_result and initial then to calculate a quantity to be sketched in logscale e.g. SSR) this was implemented to conduct pos hoc drawing of the profile likelihood by providing the new definition of likelihood after all CGNM calculations are done.
#' @return \emph{A ggplot object} including the violin plot, interquartile range and median, minimum and maximum.
#' @examples
#'\dontrun{
#'model_analytic_function=function(x){
#'
#'  observation_time=c(0.1,0.2,0.4,0.6,1,2,3,6,12)
#'  Dose=1000
#'  F=1
#'
#'  ka=x[1]
#'  V1=x[2]
#'  CL_2=x[3]
#'  t=observation_time
#'
#'  Cp=ka*F*Dose/(V1*(ka-CL_2/V1))*(exp(-CL_2/V1*t)-exp(-ka*t))
#'
#'  log10(Cp)
#'}
#'
#' observation=log10(c(4.91, 8.65, 12.4, 18.7, 24.3, 24.5, 18.4, 4.66, 0.238))
#'
#' CGNM_result=Cluster_Gauss_Newton_method(
#' nonlinearFunction=model_analytic_function,
#' targetVector = observation,
#' initial_lowerRange = c(0.1,0.1,0.1), initial_upperRange =  c(10,10,10),
#' num_iter = 10, num_minimizersToFind = 100, saveLog=TRUE)
#'
#' plot_profileLikelihood("CGNM_log")
#' }
#' @export
#' @import ggplot2
compare_profileLikelihood=function(logLocation, alpha=0.25, numBins=NA,  ParameterNames=NA, ReparameterizationDef=NA, showInitialRange=TRUE, Likelihood_function=Residual_function_def){
  plot_SSRsurface(logLocation,  alpha= alpha, profile_likelihood=TRUE, numBins=numBins, ParameterNames=ParameterNames, ReparameterizationDef=ReparameterizationDef, showInitialRange=showInitialRange, Residual_function=Likelihood_function,
                  compareBetweenModels=TRUE
  )
}


#' @title table_profileLikelihoodConfidenceInterval
#' @description
#' Make table of confidence intervals that are approximated from the profile likelihood. First inspect profile likelihood plot and make sure the plot is smooth and has good enough resolution and the initial range is appropriate. Do not report this table without checking the profile likelihood plot.
#' @param logLocation (required input) \emph{A string or a list of strings} of folder directory where CGNM computation log files exist.
#' @param alpha (default: 0.25) \emph{a number between 0 and 1} level of significance used to derive the confidence interval.
#' @param numBins (default: NA) \emph{A positive integer} SSR surface is plotted by finding the minimum SSR given one of the parameters is fixed and then repeat this for various values.  numBins specifies the number of different parameter values to fix for each parameter. (if set NA the number of bins are set as num_minimizersToFind/10)
#' @param ParameterNames (default: NA) \emph{A vector of strings} the user can supply so that these names are used when making the plot. (Note if it set as NA or vector of incorrect length then the parameters are named as theta1, theta2, ... or as in ReparameterizationDef)
#' @param ReparameterizationDef (default: NA) \emph{A vector of strings} the user can supply definition of reparameterization where each string follows R syntax
#' @param pretty (default: FALSE) \emph{TRUE or FALSE} if true then the publication ready table will be an output
#' @param Likelihood_function (default: Residual_function_def) \emph{a function} that takes CGNM_result and initial then to calculate a quantity to be sketched in logscale e.g. SSR) this was implemented to conduct pos hoc drawing of the profile likelihood by providing the new definition of likelihood after all CGNM calculations are done.
#' @param silent (default: FALSE) \emph{TRUE or FALSE} set to TRUE to suppress warning and messages.
#' @return \emph{A ggplot object} including the violin plot, interquartile range and median, minimum and maximum.
#' @examples
#'\dontrun{
#'model_analytic_function=function(x){
#'
#'  observation_time=c(0.1,0.2,0.4,0.6,1,2,3,6,12)
#'  Dose=1000
#'  F=1
#'
#'  ka=x[1]
#'  V1=x[2]
#'  CL_2=x[3]
#'  t=observation_time
#'
#'  Cp=ka*F*Dose/(V1*(ka-CL_2/V1))*(exp(-CL_2/V1*t)-exp(-ka*t))
#'
#'  log10(Cp)
#'}
#'
#' observation=log10(c(4.91, 8.65, 12.4, 18.7, 24.3, 24.5, 18.4, 4.66, 0.238))
#'
#' CGNM_result=Cluster_Gauss_Newton_method(
#' nonlinearFunction=model_analytic_function,
#' targetVector = observation,
#' initial_lowerRange = c(0.1,0.1,0.1), initial_upperRange =  c(10,10,10),
#' num_iter = 10, num_minimizersToFind = 100, saveLog=TRUE)
#'
#' table_profileLikelihoodConfidenceInterval("CGNM_log")
#' }
#' @export
table_profileLikelihoodConfidenceInterval=function(logLocation, alpha=0.25, numBins=NA, ParameterNames=NA, ReparameterizationDef=NA, pretty=FALSE, Likelihood_function=Residual_function_def, silent=FALSE){
  if(!silent){
    print("WARNING: ALWAYS first inspect the profile likelihood plot (using plot_profileLikelihood()) and then use this table, DO NOT USE this table by itself.")
  }


  boundValue=NULL
  individual=NULL
  label=NULL
  minSSR=NULL
  negative2LogLikelihood=NULL
  newvalue=NULL

  parameterName=NULL

  reparaXinit=NULL
  value=NULL

  data=makeSSRsurfaceDataset(logLocation, TRUE, numBins, NA, ParameterNames, ReparameterizationDef, FALSE, Residual_function = Likelihood_function, silent=silent)


  likelihoodSurfacePlot_df=data$likelihoodSurfacePlot_df

  minSSR=min(likelihoodSurfacePlot_df$negative2LogLikelihood)


  likelihoodSurfacePlot_df$belowSignificance=(likelihoodSurfacePlot_df$negative2LogLikelihood-qchisq(1-alpha,1)-minSSR<0)

  cutOffValue=qchisq(1-alpha,1)+minSSR

  paraKind=unique(likelihoodSurfacePlot_df$parameterName)

  upperBound_vec=c()
  lowerBound_vec=c()
  min_vec=c()

  identifiability_vec=c()
  prettyString_vec=c()

  for(para_nu in paraKind){
    dataframe_nu=subset(likelihoodSurfacePlot_df,parameterName==para_nu)

    minIndex=which(dataframe_nu$negative2LogLikelihood==minSSR)


    upperBoundIndex=which(max(dataframe_nu$value[dataframe_nu$belowSignificance])==dataframe_nu$value)
    lowerBoundIndex=which(min(dataframe_nu$value[dataframe_nu$belowSignificance])==dataframe_nu$value)

    identifiability="identifiable"

    prettyString=""

    if(length(minIndex)>1){
      min_vec=c(min_vec, median(dataframe_nu$value[minIndex] ))
#      identifiability="Not identifiable"
    }else{
      min_vec=c(min_vec, dataframe_nu$value[minIndex] )
    }

    if(upperBoundIndex==length(dataframe_nu$value)){
      upperBound=paste0(">",dataframe_nu$value[upperBoundIndex])
      identifiability="Not identifiable"
      prettyString=paste(prettyString,paste0(">",signif(dataframe_nu$value[upperBoundIndex],2),")"))

    }else{

      diffValue=(dataframe_nu$value[upperBoundIndex]-dataframe_nu$value[upperBoundIndex+1])/(dataframe_nu$negative2LogLikelihood[upperBoundIndex]-dataframe_nu$negative2LogLikelihood[upperBoundIndex+1])*(cutOffValue-dataframe_nu$negative2LogLikelihood[upperBoundIndex])
      upperBound=dataframe_nu$value[upperBoundIndex]+diffValue

      #      upperBound=(dataframe_nu$value[upperBoundIndex]+dataframe_nu$value[upperBoundIndex+1])/2
      prettyString=paste(prettyString,paste0(signif(upperBound,2),"]"))
    }

    if(lowerBoundIndex==1){
      lowerBound=paste0("<",dataframe_nu$value[lowerBoundIndex])
      identifiability="Not identifiable"

      prettyString=paste(paste0("(<",signif(dataframe_nu$value[lowerBoundIndex],2),","), prettyString)

    }else{

      diffValue=(dataframe_nu$value[lowerBoundIndex]-dataframe_nu$value[lowerBoundIndex-1])/(dataframe_nu$negative2LogLikelihood[lowerBoundIndex]-dataframe_nu$negative2LogLikelihood[lowerBoundIndex-1])*(cutOffValue-dataframe_nu$negative2LogLikelihood[lowerBoundIndex])
      lowerBound=dataframe_nu$value[lowerBoundIndex]+diffValue
      #      lowerBound=(dataframe_nu$value[lowerBoundIndex]+dataframe_nu$value[lowerBoundIndex-1])/2

      prettyString=paste(paste0("[",signif(lowerBound,2),","), prettyString)

    }

    upperBound_vec=c(upperBound_vec,upperBound)
    lowerBound_vec=c(lowerBound_vec,lowerBound)
    prettyString_vec=c(prettyString_vec,prettyString)

    identifiability_vec=c(identifiability_vec,identifiability)
  }



  out_df=data.frame(parameterName=paraKind,CI_lower=lowerBound_vec,best=min_vec,CI_upper=upperBound_vec, parameter_identifiability=identifiability_vec)



  out_df[out_df$parameter_identifiability=="Not identifiable","best"]=NA


  pretty_df=data.frame(parameterName=paraKind,value=paste(signif(as.numeric(out_df$best), digits = 3), prettyString_vec))#, parameter_identifiability=identifiability_vec)

  names(pretty_df)=c("", paste0("best-fit [",alpha*100,"percentile, ", (1-alpha)*100,"percentile ]"))#,"identifiability")

  out_df=out_df[,seq(1,4)]
  names(out_df)=c("", paste(alpha*100,"percentile"), "best-fit", paste((1-alpha)*100,"percentile"))#,"identifiability")

  if(pretty){
    return(pretty_df)
  }else{
    return(out_df)
  }
}

#' @title plot_SSRsurface
#' @description
#' Make minimum SSR v.s. parameterValue plot using the function evaluations used during CGNM computation. Note plot_SSRsurface can only be used when log is saved by setting saveLog=TRUE option when running Cluster_Gauss_Newton_method().
#' @param logLocation (required input) \emph{A string or a list of strings} of folder directory where CGNM computation log files exist.
#' @param alpha (default: 0.25) \emph{a number between 0 and 1} level of significance (used to draw horizontal line on the profile likelihood).
#' @param profile_likelihood (default: FALSE) \emph{TRUE or FALSE} If set TRUE plot profile likelihood (assuming normal distribution of residual) instead of SSR surface.
#' @param numBins (default: NA) \emph{A positive integer} SSR surface is plotted by finding the minimum SSR given one of the parameters is fixed and then repeat this for various values.  numBins specifies the number of different parameter values to fix for each parameter. (if set NA the number of bins are set as num_minimizersToFind/10)
#' @param maxSSR (default: NA) \emph{A positive number} the maximum SSR that will be plotted on SSR surface plot.  This option is used to zoom into the SSR surface near the minimum SSR.
#' @param ParameterNames (default: NA) \emph{A vector of strings} the user can supply so that these names are used when making the plot. (Note if it set as NA or vector of incorrect length then the parameters are named as theta1, theta2, ... or as in ReparameterizationDef)
#' @param ReparameterizationDef (default: NA) \emph{A vector of strings} the user can supply definition of reparameterization where each string follows R syntax
#' @param showInitialRange (default: FALSE) \emph{TRUE or FALSE} if TRUE then the initial range appears in the plot.
#' @param Residual_function (default: Residual_function_def) \emph{a function} that takes CGNM_result and initial then to calculate a quantity to be sketched in logscale e.g. SSR) this was implemented to conduct pos hoc drawing of the profile likelihood by providing the new definition of likelihood after all CGNM calculations are done.
#' @param compareBetweenModels (default: FALSE) used for the implementation of compare_profileLikelihood so do not use when calling plot_SSRsurface function.
#' @return \emph{A ggplot object} including the violin plot, interquartile range and median, minimum and maximum.
#' @examples
#'\dontrun{
#'model_analytic_function=function(x){
#'
#'  observation_time=c(0.1,0.2,0.4,0.6,1,2,3,6,12)
#'  Dose=1000
#'  F=1
#'
#'  ka=x[1]
#'  V1=x[2]
#'  CL_2=x[3]
#'  t=observation_time
#'
#'  Cp=ka*F*Dose/(V1*(ka-CL_2/V1))*(exp(-CL_2/V1*t)-exp(-ka*t))
#'
#'  log10(Cp)
#'}
#'
#' observation=log10(c(4.91, 8.65, 12.4, 18.7, 24.3, 24.5, 18.4, 4.66, 0.238))
#'
#' CGNM_result=Cluster_Gauss_Newton_method(
#' nonlinearFunction=model_analytic_function,
#' targetVector = observation,
#' initial_lowerRange = c(0.1,0.1,0.1), initial_upperRange =  c(10,10,10),
#' num_iter = 10, num_minimizersToFind = 100, saveLog=TRUE)
#'
#' plot_SSRsurface("CGNM_log") + scale_y_continuous(trans='log10')
#' }
#' @export
#' @import ggplot2
plot_SSRsurface=function(logLocation, alpha=0.25,profile_likelihood=FALSE, numBins=NA, maxSSR=NA, ParameterNames=NA, ReparameterizationDef=NA, showInitialRange=FALSE, Residual_function=Residual_function_def, compareBetweenModels=FALSE){

  boundValue=NULL
  individual=NULL
  label=NULL
  minSSR=NULL
  negative2LogLikelihood=NULL
  newvalue=NULL
  maxBoundValue=NULL
  minBoundValue=NULL
  runName=NULL


  reparaXinit=NULL
  value=NULL

  if(compareBetweenModels){
    for(i in seq(1,length(logLocation))){
      if(i==1){
        data=makeSSRsurfaceDataset(logLocation[[i]], profile_likelihood, numBins, maxSSR, ParameterNames, ReparameterizationDef, showInitialRange, Residual_function=Residual_function)
        likelihoodSurfacePlot_df=data$likelihoodSurfacePlot_df
      }else{
        data=makeSSRsurfaceDataset(logLocation[[i]], profile_likelihood, numBins, maxSSR, ParameterNames, ReparameterizationDef, showInitialRange, Residual_function=Residual_function)
        likelihoodSurfacePlot_df=rbind(likelihoodSurfacePlot_df, data$likelihoodSurfacePlot_df)
      }
    }
  }else{

    data=makeSSRsurfaceDataset(logLocation, profile_likelihood, numBins, maxSSR, ParameterNames, ReparameterizationDef, showInitialRange, Residual_function=Residual_function)
    likelihoodSurfacePlot_df=data$likelihoodSurfacePlot_df
  }

 # data=makeSSRsurfaceDataset(logLocation, profile_likelihood, numBins, maxSSR, ParameterNames, ReparameterizationDef, showInitialRange, Residual_function=Residual_function)
  #likelihoodSurfacePlot_df=data$likelihoodSurfacePlot_df
  #residual_variance=data$residual_variance
  reparaXinit=data$initialX

  ParameterNames=data$ParameterNames
  ReparameterizationDef=data$ReparameterizationDef

  likelihoodSurfacePlot_df$parameterName=factor(likelihoodSurfacePlot_df$parameterName, levels=ParameterNames)

  if(profile_likelihood){
    #    likelihoodSurfacePlot_df=subset(likelihoodSurfacePlot_df, negative2LogLikelihood<=(2*(3.84+min(likelihoodSurfacePlot_df$negative2LogLikelihood))))

    if(compareBetweenModels){
      g=ggplot2::ggplot(likelihoodSurfacePlot_df, ggplot2::aes(x=value,y=negative2LogLikelihood, colour=runName))+ggplot2::geom_point()+ggplot2::geom_line()+
        ggplot2::coord_cartesian(ylim=c(-3.84+min(likelihoodSurfacePlot_df$negative2LogLikelihood, na.rm=TRUE), ((4*3.84+min(likelihoodSurfacePlot_df$negative2LogLikelihood,na.rm=TRUE)))))

    }else{
      g=ggplot2::ggplot(likelihoodSurfacePlot_df, ggplot2::aes(x=value,y=negative2LogLikelihood))+ggplot2::geom_point()+ggplot2::geom_line()+
        ggplot2::coord_cartesian(ylim=c(-3.84+min(likelihoodSurfacePlot_df$negative2LogLikelihood, na.rm=TRUE), ((4*3.84+min(likelihoodSurfacePlot_df$negative2LogLikelihood,na.rm=TRUE)))))

    }

  }else{
    g=ggplot2::ggplot(likelihoodSurfacePlot_df, ggplot2::aes(x=value,y=minSSR))+ggplot2::geom_point()+ggplot2::geom_line()

  }

  #  g=ggplot2::ggplot(likelihoodSurfacePlot_df, ggplot2::aes(x=value,y=minSSR))+ggplot2::geom_point()+ggplot2::geom_line()


  if(showInitialRange){
    Xinit_min=c()
    Xinit_max=c()

    for(i in seq(1,dim(reparaXinit)[2])){
      Xinit_min=c(Xinit_min,min(reparaXinit[,i]))
      Xinit_max=c(Xinit_max,max(reparaXinit[,i]))
    }

    g=g+ggplot2:: geom_rect(data=data.frame(minBoundValue=Xinit_min, maxBoundValue=Xinit_max,parameterName=factor(ParameterNames, levels=ParameterNames)), ggplot2::aes(xmin=minBoundValue, xmax=maxBoundValue, ymin=-Inf, ymax=Inf), alpha=0.3, fill="grey", inherit.aes = FALSE)
    g=g+ggplot2::geom_vline(data=data.frame(boundValue=Xinit_min, parameterName=factor(ParameterNames, levels=ParameterNames)), ggplot2::aes(xintercept=boundValue), colour="gray", inherit.aes = FALSE)
    g=g+ggplot2::geom_vline(data=data.frame(boundValue=Xinit_max, parameterName=factor(ParameterNames, levels=ParameterNames)), ggplot2::aes(xintercept=boundValue), colour="gray", inherit.aes = FALSE)
  }

  if(profile_likelihood){
    g=g+ggplot2::facet_wrap(.~parameterName,scales = "free_x")+ggplot2::ylab("-2log likelihood")+ggplot2::xlab("Parameter Value")
    g=g+ggplot2::geom_hline(yintercept = qchisq(1-alpha,1)+min(likelihoodSurfacePlot_df$negative2LogLikelihood), colour="grey")+
      ggplot2::labs(caption = paste0("Horizontal grey line is the threshold for confidence intervals corresponding to a confidence level of ",(1-alpha)*100,"%."))
  }else{
    g=g+ggplot2::facet_wrap(.~parameterName,scales = "free_x")+ggplot2::ylab("SSR")+ggplot2::xlab("Parameter Value")

  }
  return(g)
}

Residual_function_def=function(CGNM_result, initial=FALSE){
  matlabSum((Residual_matrix_function(CGNM_result, initial))^2,2)
}


Residual_matrix_function=function(CGNM_result, initial=FALSE){
  if(initial){
    out=CGNM_result$initialY[,!is.na(CGNM_result$runSetting$targetVector)]-repmat(matrix(CGNM_result$runSetting$targetVector[!is.na(CGNM_result$runSetting$targetVector)],nrow=1), dim(CGNM_result$initialY)[1],1)
  }else{
    out=CGNM_result$Y[,!is.na(CGNM_result$runSetting$targetVector)]-repmat(matrix(CGNM_result$runSetting$targetVector[!is.na(CGNM_result$runSetting$targetVector)],nrow=1), dim(CGNM_result$Y)[1],1)
  }
  return(out)
}

prepSSRsurfaceData=function(logLocation, ParameterNames=NA, ReparameterizationDef=NA,numBins=NA, Residual_function=Residual_function_def,silent=FALSE){

  CGNM_result=NULL
  minSSR=NULL
  negative2LogLikelihood=NULL
  value=NULL
  boundValue=NULL

  ReReparameterise=TRUE

  Re_ReparameterizationDef=ReparameterizationDef
  Re_ParameterNames=ParameterNames

  if(length(Re_ParameterNames)!=length(Re_ReparameterizationDef)|is.na(Re_ParameterNames)[1]){
    Re_ParameterNames=Re_ReparameterizationDef
  }


  if(is.na(Re_ReparameterizationDef)[1]){
    ReReparameterise=FALSE
  }


  if(typeof(logLocation)=="character"){
    DirectoryName_vec=c(logLocation)

  }else if(typeof(logLocation)=="list"){
    if(logLocation$runSetting$runName==""){
      DirectoryName_vec=c("CGNM_log","CGNM_log_bootstrap")
    }else{
      DirectoryName_vec=c(paste0(logLocation$runSetting$runName, "_CGNM_log"),paste0(logLocation$runSetting$runName, "_CGNM_log_bootstrap"))
    }



  }else{
    warning("DirectoryName need to be either the CGNM_result object or a string")
  }

  load(paste0(DirectoryName_vec[1],"/iteration_1.RDATA"))


  if(is.null(CGNM_result$runSetting$ReparameterizationDef)){
    ReparameterizationDef=paste0("x", seq(1,length(CGNM_result$runSetting$initial_lowerRange)))
  }else{
    ReparameterizationDef=CGNM_result$runSetting$ReparameterizationDef
  }

  if(length(CGNM_result$runSetting$ParameterNames)!=length(ReparameterizationDef)){
    ParameterNames=ReparameterizationDef
  }else{
    ParameterNames=CGNM_result$runSetting$ParameterNames
  }
  #
  #   if(is.na(ParameterNames)[1]&!is.null(CGNM_result$runSetting$ParameterNames)){
  #     ParameterNames=CGNM_result$runSetting$ParameterNames
  #   }
  #
  #
  #   if(is.na(ReparameterizationDef)[1]&!is.null(CGNM_result$runSetting$ReparameterizationDef)){
  #     ReparameterizationDef=CGNM_result$runSetting$ReparameterizationDef
  #   }
  #
  #
  #
  # if(is.na(ParameterNames[1])|length(ParameterNames)!=length(ReparameterizationDef)){
  #   ParameterNames=ReparameterizationDef
  # }

  initialX=CGNM_result$initialX

  numPara=dim(initialX)[2]

  if(is.na(ParameterNames[1])){
    ParameterNames=paste0("x",seq(1,dim(CGNM_result$initialX)[2]))
    ReparameterizationDef=ParameterNames
  }




  rawX_nu=CGNM_result$initialX
  #Residual_matrix=Residual_matrix_function(CGNM_result, initial=TRUE)
  R_nu=Residual_function(CGNM_result, initial=TRUE)#matlabSum((Residual_matrix)^2,2)

  #R_nu=CGNM_result$residual_history[,1]

  # if(is.null(CGNM_result$initialY))
  # {
  #   ResVar_nu=R_nu*NA
  # }else{
  #   Residual_matrix=Residual_matrix_function(CGNM_result)
  #   tempVariance=c()
  #   for(j in seq(1,dim(Residual_matrix)[1])){
  #     tempVariance=c(tempVariance, var(Residual_matrix[j,])*(numPara-1)/numPara)
  #   }
  #   ResVar_nu=tempVariance
  # }

  initiLNLfunc=CGNM_result$runSetting$nonlinearFunction
  rawXinit=CGNM_result$initialX
  runName= CGNM_result$runSetting$runName

  residual_variance_vec=c()

  for(DirectoryName in DirectoryName_vec){

    checkNL=TRUE
    for(i in seq(1,CGNM_result$runSetting$num_iteration)){
      fileName=paste0(DirectoryName,"/iteration_",i,".RDATA")
      if (file.exists(fileName)){
        load(fileName)
        rawX_nu=rbind(rawX_nu,CGNM_result$X)

        #Residual_matrix=Residual_matrix_function(CGNM_result)#$Y[,!is.na(CGNM_result$runSetting$targetVector)]-repmat(matrix(CGNM_result$runSetting$targetVector[!is.na(CGNM_result$runSetting$targetVector)],nrow=1), dim(CGNM_result$Y)[1],1)
        # tempVariance=c()
        # for(j in seq(1,dim(Residual_matrix)[1])){
        #   tempVariance=c(tempVariance,var(Residual_matrix[j,])*(numPara-1)/numPara)
        # }
        #
        # ResVar_nu=c(ResVar_nu,tempVariance)
        tempR=Residual_function(CGNM_result, initial=FALSE)#matlabSum((Residual_matrix)^2,2)
        R_nu=c(R_nu, tempR)

        if(checkNL){
          if(!silent){
            print(paste0("log saved in ",getwd(),"/",DirectoryName," is used to draw SSR/likelihood surface"))
          }
          # if(!identical(initiLNLfunc, CGNM_result$runSetting$nonlinearFunction)){
          #   warning(paste0("the nonlinear function used in this log in ",getwd(),"/",DirectoryName," is not the same as ",getwd(),"/",DirectoryName_vec[1]))
          # }
          checkNL=FALSE
        }

      }
    }

    #  residual_variance_vec=c(residual_variance_vec, min(R_nu,na.rm = TRUE)/(sum(!is.na(CGNM_result$runSetting$targetVector))-1))

  }
  #  residual_variance=min(residual_variance_vec, na.rm = TRUE)

  reparaXinit=data.frame(row.names = seq(1,dim(rawXinit)[1]))
  temp_reparaXinit=data.frame(row.names = seq(1,dim(rawXinit)[1]))

  colnames(rawXinit)=paste0("x",seq(1,dim(rawXinit)[2]))
  rawXinit=data.frame(rawXinit)

  for(i in seq(1,length(ParameterNames))){
    temp_reparaXinit[,ParameterNames[i]]=with(rawXinit, eval(parse(text=ReparameterizationDef[i])))
  }
  if(ReReparameterise){
    for(i in seq(1,length(Re_ParameterNames))){
      reparaXinit[,Re_ParameterNames[i]]=with(temp_reparaXinit, eval(parse(text=Re_ReparameterizationDef[i])))
    }
  }else{
    reparaXinit=temp_reparaXinit
  }

  #
  #   for(i in seq(1,CGNM_result$runSetting$num_iteration)){
  #     fileName=paste0(DirectoryName,"bootstrap/iteration_",i,".RDATA")
  #
  #
  #     if (file.exists(fileName)){
  #       load(fileName)
  #       rawX_nu=rbind(rawX_nu,CGNM_result$X)
  #
  #       R_temp=c()
  #       for(j in seq(1, dim(CGNM_result$Y)[1])){
  #         R_temp=c(R_temp,sum((CGNM_result$Y[j,]-CGNM_result$runSetting$targetVector)^2,na.rm = TRUE))
  #       }
  #
  #       R_nu=c(R_nu, R_temp)
  #     }
  #   }



  tempMatrix=unique(cbind(rawX_nu,R_nu))
  rawX_nu=tempMatrix[,seq(1,dim(rawX_nu)[2])]
  R_nu=tempMatrix[,dim(tempMatrix)[2]]

  # tempMatrix=unique(cbind(rawX_nu,R_nu,ResVar_nu))
  # rawX_nu=tempMatrix[,seq(1,dim(rawX_nu)[2])]
  # R_nu=tempMatrix[,dim(tempMatrix)[2]-1]
  # ResVar_nu=tempMatrix[,dim(tempMatrix)[2]]

  negative2LogLikelihood_nu=numPara*log(R_nu)#R_nu/ResVar_nu+numPara*log(ResVar_nu)
  #negative2LogLikelihood_nu=R_nu/min(ResVar_nu,na.rm = TRUE)#R_nu/ResVar_nu+numPara*log(ResVar_nu)

  colnames(rawX_nu)=paste0("x",seq(1,dim(rawX_nu)[2]))
  rawX_nu=data.frame(rawX_nu)

  X_nu=data.frame(row.names = seq(1,dim(rawX_nu)[1]))
  temp_X_nu=data.frame(row.names = seq(1,dim(rawX_nu)[1]))
  for(i in seq(1,length(ParameterNames))){
    temp_X_nu[,ParameterNames[i]]=with(rawX_nu, eval(parse(text=ReparameterizationDef[i])))
  }


  if(ReReparameterise){
    for(i in seq(1,length(Re_ParameterNames))){
      X_nu[,Re_ParameterNames[i]]=with(temp_X_nu, eval(parse(text=Re_ReparameterizationDef[i])))
    }
  }else{
    X_nu=temp_X_nu
  }

  if(is.na(numBins)){
    numBins=round(dim(X_nu)[1]/100)
  }
  #  residual_variance=residual_variance
  #  out=list(X_matrix=X_nu, R_vec=R_nu, ResVar_vec=ResVar_nu,negative2LogLikelihood=negative2LogLikelihood_nu,ParameterNames=names(X_nu),numBins=numBins, initialX=reparaXinit, ReparameterizationDef=Re_ReparameterizationDef)
  out=list(X_matrix=X_nu, R_vec=R_nu,negative2LogLikelihood=negative2LogLikelihood_nu,ParameterNames=names(X_nu),numBins=numBins, initialX=reparaXinit, ReparameterizationDef=Re_ReparameterizationDef, runName=runName)

  return(out)
}

#' @title plot_2DprofileLikelihood
#' @description
#' Make likelihood related values v.s. parameterValues plot using the function evaluations used during CGNM computation. Note plot_SSRsurface can only be used when log is saved by setting saveLog=TRUE option when running Cluster_Gauss_Newton_method().
#' @param logLocation (required input) \emph{A string or a list of strings} of folder directory where CGNM computation log files exist.
#' @param index_x (default: NA) \emph{A vector of strings or numbers} List parameter names or indices used for the surface plot. (if NA all parameters are used)
#' @param index_y (default: NA) \emph{A vector of strings or numbers} List parameter names or indices used for the surface plot. (if NA all parameters are used)
#' @param plotType (default: 2) \emph{A number 0,1,2,3, or 4} 0: number of model evaluations done, 1: (1-alpha) where alpha is the significance level, this plot is recommended for the ease of visualization as it ranges from 0 to 1. 2: -2log likelihood. 3: SSR. 4: all points within 1-alpha confidence region
#' @param plotMax (default: NA) \emph{A number} the maximum value that will be plotted on surface plot. (If NA all values are included in the plot, note SSR or likelihood can range many orders of magnitudes fo may want to restrict when plotting them)
#' @param numBins (default: NA) \emph{A positive integer} 2D profile likelihood surface is plotted by finding the minimum SSR given two of the parameters are fixed and then repeat this for various values.  numBins specifies the number of different parameter values to fix for each parameter. (if set NA the number of bins are set as num_minimizersToFind/10)
#' @param ParameterNames (default: NA) \emph{A vector of strings} the user can supply so that these names are used when making the plot. (Note if it set as NA or vector of incorrect length then the parameters are named as theta1, theta2, ... or as in ReparameterizationDef)
#' @param ReparameterizationDef (default: NA) \emph{A vector of strings} the user can supply definition of reparameterization where each string follows R syntax
#' @param showInitialRange (default: TRUE) \emph{TRUE or FALSE} if TRUE then the initial range appears in the plot.
#' @param alpha (default: 0.25) \emph{a number between 0 and 1} level of significance (all the points outside of this significance level will not be plotted when plot tyoe 1,2 or 4 are chosen).
#' @param Likelihood_function (default: Residual_function_def) \emph{a function} that takes CGNM_result and initial then to calculate a quantity to be sketched in logscale e.g. SSR) this was implemented to conduct pos hoc drawing of the profile likelihood by providing the new definition of likelihood after all CGNM calculations are done.
#' @param dependentVariable (default: NA) \emph{NA or a vector} can be specified as dependent variable for the parameter-parameter correlation plot (i.e., for the default case it is set as the -2 log-likelihood to plot 2D profile likelihood)
#' @return \emph{A ggplot object} including the violin plot, interquartile range and median, minimum and maximum.
#' @examples
#'\dontrun{
#'model_analytic_function=function(x){
#'
#'  observation_time=c(0.1,0.2,0.4,0.6,1,2,3,6,12)
#'  Dose=1000
#'  F=1
#'
#'  ka=10^x[1]
#'  V1=10^x[2]
#'  CL_2=10^x[3]
#'  t=observation_time
#'
#'  Cp=ka*F*Dose/(V1*(ka-CL_2/V1))*(exp(-CL_2/V1*t)-exp(-ka*t))
#'
#'  log10(Cp)
#'}
#'
#' observation=log10(c(4.91, 8.65, 12.4, 18.7, 24.3, 24.5, 18.4, 4.66, 0.238))
#'
#' CGNM_result=Cluster_Gauss_Newton_method(
#' nonlinearFunction=model_analytic_function,
#' targetVector = observation,
#' initial_lowerRange = c(-1,-1,-1), initial_upperRange =  c(1,1,1),
#' num_iter = 10, num_minimizersToFind = 500, saveLog=TRUE)
#'
#' ## the minimum example
#' plot_2DprofileLikelihood("CGNM_log")
#'
#' ## we can draw profilelikelihood also including bootstrap result
#'CGNM_result=Cluster_Gauss_Newton_Bootstrap_method(CGNM_result,
#'                       nonlinearFunction = model_analytic_function)
#'
#' ## example with various options
#' plot_2DprofileLikelihood(c("CGNM_log","CGNM_log_bootstrap"),
#'  showInitialRange = TRUE,index_x = c("ka","V1"))
#'  }
#' @export
#' @import ggplot2
plot_2DprofileLikelihood=function(logLocation, index_x=NA, index_y=NA, plotType=2,plotMax=NA, ParameterNames=NA, ReparameterizationDef=NA,numBins=NA, showInitialRange=TRUE, alpha=0.25, Likelihood_function=Residual_function_def, dependentVariable=NA){

  x_axis=NULL
  y_axis=NULL
  value=NULL
  x_min=NULL
  x_max=NULL
  y_min=NULL
  y_max=NULL

  if(plotType==1|plotType=="1-alpha"){
    plotType="1-alpha"
  }else if(plotType==2|plotType=="-2logLikelihood"){
    plotType="-2logLikelihood"
  }else if(plotType==3|plotType=="SSR"){
    plotType="SSR"
  }else if(plotType=="count"|plotType==0){
    plotType="count"
  }else if(plotType=="statSignificant"|plotType==4){
    plotType="statSignificant"
  }

  preppedDataset_list=prepSSRsurfaceData(logLocation, ParameterNames, ReparameterizationDef,numBins, Residual_function=Likelihood_function)

  Xinit_min=c()
  Xinit_max=c()
  for(i in seq(1,dim(preppedDataset_list$initialX)[2])){
    Xinit_min=c(Xinit_min,min(preppedDataset_list$initialX[,i]))
    Xinit_max=c(Xinit_max,max(preppedDataset_list$initialX[,i]))
  }

  index_x_vec=c()

  if(is.na(index_x)){
    index_x_vec=seq(1,length(preppedDataset_list$ParameterNames))
  }else if(is.numeric(index_x)){
    index_x_vec=index_x
  }else{
    for(index_nu in index_x){
      index_x_vec=c(index_x_vec, which(preppedDataset_list$ParameterNames==index_nu))
    }
  }

  index_y_vec=c()
  if(is.na(index_y)){
    index_y_vec=seq(1,length(preppedDataset_list$ParameterNames))
  }else if(is.numeric(index_y)){
    index_y_vec=index_y
  }else{
    for(index_nu in index_y){
      index_y_vec=c(index_y_vec, which(preppedDataset_list$ParameterNames==index_nu))
    }
  }

  X_val=preppedDataset_list$X_matrix
  numBins=preppedDataset_list$numBins

  SSR_vec=preppedDataset_list$R_vec
  neg2likelihood_vec=preppedDataset_list$negative2LogLikelihood


  if(plotType=="1-alpha"){
    z_value=pchisq(neg2likelihood_vec-min(neg2likelihood_vec),df=2)
  }else if(plotType=="-2logLikelihood"){
    z_value=neg2likelihood_vec
  }else if(plotType=="SSR"){
    z_value=SSR_vec
  }else if(plotType=="statSignificant"){
    z_value=as.numeric((neg2likelihood_vec-min(neg2likelihood_vec)-qchisq(0.95,1))<0)
  }


  n2likelihood=preppedDataset_list$negative2LogLikelihood

  if(!is.na(plotMax)){
    useIndex=(((n2likelihood-min(n2likelihood))<qchisq(1-alpha,2))&(z_value<plotMax))
  }else{
    useIndex=((n2likelihood-min(n2likelihood))<qchisq(1-alpha,2))
  }


  X_val=X_val[useIndex,]

  SSR_vec=SSR_vec[useIndex]
  neg2likelihood_vec=neg2likelihood_vec[useIndex]
  X_rounded=X_val


  for(i in seq(1,dim(X_rounded)[2])){

    # for(j in seq(1,numBins)){
    #   binUpper=quantile(X_val[,i], probs = c((j-1)*(1/numBins),(j)*(1/numBins)))[2]
    #   binLower=quantile(X_val[,i], probs = c((j-1)*(1/numBins),(j)*(1/numBins)))[1]
    #   indexToUse=((X_val[,i]<=binUpper)&(X_val[,i]>binLower))
    #
    #   X_rounded[indexToUse,i]=(binUpper+binLower)/2
    # }

    for(j in seq(1,numBins)){

      X_rounded[,i]=(X_val[,i]-min(X_val[,i]))/(max(X_val[,i])-min(X_val[,i]))
      X_rounded[,i]=round(X_rounded[,i]*numBins)/numBins
      X_rounded[,i]=X_rounded[,i]*(max(X_val[,i])-min(X_val[,i]))+min(X_val[,i])
    }
  }

  aggData_comnbined_df=data.frame()
  boundData_df=data.frame()

  X_forPlot4=unique(X_rounded)

  out <- tryCatch(
    {
      optimal_kmeans((X_forPlot4))

    },
    error=function(cond) {
      rep(1,dim(X_forPlot4)[1])
    },
    warning=function(cond) {
      optimal_kmeans((X_forPlot4))
    },
    finally={

    }
  )

  cluster_index=as.factor(out$cluster)

  for(index_x in index_x_vec){
    for(index_y in index_y_vec){

      if(index_x!=index_y){

        boundData_df=rbind(boundData_df, data.frame(x_min=Xinit_min[index_x],x_max=Xinit_max[index_x],y_min=Xinit_min[index_y],y_max=Xinit_max[index_y],x_lab=preppedDataset_list$ParameterNames[index_x],y_lab=preppedDataset_list$ParameterNames[index_y]))
        if(length(dependentVariable)==dim(X_forPlot4)[1]){
          aggData_df=data.frame(x_axis=X_forPlot4[,index_x],y_axis=X_forPlot4[,index_y],value=dependentVariable)

        }else if(plotType=="count"){

          aggData_df=data.frame(aggregate(SSR_vec, by=list(X_rounded[,index_x],X_rounded[,index_y]), FUN="length"))

          names(aggData_df)=c("x_axis","y_axis", "value")

        }else if(plotType=="statSignificant"){

          aggData_df=data.frame(x_axis=X_forPlot4[,index_x],y_axis=X_forPlot4[,index_y],value=cluster_index)

        }else if(plotType=="SSR"){


          aggData_df=data.frame(aggregate(SSR_vec, by=list(X_rounded[,index_x],X_rounded[,index_y]), FUN="min"))

          names(aggData_df)=c("x_axis","y_axis","SSR")
          aggData_df$value=aggData_df$SSR
        }else{

          aggData_df=data.frame(aggregate(neg2likelihood_vec, by=list(X_rounded[,index_x],X_rounded[,index_y]), FUN="min"))

          names(aggData_df)=c("x_axis","y_axis","neg2likelihood")

          if(plotType=="1-alpha"){
            aggData_df$value=pchisq(aggData_df$neg2likelihood-min(aggData_df$neg2likelihood),df=2)
          }else if(plotType=="-2logLikelihood"){

            aggData_df$value=aggData_df$neg2likelihood

            aggData_df$value[pchisq(aggData_df$neg2likelihood-min(aggData_df$neg2likelihood),df=2)>(1-alpha)]=NA

          }

        }

        aggData_df$x_lab=factor(preppedDataset_list$ParameterNames[index_x],levels=preppedDataset_list$ParameterNames[index_x])
        aggData_df$y_lab=factor(preppedDataset_list$ParameterNames[index_y],levels=preppedDataset_list$ParameterNames[index_y])

        aggData_comnbined_df=rbind(aggData_comnbined_df,aggData_df)
      }

    }
  }
  aggData_comnbined_df=aggData_comnbined_df[order(-aggData_comnbined_df$value),]


  if(plotType=="count"){

    g=ggplot2::ggplot(aggData_comnbined_df, ggplot2::aes(x=x_axis, y=y_axis, alpha=log10(value)))+ggplot2::geom_point()+ggplot2::xlab("")+ggplot2::ylab("")+ggplot2::labs(alpha="log10(count)")+ggplot2::facet_grid(y_lab~x_lab,scales = "free")

  }else if(plotType=="1-alpha"){
    g=ggplot2::ggplot(subset(aggData_comnbined_df, value<(1-alpha)), ggplot2::aes(x=x_axis, y=y_axis, z=value))+ggplot2::xlab("")+ggplot2::ylab("")+ggplot2::labs(fill="1-alpha")+ggplot2::facet_grid(y_lab~x_lab,scales = "free")+ggplot2::geom_contour_filled(bins=20)

  }else if(plotType=="statSignificant"){

    g=ggplot2::ggplot(aggData_comnbined_df, ggplot2::aes(x=x_axis, y=y_axis, colour=value))+ggplot2::xlab("")+ggplot2::ylab("")+ggplot2::geom_point()+ggplot2::facet_grid(y_lab~x_lab,scales = "free")

  }else{
    #g=ggplot2::ggplot(aggData_comnbined_df, ggplot2::aes(x=x_axis, y=y_axis, z=value))+ggplot2::xlab("")+ggplot2::ylab("")+ggplot2::labs(fill=plotType)+ggplot2::facet_grid(y_lab~x_lab,scales = "free")+ggplot2::geom_contour_filled(bins=20)

    g=ggplot2::ggplot(aggData_comnbined_df, ggplot2::aes(x=x_axis, y=y_axis, colour=value))+ggplot2::xlab("")+ggplot2::ylab("")+ggplot2::labs(colour=plotType)+ggplot2::facet_grid(y_lab~x_lab,scales = "free")+ggplot2::geom_point()

  }


  if(showInitialRange){

    g=g+ggplot2::geom_vline(data=boundData_df, ggplot2::aes(xintercept=x_min), colour="gray")
    g=g+ggplot2::geom_vline(data=boundData_df, ggplot2::aes(xintercept=x_max), colour="gray")

    g=g+ggplot2::geom_hline(data=boundData_df, ggplot2::aes(yintercept=y_min), colour="gray")
    g=g+ggplot2::geom_hline(data=boundData_df, ggplot2::aes(yintercept=y_max), colour="gray")
  }

  return(g)
}



makeSSRsurfaceDataset=function(logLocation, profile_likelihood=FALSE, numBins=NA, maxSSR=NA, ParameterNames=NA, ReparameterizationDef=NA, showInitialRange=FALSE, Residual_function=Residual_function_def, silent=FALSE){

  CGNM_result=NULL
  minSSR=NULL
  negative2LogLikelihood=NULL
  value=NULL
  boundValue=NULL
  initialX=NULL


  temp_list=prepSSRsurfaceData(logLocation, ParameterNames, ReparameterizationDef, numBins, Residual_function=Residual_function, silent=silent)
  X_nu=temp_list$X_matrix
  R_nu=temp_list$R_vec
  #  ResVar_nu=temp_list$ResVar_vec
  #  residual_variance=temp_list$residual_variance
  negative2LogLikelihood=temp_list$negative2LogLikelihood
  ParameterNames=temp_list$ParameterNames
  numBins=temp_list$numBins
  initialX=temp_list$initialX
  ReparameterizationDef=temp_list$ReparameterizationDef
  runName=temp_list$runName

  likelihoodSurfacePlot_df=data.frame()

  indexToUseForRefine_df=data.frame()

  for(j in seq(1,length(ParameterNames))){
    for(i in seq(1,numBins)){
      binUpper=quantile(X_nu[,j], probs = c((i-1)*(1/numBins),(i)*(1/numBins)))[2]
      binLower=quantile(X_nu[,j], probs = c((i-1)*(1/numBins),(i)*(1/numBins)))[1]
      indexToUse=((X_nu[,j]<=binUpper)&(X_nu[,j]>binLower))

      n2l_toUse=negative2LogLikelihood[indexToUse]
      R_toUse=R_nu[indexToUse]

      #plot_index=which(indexToUse&R_nu==min(R_toUse,na.rm = TRUE))
      plot_index=which(indexToUse&negative2LogLikelihood==min(n2l_toUse,na.rm = TRUE))

      if(length(plot_index)>0){
        plot_index=plot_index[1]
        indexToUseForRefine_df=rbind(indexToUseForRefine_df,data.frame(index= which(indexToUse)[order(R_toUse)[seq(1, min(length(R_nu),length(ParameterNames)+1))]], paraIndex=j))
        likelihoodSurfacePlot_df=rbind(likelihoodSurfacePlot_df, data.frame(value=X_nu[plot_index,j], negative2LogLikelihood=negative2LogLikelihood[plot_index], minSSR=min(R_nu[indexToUse],na.rm = TRUE), parameterName=ParameterNames[j]))

      }


    }
  }

  likelihoodSurfacePlot_df$runName=runName


  #  ggplot(likelihoodSurfacePlot_df, aes(x=value,y=minSSR))+geom_point()+facet_wrap(.~parameterName,scales = "free")+geom_smooth()+ylim(0.5,1)

  if(!is.na(maxSSR)){
    likelihoodSurfacePlot_df=subset(likelihoodSurfacePlot_df, minSSR<=maxSSR)
  }


  out=list(likelihoodSurfacePlot_df=likelihoodSurfacePlot_df, indexToUseForRefine_df=indexToUseForRefine_df, X_nu=X_nu,initialX=initialX, ParameterNames=ParameterNames, ReparameterizationDef=ReparameterizationDef)

  return(out)
}




makeInitialClusterForPLrefinment=function(logLocation, paraIndex, range=NA,numPerBin=10, profile_likelihood=FALSE, numBins=NA, maxSSR=NA, ParameterNames=NA, ReparameterizationDef=NA, showInitialRange=FALSE, Residual_function=Residual_function_def){

  CGNM_result=NULL
  minSSR=NULL
  negative2LogLikelihood=NULL
  value=NULL
  boundValue=NULL
  initialX=NULL


  temp_list=prepSSRsurfaceData(logLocation, ParameterNames, ReparameterizationDef, numBins, Residual_function=Residual_function)
  X_nu=temp_list$X_matrix
  R_nu=temp_list$R_vec
  #  residual_variance=temp_list$residual_variance
  ParameterNames=temp_list$ParameterNames
  numBins=temp_list$numBins
  initialX=temp_list$initialX
  ReparameterizationDef=temp_list$ReparameterizationDef


  if(!is.na(range)[1]){
    OK_index=X_nu[,paraIndex]<range[2]&X_nu[,paraIndex]>range[1]
  }

  X_nu=X_nu[OK_index,]
  R_nu=R_nu[OK_index]

  likelihoodSurfacePlot_df=data.frame()

  indexToUseForRefine_df=data.frame()

  j=paraIndex

  out_df=data.frame()
  for(i in seq(1,numBins)){

    binLower=(max(X_nu[,j])-min(X_nu[,j]))/numBins*(i-1)+min(X_nu[,j])
    binUpper=(max(X_nu[,j])-min(X_nu[,j]))/numBins*i+min(X_nu[,j])

    #    binUpper=quantile(X_nu[,j], probs = c((i-1)*(1/numBins),(i)*(1/numBins)))[2]
    #    binLower=quantile(X_nu[,j], probs = c((i-1)*(1/numBins),(i)*(1/numBins)))[1]
    indexToUse=((X_nu[,j]<=binUpper)&(X_nu[,j]>binLower))

    R_toUse=R_nu[indexToUse]

    plot_index=which(indexToUse&R_nu<=sort(R_toUse, na.last= TRUE)[numPerBin])

    tempDf=X_nu[plot_index,]

    best_index=which(indexToUse&R_nu==min(R_toUse))

    tempDf[,j]=X_nu[best_index,j]

    out_df=rbind(out_df, tempDf)

  }

  return(out_df)
}

#' @title plot_SSR_parameterValue
#' @description
#' Make SSR v.s. parameterValue plot of the accepted approximate minimizers found by the CGNM. Bars in the violin plots indicates the interquartile range.
#' @param CGNM_result (required input) \emph{A list} stores the computational result from Cluster_Gauss_Newton_method() function in CGNM package.
#' @param indicesToInclude (default: NA) \emph{A vector of integers} indices to include in the plot (if NA, use indices chosen by the acceptedIndices() function with default setting).
#' @param ParameterNames (default: NA) \emph{A vector of strings} the user can supply so that these names are used when making the plot. (Note if it set as NA or vector of incorrect length then the parameters are named as theta1, theta2, ... or as in ReparameterizationDef)
#' @param ReparameterizationDef (default: NA) \emph{A vector of strings} the user can supply definition of reparameterization where each string follows R syntax
#' @param showInitialRange (default: TRUE) \emph{TRUE or FALSE} if TRUE then the initial range appears in the plot.
#' @return \emph{A ggplot object} including the violin plot, interquartile range and median, minimum and maximum.
#' @examples
#'
#'model_analytic_function=function(x){
#'
#'  observation_time=c(0.1,0.2,0.4,0.6,1,2,3,6,12)
#'  Dose=1000
#'  F=1
#'
#'  ka=x[1]
#'  V1=x[2]
#'  CL_2=x[3]
#'  t=observation_time
#'
#'  Cp=ka*F*Dose/(V1*(ka-CL_2/V1))*(exp(-CL_2/V1*t)-exp(-ka*t))
#'
#'  log10(Cp)
#'}
#'
#' observation=log10(c(4.91, 8.65, 12.4, 18.7, 24.3, 24.5, 18.4, 4.66, 0.238))
#'
#' CGNM_result=Cluster_Gauss_Newton_method(
#' nonlinearFunction=model_analytic_function,
#' targetVector = observation,
#' initial_lowerRange = c(0.1,0.1,0.1), initial_upperRange =  c(10,10,10),
#' num_iter = 10, num_minimizersToFind = 100, saveLog = FALSE)
#'
#' plot_SSR_parameterValue(CGNM_result)
#' @export
#' @import ggplot2
plot_SSR_parameterValue=function(CGNM_result, indicesToInclude=NA, ParameterNames=NA, ReparameterizationDef=NA, showInitialRange=TRUE){

  Kind_iter=NULL
  SSR=NULL
  X_value=NULL
  cluster=NULL
  boundValue=NULL

  cutoff_pvalue=0.05
  numParametersIncluded=NA
  useAcceptedApproximateMinimizers=TRUE
  #   if(is.na(ParameterNames)[1]&!is.null(CGNM_result$runSetting$ParameterNames)){
  #     ParameterNames=CGNM_result$runSetting$ParameterNames
  # }
  #   if(is.na(ReparameterizationDef)[1]&!is.null(CGNM_result$runSetting$ReparameterizationDef)){
  #     ReparameterizationDef=CGNM_result$runSetting$ReparameterizationDef
  #   }


  freeParaValues=makeParaDistributionPlotDataFrame(CGNM_result, indicesToInclude, cutoff_pvalue, numParametersIncluded, ParameterNames, ReparameterizationDef, useAcceptedApproximateMinimizers)

  plotdata_df=subset(freeParaValues, Kind_iter=="Final Accepted")

  kmeans_result=optimal_kmeans(matrix(plotdata_df$X_value, ncol=dim(CGNM_result$X)[2]))

  plotdata_df$cluster=as.factor(kmeans_result$cluster)

  g=ggplot2::ggplot(plotdata_df, ggplot2::aes(y=SSR, x=X_value, colour=cluster))+ggplot2::geom_point(alpha=0.3)

  if(showInitialRange&&is.na(ParameterNames)[1]){
    ParameterNames=paste0("x",seq(1,dim(CGNM_result$X)[2]))
    g=g+ggplot2::geom_vline(data=data.frame(boundValue=CGNM_result$runSetting$initial_lowerRange, Name=ParameterNames), ggplot2::aes(xintercept=boundValue), colour="gray")
    g=g+ggplot2::geom_vline(data=data.frame(boundValue=CGNM_result$runSetting$initial_upperRange, Name=ParameterNames), ggplot2::aes(xintercept=boundValue), colour="gray")
  }

  g+ggplot2::facet_wrap(.~Name,scales="free")
}



#' @title plot_parameterValue_scatterPlots
#' @description
#' Make scatter plots of the accepted approximate minimizers found by the CGNM. Bars in the violin plots indicates the interquartile range.
#' @param CGNM_result (required input) \emph{A list} stores the computational result from Cluster_Gauss_Newton_method() function in CGNM package.
#' @param indicesToInclude (default: NA) \emph{A vector of integers} indices to include in the plot (if NA, use indices chosen by the acceptedIndices() function with default setting).
#' @return \emph{A ggplot object} including the violin plot, interquartile range and median, minimum and maximum.
#' @examples
#'
#'model_analytic_function=function(x){
#'
#'  observation_time=c(0.1,0.2,0.4,0.6,1,2,3,6,12)
#'  Dose=1000
#'  F=1
#'
#'  ka=x[1]
#'  V1=x[2]
#'  CL_2=x[3]
#'  t=observation_time
#'
#'  Cp=ka*F*Dose/(V1*(ka-CL_2/V1))*(exp(-CL_2/V1*t)-exp(-ka*t))
#'
#'  log10(Cp)
#'}
#'
#' observation=log10(c(4.91, 8.65, 12.4, 18.7, 24.3, 24.5, 18.4, 4.66, 0.238))
#'
#' CGNM_result=Cluster_Gauss_Newton_method(
#' nonlinearFunction=model_analytic_function,
#' targetVector = observation,
#' initial_lowerRange = c(0.1,0.1,0.1), initial_upperRange =  c(10,10,10),
#' num_iter = 10, num_minimizersToFind = 100, saveLog = FALSE)
#'
#' plot_parameterValue_scatterPlots(CGNM_result)
#' @export
#' @import ggplot2
plot_parameterValue_scatterPlots=function(CGNM_result, indicesToInclude=NA){

  cutoff_pvalue=0.05
  numParametersIncluded=NA
  ParameterNames=NA
  useAcceptedApproximateMinimizers=TRUE

  SSR=NULL
  X_value=NULL
  cluster=NULL
  Y_value=NULL
  cluster=NULL

  lastIter=dim(CGNM_result$residual_history)[2]

  plot_df=data.frame()

  if(!is.na(indicesToInclude[1])){
    useIndecies_b=seq(1,dim(CGNM_result$X)[1])%in%indicesToInclude
  }else{
    useIndecies_b=acceptedIndices_binary(CGNM_result,cutoff_pvalue, numParametersIncluded, useAcceptedApproximateMinimizers)
  }

  if(length(ParameterNames)==dim(CGNM_result$initialX)[2]){

  }else{
    ParameterNames=paste0("x_",seq(1,dim(CGNM_result$initialX)[2]))
  }


  kmeans_result=optimal_kmeans(cbind(CGNM_result$X[useIndecies_b,], plot_df$SSR))
  cluster_vec=as.factor(kmeans_result$cluster)

  for(i in seq(1,dim(CGNM_result$X)[2])){
    for(j in seq(1,dim(CGNM_result$X)[2])){
      plot_df=rbind(plot_df, data.frame(xName=ParameterNames[i],X_value=CGNM_result$X[useIndecies_b,i],yName=ParameterNames[j],Y_value=CGNM_result$X[useIndecies_b,j],cluster=cluster_vec))
    }
  }



  ggplot2::ggplot(plot_df, ggplot2::aes(y=Y_value, x=X_value, colour=cluster))+ggplot2::geom_point(alpha=0.3)+ggplot2::facet_wrap(xName~yName,scales="free")+
    ggplot2::theme(
      strip.background = element_blank(),
      strip.text.x = element_blank()
    )

}


#' @title bestApproximateMinimizers
#' @description
#' Returns the approximate minimizers with minimum SSR found by CGNM.
#' @param CGNM_result (required input) \emph{A list} stores the computational result from Cluster_Gauss_Newton_method() function in CGNM package.
#' @param numParameterSet (default 1) \emph{A natural number} number of parameter sets to output (chosen from the smallest SSR to numParameterSet-th smallest SSR) .
#' @param ParameterNames (default: NA) \emph{A vector of strings} the user can supply so that these names are used when making the plot. (Note if it set as NA or vector of incorrect length then the parameters are named as theta1, theta2, ... or as in ReparameterizationDef)
#' @param ReparameterizationDef (default: NA) \emph{A vector of strings} the user can supply definition of reparameterization where each string follows R syntax
#' @return \emph{A vector} a vector of accepted approximate minimizers with minimum SSR found by CGNM.
#' @examples
#'
#'model_analytic_function=function(x){
#'
#'  observation_time=c(0.1,0.2,0.4,0.6,1,2,3,6,12)
#'  Dose=1000
#'  F=1
#'
#'  ka=x[1]
#'  V1=x[2]
#'  CL_2=x[3]
#'  t=observation_time
#'
#'  Cp=ka*F*Dose/(V1*(ka-CL_2/V1))*(exp(-CL_2/V1*t)-exp(-ka*t))
#'
#'  log10(Cp)
#'}
#'
#' observation=log10(c(4.91, 8.65, 12.4, 18.7, 24.3, 24.5, 18.4, 4.66, 0.238))
#'
#' CGNM_result=Cluster_Gauss_Newton_method(
#' nonlinearFunction=model_analytic_function,
#' targetVector = observation,
#' initial_lowerRange = c(0.1,0.1,0.1), initial_upperRange =  c(10,10,10),
#' num_iter = 10, num_minimizersToFind = 100, saveLog = FALSE)
#'
#' bestApproximateMinimizers(CGNM_result,10)
#' @export
bestApproximateMinimizers=function(CGNM_result, numParameterSet=1,ParameterNames=NA, ReparameterizationDef=NA){
  SSRorder=order(CGNM_result$residual_history[,dim(CGNM_result$residual_history)[2]])

  out=CGNM_result$X[SSRorder,]

  out=data.frame(out)


  out=out[1:numParameterSet,]

  colnames(out)=paste0("x",seq(1,dim(CGNM_result$initialX)[2]))


  if(is.na(ParameterNames)[1]&!is.null(CGNM_result$runSetting$ParameterNames)){
    ParameterNames=CGNM_result$runSetting$ParameterNames
  }


  if(is.na(ReparameterizationDef)[1]&!is.null(CGNM_result$runSetting$ReparameterizationDef)){
    ReparameterizationDef=CGNM_result$runSetting$ReparameterizationDef
  }



  if(is.na(ParameterNames[1])|length(ParameterNames)!=length(ReparameterizationDef)){
    ParameterNames=ReparameterizationDef
  }

  if(is.na(ParameterNames[1])){
    ParameterNames=paste0("x",seq(1,dim(CGNM_result$initialX)[2]))
    ReparameterizationDef=ParameterNames
  }

  repara_finalX_df=data.frame(row.names = seq(1,dim(out)[1]))

  for(i in seq(1,length(ParameterNames))){
    repara_finalX_df[,ParameterNames[i]]=with(out, eval(parse(text=ReparameterizationDef[i])))
  }

  return(repara_finalX_df)


}

#' @title plot_goodnessOfFit
#' @description
#' Make goodness of fit plots to assess the model-fit and bias in residual distribution. The linear model is fit to the residual and plotted using geom_smooth(method=lm) in ggplot.\cr\cr
#' Explanation of the terminologies in terms of PBPK model fitting to the time-course drug concentration measurements:
#' \cr "independent variable" is time
#' \cr "dependent variable" is the concentration.
#' \cr "Residual" is the difference between the measured concentration and the model simulation with the parameter fond by the CGNM.
#' \cr "m" is number of observations
#' @param CGNM_result (required input) \emph{A list} stores the computational result from Cluster_Gauss_Newton_method() function in CGNM package.
#' @param plotType (default: 1) \emph{1,2 or 3}\cr specify the kind of goodness of fit plot to create
#' @param plotRank (default: c(1)) \emph{a vector of integers}\cr Specify which rank of the parameter to use for the goodness of fit plots. (e.g., if one wishes to use rank 1 to 100 then set it to be seq(1,100), or if one wish to use 88th rank parameters then set this as 88.)
#' @param independentVariableVector (default: NA) \emph{a vector of numerics of length m} \cr set independent variables that target values are associated with (e.g., time of the drug concentration measurement one is fitting PBPK model to) \cr(when this variable is set to NA, seq(1,m) will be used as independent variable when appropriate).
#' @param dependentVariableTypeVector (default: NA) \emph{a vector of text of length m} \cr when this variable is set (i.e., not NA) then the goodness of fit analyses is done for each variable type.  For example, if we are fitting the PBPK model to data with multiple dose arms, one can see the goodness of fit for each dose arm by specifying which dose group the observations are from.
#' @param absResidual (default: FALSE)  \emph{TRUE or FALSE} If TRUE plot absolute values of the residual.
#' @return \emph{A ggplot object} of the goodness of fit plot.
#' @examples
#'
#'model_analytic_function=function(x){
#'
#'  observation_time=c(0.1,0.2,0.4,0.6,1,2,3,6,12)
#'  Dose=1000
#'  F=1
#'
#'  ka=10^x[1]
#'  V1=10^x[2]
#'  CL_2=10^x[3]
#'  t=observation_time
#'
#'  Cp=ka*F*Dose/(V1*(ka-CL_2/V1))*(exp(-CL_2/V1*t)-exp(-ka*t))
#'
#'  log10(Cp)
#'}
#'
#' observation=log10(c(4.91, 8.65, 12.4, 18.7, 24.3, 24.5, 18.4, 4.66, 0.238))
#'
#' CGNM_result=Cluster_Gauss_Newton_method(
#' nonlinearFunction=model_analytic_function,
#' targetVector = observation,
#' initial_lowerRange = rep(0.01,3), initial_upperRange =  rep(100,3),
#' lowerBound=rep(0,3), ParameterNames = c("Ka","V1","CL"),
#' num_iter = 10, num_minimizersToFind = 100, saveLog = FALSE)
#'
#' plot_goodnessOfFit(CGNM_result)
#' plot_goodnessOfFit(CGNM_result,
#'      independentVariableVector=c(0.1,0.2,0.4,0.6,1,2,3,6,12))
#' @export
#' @import ggplot2

plot_goodnessOfFit=function(CGNM_result, plotType=1, plotRank=c(1), independentVariableVector=NA, dependentVariableTypeVector=NA, absResidual=FALSE){

  CGNM_result$runSetting$targetVector=as.numeric(CGNM_result$runSetting$targetVector)

  independent_variable=NULL
  lower=NULL
  upper=NULL
  target=NULL
  model_fit=NULL
  ind=NULL
  residual=NULL
  SSR=NULL
  SD=NULL


  independentVariableVector_in=independentVariableVector

  SSR_vec=CGNM_result$residual_history[,dim(CGNM_result$residual_history)[2]]
  indexToInclude=seq(1,length(SSR_vec))[rank(SSR_vec,na.last = TRUE, ties.method = "last")%in%plotRank]

  if(is.na(independentVariableVector[1])||length(independentVariableVector)!=length(CGNM_result$runSetting$targetVector)){
    independentVariableVector=seq(1, length(CGNM_result$runSetting$targetVector))
  }else{
    independentVariableVector=as.numeric(independentVariableVector)
  }

  if(is.na(dependentVariableTypeVector[1])||length(dependentVariableTypeVector)!=length(CGNM_result$runSetting$targetVector)){
    dependentVariableTypeVector=NA
  }else{
    dependentVariableTypeVector=as.factor(dependentVariableTypeVector)
  }

  residualPlot_df=data.frame()

  for(index in indexToInclude){
    residualPlot_df=rbind(residualPlot_df, data.frame(independent_variable=independentVariableVector,
                                                      residual=CGNM_result$Y[index,]-CGNM_result$runSetting$targetVector,
                                                      target=CGNM_result$runSetting$targetVector,
                                                      model_fit=CGNM_result$Y[index,],
                                                      dependent_variable_type=dependentVariableTypeVector,
                                                      square_residual=(CGNM_result$Y[index,]-CGNM_result$runSetting$targetVector)^2,
                                                      ind=index))
  }
  residualPlot_df$ind=as.factor(residualPlot_df$ind)


  if(absResidual){
    residualPlot_df$residual=abs(residualPlot_df$residual)
  }

  residualPlot_df=residualPlot_df[!is.na(residualPlot_df$residual),]

  withBootstrap=!is.null(CGNM_result$bootstrapY)

  median_vec=c()
  percentile5_vec=c()
  percentile95_vec=c()
  uncertaintyBound_df=data.frame()

  if(withBootstrap){
    for(i in seq(1, dim(CGNM_result$bootstrapY)[2])){
      tempQ=quantile(CGNM_result$bootstrapY[,i], prob=c(0.05,0.5,0.95), na.rm=TRUE)
      median_vec=c(median_vec, tempQ[2])
      percentile5_vec=c(percentile5_vec, tempQ[1])
      percentile95_vec=c(percentile95_vec, tempQ[3])
    }
    uncertaintyBound_df=data.frame(independent_variable=independentVariableVector, dependent_variable_type=dependentVariableTypeVector, median=median_vec, lower=percentile5_vec, upper=percentile95_vec, target=CGNM_result$runSetting$targetVector)
  }

  uncertaintyBound_df=uncertaintyBound_df[!is.na(uncertaintyBound_df$target),]

  if(plotType==1){ #dependent variable v.s. independent variable (e.g., concentration-time profile)
    if(withBootstrap){
      p<-ggplot2::ggplot(uncertaintyBound_df,ggplot2::aes(x=independent_variable, y=median))+ggplot2::geom_line(colour="blue")+ggplot2::geom_ribbon(ggplot2::aes(ymin =lower, ymax =upper, x=independent_variable), alpha=0.2,fill="blue")+ggplot2::geom_point(ggplot2::aes(x=independent_variable, y=target), colour="red")

    }else if(length(unique(residualPlot_df$ind))==1){
      p<-ggplot2::ggplot(residualPlot_df,ggplot2::aes(x=independent_variable, y=model_fit))+ggplot2::geom_line()+ggplot2::geom_point(ggplot2::aes(x=independent_variable, y=target), colour="red")

    }else{
      p<-ggplot2::ggplot(residualPlot_df,ggplot2::aes(x=independent_variable, y=model_fit, group=ind))+ggplot2::geom_line()+ggplot2::geom_point(ggplot2::aes(x=independent_variable, y=target), colour="red")

    }

    p=p+ggplot2::xlab("Independent Variable")
    p=p+ggplot2::ylab("Dependent Variable")

    if(!is.na(dependentVariableTypeVector[1])){
      p=p+ggplot2::facet_grid(.~dependent_variable_type, scales = "free")
    }
  }else if(plotType==2){ #residual v.s. dependent variable (e.g., residual-fitted model profile)
    p<-ggplot2::ggplot(residualPlot_df,ggplot2::aes(x=model_fit, y=residual))+ggplot2::geom_point()
    p=p+ggplot2::geom_smooth(method = lm)
    p=p+ggplot2::xlab("Dependent Variable (simulation)")
    if(absResidual){
      p=p+ggplot2::ylab("Residual (absolute value)")
    }else{
      p=p+ggplot2::ylab("Residual")
    }

    if(!is.na(dependentVariableTypeVector[1])){
      p=p+ggplot2::facet_grid(.~dependent_variable_type, scales = "free")
    }

    p=p+geom_hline(yintercept = 0)
  }else if(plotType==4){ #residual v.s. dependent variable (e.g., residual-fitted model profile)
    p<-ggplot2::ggplot(residualPlot_df,ggplot2::aes(x=target, y=residual))+ggplot2::geom_point()
    p=p+ggplot2::geom_smooth(method = lm)
    p=p+ggplot2::xlab("Dependent Variable (target)")
    if(absResidual){
      p=p+ggplot2::ylab("Residual (absolute value)")
    }else{
      p=p+ggplot2::ylab("Residual")
    }
    if(!is.na(dependentVariableTypeVector[1])){
      p=p+ggplot2::facet_grid(.~dependent_variable_type, scales = "free")
    }

    p=p+geom_hline(yintercept = 0)
  }else if(plotType==3){ #residual v.s. independent variable (e.g., residual-time profile)

    if(is.na(independentVariableVector_in[1])){
      p<-ggplot2::ggplot(residualPlot_df,ggplot2::aes(x=residual))+ggplot2::geom_histogram()
      if(absResidual){
        p=p+ggplot2::xlab("Residual (absolute value)")
      }else{
        p=p+ggplot2::xlab("Residual")
      }

      if(!is.na(dependentVariableTypeVector[1])){
        p=p+ggplot2::facet_grid(.~dependent_variable_type, scales = "free")
      }
    }else{
      p<-ggplot2::ggplot(residualPlot_df,ggplot2::aes(x=independent_variable, y=residual))+ggplot2::geom_point()
      p=p+ggplot2::geom_smooth(method = lm)
      p=p+ggplot2::xlab("Independent Variable")
      if(absResidual){
        p=p+ggplot2::ylab("Residual (absolute value)")
      }else{
        p=p+ggplot2::ylab("Residual")
      }
      if(!is.na(dependentVariableTypeVector[1])){
        p=p+ggplot2::facet_grid(.~dependent_variable_type, scales = "free")
      }

      p=p+geom_hline(yintercept = 0)
    }

  }else{
    p=ggplot2::ggplot(residualPlot_df)
  }


  if(plotType==1&&(length(unique(dependentVariableTypeVector ))>1&length(plotRank)==1)){
    SSR_byVariable=aggregate(residualPlot_df$square_residual, by=list(dependent_variable_type=residualPlot_df$dependent_variable_type), FUN=sum)
    SSR_byVariable$SSR=formatC(SSR_byVariable$x, format = "g", digits = 3)

    p=p+
      geom_text(
        size    = 5,
        data    = SSR_byVariable,
        mapping = aes(x = Inf, y = Inf, label = paste0("SSR=",SSR)),
        hjust   = 1.05,
        vjust   = 1.5,
        inherit.aes = FALSE
      )
  }else if ((length(unique(dependentVariableTypeVector ))>1&length(plotRank)==1)){

    SD_byVariable=aggregate(residualPlot_df$residual, by=list(dependent_variable_type=residualPlot_df$dependent_variable_type), FUN=sd)
    SD_byVariable$SD=formatC(SD_byVariable$x, format = "g", digits = 3)

    p=p+
      geom_text(
        size    = 5,
        data    = SD_byVariable,
        mapping = aes(x = Inf, y = Inf, label = paste0("StandardDiv.=",SD)),
        hjust   = 1.05,
        vjust   = 1.5,
        inherit.aes = FALSE
      )
  }

  p
}




#' @title table_parameterSummary
#' @description
#' Make summary table of the approximate local minimizers found by CGNM.  If bootstrap analysis result is available, relative standard error (RSE: standard deviation/mean) will also be included in the table.
#' @param CGNM_result (required input) \emph{A list} stores the computational result from Cluster_Gauss_Newton_method() function in CGNM package.
#' @param indicesToInclude (default: NA) \emph{A vector of integers} indices to include in the plot (if NA, use indices chosen by the acceptedIndices() function with default setting).
#' @param ParameterNames (default: NA) \emph{A vector of strings} the user can supply so that these names are used when making the plot. (Note if it set as NA or vector of incorrect length then the parameters are named as theta1, theta2, ... or as in ReparameterizationDef)
#' @param ReparameterizationDef (default: NA) \emph{A vector of strings} the user can supply definition of reparameterization where each string follows R syntax.
#' @return \emph{A ggplot object} including the violin plot, interquartile range and median, minimum and maximum.
#' @examples
#'
#'model_analytic_function=function(x){
#'
#'  observation_time=c(0.1,0.2,0.4,0.6,1,2,3,6,12)
#'  Dose=1000
#'  F=1
#'
#'  ka=x[1]
#'  V1=x[2]
#'  CL_2=x[3]
#'  t=observation_time
#'
#'  Cp=ka*F*Dose/(V1*(ka-CL_2/V1))*(exp(-CL_2/V1*t)-exp(-ka*t))
#'
#'  log10(Cp)
#'}
#'
#' observation=log10(c(4.91, 8.65, 12.4, 18.7, 24.3, 24.5, 18.4, 4.66, 0.238))
#'
#'
#'
#' CGNM_result=Cluster_Gauss_Newton_method(
#' nonlinearFunction=model_analytic_function,
#' targetVector = observation,
#' initial_lowerRange = rep(0.01,3), initial_upperRange =  rep(100,3),
#' lowerBound=rep(0,3), ParameterNames = c("Ka","V1","CL"),
#' num_iter = 10, num_minimizersToFind = 100, saveLog = FALSE)
#'
#' table_parameterSummary(CGNM_result)
#' table_parameterSummary(CGNM_result,
#'      ReparameterizationDef=c("log10(Ka)","log10(V1)","log10(CL)"))
#'
#' @export
#' @import ggplot2

table_parameterSummary=function(CGNM_result, indicesToInclude=NA, ParameterNames=NA, ReparameterizationDef=NA){

  #   if(is.na(ParameterNames)[1]&!is.null(CGNM_result$runSetting$ParameterNames)){
  #
  #     ParameterNames=CGNM_result$runSetting$ParameterNames
  #
  # }
  #
  # if(is.na(ReparameterizationDef)[1]&!is.null(CGNM_result$runSetting$ReparameterizationDef)){
  #   ReparameterizationDef=CGNM_result$runSetting$ReparameterizationDef
  # }

  Name=NULL
  Kind_iter=NULL

  cutoff_pvalue=0.05
  numParametersIncluded=NA
  useAcceptedApproximateMinimizers=TRUE


  freeParaValues=makeParaDistributionPlotDataFrame(CGNM_result, indicesToInclude, cutoff_pvalue, numParametersIncluded, ParameterNames, ReparameterizationDef, useAcceptedApproximateMinimizers)
  summaryData_df=data.frame()
  variableNames=unique(freeParaValues$Name)


  if(!is.null(CGNM_result$bootstrapX)){

    for(vName in variableNames){

      tempVec=quantile(subset(freeParaValues, Name==vName&Kind_iter=="Bootstrap")$X_value, probs = c(0,0.25,0.5,0.75,1), na.rm = TRUE)

      bootstrapDS=subset(freeParaValues, Name==vName&Kind_iter=="Bootstrap")$X_value
      tempVec=c(tempVec, sd(bootstrapDS,na.rm = TRUE)/abs(mean(bootstrapDS,na.rm = TRUE))*100)

      summaryData_df=rbind(summaryData_df,as.numeric(tempVec))
    }

    colnames(summaryData_df)=c("CGNM Bootstrap: Minimum","25 percentile", "Median","75 percentile", "Maximum", "RSE (%)")

  }else{


    for(vName in variableNames){

      tempVec=quantile(subset(freeParaValues, Name==vName&Kind_iter=="Final Accepted")$X_value, probs = c(0,0.25,0.5,0.75,1), na.rm = TRUE)

      summaryData_df=rbind(summaryData_df,as.numeric(tempVec))
    }

    colnames(summaryData_df)=c("CGNM: Minimum","25 percentile", "Median","75 percentile", "Maximum")

  }

  rownames(summaryData_df)=variableNames


  return(summaryData_df)
}



#' @title suggestInitialLowerRange
#' @description
#' Suggest initial lower range based on the profile likelihood. The user can re-run CGNM with this suggested initial range so that to improve the convergence.
#' @param logLocation (required input) \emph{A string or a list of strings} of folder directory where CGNM computation log files exist.
#' @param alpha (default: 0.25) \emph{a number between 0 and 1} level of significance used to derive the confidence interval.
#' @param numBins (default: NA) \emph{A positive integer} SSR surface is plotted by finding the minimum SSR given one of the parameters is fixed and then repeat this for various values.  numBins specifies the number of different parameter values to fix for each parameter. (if set NA the number of bins are set as num_minimizersToFind/10)
#' @return \emph{A numerical vector} of suggested initial lower range based on profile likelihood.
#' @examples
#'\dontrun{
#'model_analytic_function=function(x){
#'
#'  observation_time=c(0.1,0.2,0.4,0.6,1,2,3,6,12)
#'  Dose=1000
#'  F=1
#'
#'  ka=x[1]
#'  V1=x[2]
#'  CL_2=x[3]
#'  t=observation_time
#'
#'  Cp=ka*F*Dose/(V1*(ka-CL_2/V1))*(exp(-CL_2/V1*t)-exp(-ka*t))
#'
#'  log10(Cp)
#'}
#'
#' observation=log10(c(4.91, 8.65, 12.4, 18.7, 24.3, 24.5, 18.4, 4.66, 0.238))
#'
#' CGNM_result=Cluster_Gauss_Newton_method(
#' nonlinearFunction=model_analytic_function,
#' targetVector = observation,
#' initial_lowerRange = c(0.1,0.1,0.1), initial_upperRange =  c(10,10,10),
#' num_iter = 10, num_minimizersToFind = 100, saveLog=TRUE)
#'
#' suggestInitialLowerRange("CGNM_log")
#' }
#' @export
suggestInitialLowerRange=function(logLocation, alpha=0.25, numBins=NA){
  Residual_function=Residual_function_def
  ParameterNames=NA
  ReparameterizationDef=NA
  CGNM_result=NULL

  boundValue=NULL
  individual=NULL
  label=NULL
  minSSR=NULL
  negative2LogLikelihood=NULL
  newvalue=NULL

  parameterName=NULL

  reparaXinit=NULL
  value=NULL


  if(typeof(logLocation)=="character"){

    DirectoryName_vec=c(logLocation)

  }else if(typeof(logLocation)=="list"){
    if(logLocation$runSetting$runName==""){
      DirectoryName_vec=c("CGNM_log","CGNM_log_bootstrap")
    }else{
      DirectoryName_vec=c(paste0(logLocation$runSetting$runName,"_CGNM_log"),paste0(logLocation$runSetting$runName,"_CGNM_log_bootstrap"))
    }

  }else{
    warning("DirectoryName need to be either the CGNM_result object or a string")
  }

  load(paste0(DirectoryName_vec[1],"/iteration_1.RDATA"))


  data=makeSSRsurfaceDataset(logLocation, TRUE, numBins, NA, ParameterNames, ReparameterizationDef, FALSE,Residual_function=Residual_function)


  likelihoodSurfacePlot_df=data$likelihoodSurfacePlot_df

  minSSR=min(likelihoodSurfacePlot_df$negative2LogLikelihood)

  likelihoodSurfacePlot_df$belowSignificance=(likelihoodSurfacePlot_df$negative2LogLikelihood-qchisq(1-alpha,1)-minSSR<0)

  paraKind=unique(likelihoodSurfacePlot_df$parameterName)

  upperBound_vec=c()
  lowerBound_vec=c()

  for(para_nu in paraKind){
    dataframe_nu=subset(likelihoodSurfacePlot_df,parameterName==para_nu)

    upperBoundIndex=which(max(dataframe_nu$value[dataframe_nu$belowSignificance])==dataframe_nu$value)
    lowerBoundIndex=which(min(dataframe_nu$value[dataframe_nu$belowSignificance])==dataframe_nu$value)

    theoreticalLB=CGNM_result$runSetting$lowerBound[CGNM_result$runSetting$ParameterNames==para_nu]
    theoreticalUB=CGNM_result$runSetting$upperBound[CGNM_result$runSetting$ParameterNames==para_nu]

    if(upperBoundIndex==length(dataframe_nu$value)){
      if(!is.na(theoreticalUB)){
        upperBound=(dataframe_nu$value[upperBoundIndex]+theoreticalUB)/2

      }else{
        upperBound=dataframe_nu$value[upperBoundIndex]*10

      }
    }else{
      upperBound=(dataframe_nu$value[upperBoundIndex]+dataframe_nu$value[upperBoundIndex+1])/2
    }

    if(lowerBoundIndex==1){
      if(!is.na(theoreticalLB)){
        lowerBound_suggest=(dataframe_nu$value[lowerBoundIndex]+theoreticalLB)/2

      }else{
        lowerBound_suggest=10^floor(log10(dataframe_nu$value[lowerBoundIndex]))/10
      }

    }else if(lowerBoundIndex>2){
      lowerBound_suggest=dataframe_nu$value[lowerBoundIndex-2]
    }else{

      if(!is.na(theoreticalLB)){

        lowerBound=(dataframe_nu$value[lowerBoundIndex]+dataframe_nu$value[lowerBoundIndex-1])/2
        width=upperBound-lowerBound

        if((lowerBound-width)>theoreticalLB){
          lowerBound_suggest=(lowerBound-width)

        }else{
          lowerBound_suggest=(lowerBound+theoreticalLB)/2
        }

      }else{
        lowerBound=(dataframe_nu$value[lowerBoundIndex]+dataframe_nu$value[lowerBoundIndex-1])/2
        lowerBound_suggest=10^floor(log10(lowerBound))/10

      }
    }

    upperBound_vec=c(upperBound_vec,upperBound)
    lowerBound_vec=c(lowerBound_vec,lowerBound_suggest)
  }

  return(lowerBound_vec)
}


#' @title suggestInitialUpperRange
#' @description
#' Suggest initial upper range based on the profile likelihood. The user can re-run CGNM with this suggested initial range so that to improve the convergence.
#' @param logLocation (required input) \emph{A string or a list of strings} of folder directory where CGNM computation log files exist.
#' @param alpha (default: 0.25) \emph{a number between 0 and 1} level of significance used to derive the confidence interval.
#' @param numBins (default: NA) \emph{A positive integer} SSR surface is plotted by finding the minimum SSR given one of the parameters is fixed and then repeat this for various values.  numBins specifies the number of different parameter values to fix for each parameter. (if set NA the number of bins are set as num_minimizersToFind/10)
#' @return \emph{A numerical vector} of suggested initial upper range based on profile likelihood.
#' @examples
#'\dontrun{
#'model_analytic_function=function(x){
#'
#'  observation_time=c(0.1,0.2,0.4,0.6,1,2,3,6,12)
#'  Dose=1000
#'  F=1
#'
#'  ka=x[1]
#'  V1=x[2]
#'  CL_2=x[3]
#'  t=observation_time
#'
#'  Cp=ka*F*Dose/(V1*(ka-CL_2/V1))*(exp(-CL_2/V1*t)-exp(-ka*t))
#'
#'  log10(Cp)
#'}
#'
#' observation=log10(c(4.91, 8.65, 12.4, 18.7, 24.3, 24.5, 18.4, 4.66, 0.238))
#'
#' CGNM_result=Cluster_Gauss_Newton_method(
#' nonlinearFunction=model_analytic_function,
#' targetVector = observation,
#' initial_lowerRange = c(0.1,0.1,0.1), initial_upperRange =  c(10,10,10),
#' num_iter = 10, num_minimizersToFind = 100, saveLog=TRUE)
#'
#' suggestInitialLowerRange("CGNM_log")
#' }
#' @export
suggestInitialUpperRange=function(logLocation, alpha=0.25, numBins=NA){
  Residual_function=Residual_function_def
  ParameterNames=NA
  ReparameterizationDef=NA
  CGNM_result=NULL

  boundValue=NULL
  individual=NULL
  label=NULL
  minSSR=NULL
  negative2LogLikelihood=NULL
  newvalue=NULL

  parameterName=NULL

  reparaXinit=NULL
  value=NULL


  if(typeof(logLocation)=="character"){

    DirectoryName_vec=c(logLocation)

  }else if(typeof(logLocation)=="list"){
    if(logLocation$runSetting$runName==""){
      DirectoryName_vec=c("CGNM_log","CGNM_log_bootstrap")
    }else{
      DirectoryName_vec=c(paste0(logLocation$runSetting$runName,"_CGNM_log"),paste0(logLocation$runSetting$runName,"_CGNM_log_bootstrap"))
    }

  }else{
    warning("DirectoryName need to be either the CGNM_result object or a string")
  }

  load(paste0(DirectoryName_vec[1],"/iteration_1.RDATA"))


  data=makeSSRsurfaceDataset(logLocation, TRUE, numBins, NA, ParameterNames, ReparameterizationDef, FALSE, Residual_function=Residual_function)


  likelihoodSurfacePlot_df=data$likelihoodSurfacePlot_df

  minSSR=min(likelihoodSurfacePlot_df$negative2LogLikelihood)

  likelihoodSurfacePlot_df$belowSignificance=(likelihoodSurfacePlot_df$negative2LogLikelihood-qchisq(1-alpha,1)-minSSR<0)

  paraKind=unique(likelihoodSurfacePlot_df$parameterName)

  upperBound_vec=c()
  lowerBound_vec=c()

  for(para_nu in paraKind){
    dataframe_nu=subset(likelihoodSurfacePlot_df,parameterName==para_nu)

    upperBoundIndex=which(max(dataframe_nu$value[dataframe_nu$belowSignificance])==dataframe_nu$value)
    lowerBoundIndex=which(min(dataframe_nu$value[dataframe_nu$belowSignificance])==dataframe_nu$value)

    theoreticalLB=CGNM_result$runSetting$lowerBound[CGNM_result$runSetting$ParameterNames==para_nu]
    theoreticalUB=CGNM_result$runSetting$upperBound[CGNM_result$runSetting$ParameterNames==para_nu]

    if(lowerBoundIndex==1){
      if(!is.na(theoreticalLB)){
        lowerBound=(dataframe_nu$value[lowerBoundIndex]+theoreticalLB)/2

      }else{
        lowerBound=dataframe_nu$value[lowerBoundIndex]*10

      }
    }else{
      lowerBound=(dataframe_nu$value[lowerBoundIndex]+dataframe_nu$value[lowerBoundIndex+1])/2
    }



    if(upperBoundIndex==length(dataframe_nu$value)){
      if(!is.na(theoreticalUB)){
        upperBound_suggest=(dataframe_nu$value[upperBoundIndex]+theoreticalUB)/2

      }else{
        upperBound_suggest=10^ceiling(log10(dataframe_nu$value[upperBoundIndex]))*10
      }


    }else if((upperBoundIndex+2)<=length(dataframe_nu$value)){

      upperBound_suggest=dataframe_nu$value[upperBoundIndex+2]

    }else{

      if(!is.na(theoreticalUB)){


        upperBound=(dataframe_nu$value[upperBoundIndex]+dataframe_nu$value[upperBoundIndex-1])/2
        width=upperBound-lowerBound

        if((upperBound+width)<theoreticalUB){
          upperBound_suggest=(upperBound+width)

        }else{
          upperBound_suggest=(upperBound+theoreticalUB)/2

        }

      }else{
        upperBound=(dataframe_nu$value[upperBoundIndex]+dataframe_nu$value[upperBoundIndex-1])/2
        upperBound_suggest=10^ceiling(log10(upperBound))*10

      }
    }

    upperBound_vec=c(upperBound_vec,upperBound_suggest)
  }

  return(upperBound_vec)
}




#' @title plot_simulationWithCI
#' @description
#' Plot model simulation where the various parameter combinations are provided and conduct simulations and then the confidence interval (or more like a confidence region) is plotted.
#' @param simulationFunction (required input) \emph{A function} that maps the parameter vector to the simulation.
#' @param parameter_matrix (required input) \emph{A matrix of numbers} where each row contains the parameter combination that will be used for the simulations.
#' @param independentVariableVector (default: NA) \emph{A vector of numbers} that represents the independent variables of each points of the simulation (e.g., observation time) where used for the values of x-axis when plotting. If set at NA then sequence of 1,2,3,... will be used.
#' @param dependentVariableTypeVector (default: NA) \emph{A vector of strings} specify the kind of variable the simulationFunction simulate out. (i.e., if it simulate both PK and PD then indicate which simulation output is PK and which is PD).
#' @param confidenceLevels (default: c(25,75)) \emph{A vector of two numbers between 0 and 1} set the confidence interval that will be used for the plot.  Default is inter-quartile range.
#' @param observationVector (default: NA) \emph{A vector of numbers} used when wishing to overlay the plot of observations to the simulation.
#' @param observationIndpendentVariableVector (default: NA) \emph{A vector of numbers} used when wishing to overlay the plot of observations to the simulation.
#' @param observationDependentVariableTypeVector (default: NA) \emph{A vector of numbers} used when wishing to overlay the plot of observations to the simulation.
#' @param overLay (default: FALSE) \emph{TRUE or FALSE} if TRUE all variable types are overlayed in one plot, if FALSE the plot will be faceted by the variable type.
#' @return \emph{A list} including ggplot object ($plot) and simulated data matrix ($plotData_matrix).
#' @examples
#'\dontrun{
#'model_analytic_function=function(x){
#'
#'  observation_time=c(0.1,0.2,0.4,0.6,1,2,3,6,12)
#'  Dose=1000
#'  F=1
#'
#'  ka=x[1]
#'  V1=x[2]
#'  CL_2=x[3]
#'  t=observation_time
#'
#'  Cp=ka*F*Dose/(V1*(ka-CL_2/V1))*(exp(-CL_2/V1*t)-exp(-ka*t))
#'
#'  (Cp)
#'}
#'
#' observation=(c(4.91, 8.65, 12.4, 18.7, 24.3, 24.5, 18.4, 4.66, 0.238))
#'
#' CGNM_result=Cluster_Gauss_Newton_method(
#' nonlinearFunction=model_analytic_function,
#' targetVector = observation, num_iteration = 10, num_minimizersToFind = 100,
#' initial_lowerRange = c(0.1,0.1,0.1), initial_upperRange =  c(10,10,10),
#' lowerBound=rep(0,3), ParameterNames=c("Ka","V1","CL_2"), saveLog = FALSE)
#'
#' CGNM_bootstrap=Cluster_Gauss_Newton_Bootstrap_method(CGNM_result,
#'      nonlinearFunction=model_analytic_function, num_bootstrapSample=100)
#'
#' plot_simulationWithCI(model_analytic_function, as.matrix(CGNM_result$bootstrapTheta),
#' independentVariableVector=observation_time, observationVector=observation)
#' }
#' @export
#' @import ggplot2
plot_simulationWithCI=function(simulationFunction, parameter_matrix,  independentVariableVector=NA, dependentVariableTypeVector=NA, confidenceLevels=c(0.25,0.75), observationVector=NA, observationIndpendentVariableVector=NA, observationDependentVariableTypeVector=NA,overLay=FALSE){

  lengthSimulation=length(simulationFunction(as.numeric(parameter_matrix[1,])))

  plotData_matrix=matrix(NA, nrow = dim(parameter_matrix)[1], ncol = lengthSimulation)



  for(i in seq(1,dim(parameter_matrix)[1] )){
    plotData_matrix[i,]=simulationFunction(as.numeric(parameter_matrix[i,]))#, independentVariable=independentVariableVector, dependentVariableType=dependentVariableTypeVector))
  }

  g=plot_simulationMatrixWithCI(plotData_matrix,  independentVariableVector=independentVariableVector, dependentVariableTypeVector=dependentVariableTypeVector, confidenceLevels=confidenceLevels, observationVector=observationVector, observationIndpendentVariableVector=observationIndpendentVariableVector, observationDependentVariableTypeVector=observationDependentVariableTypeVector, overLay)


 # out=list(plot=g, plotData_matrix=plotData_matrix)
  return(g)
}

# plot_simulationWithCI=function(simulationFunction, parameter_matrix,  independentVariableVector=NA, dependentVariableTypeVector=NA, confidenceLevels=c(0.25,0.75), observationVector=NA, observationIndpendentVariableVector=NA, observationDependentVariableTypeVector=NA){
#
#   CGNM_result=NULL
#   independentVariable=NULL
#   dependentVariableType=NULL
#   lower_percentile=NULL
#   upper_percentile=NULL
#   observation=NULL
#
#   lengthSimulation=length(simulationFunction(as.numeric(parameter_matrix[1,])))
#
#   ValidIndependentVariableInput=TRUE
#
#   if(is.na(independentVariableVector[1])){
#     ValidIndependentVariableInput=FALSE
#     warning("independentVariableVector was not provided so replaced with seq(1,length of simulation)")
#
#   }else if(length(independentVariableVector)!=lengthSimulation){
#     ValidIndependentVariableInput=FALSE
#     warning("length of independentVariableVector was not the same length as the output of the simulationFunction so replaced with seq(1,length of output of simulation)")
#
#   }
#
#   if(!ValidIndependentVariableInput){
#     independentVariableVector=seq(1,lengthSimulation)
#   }
#
#   ValidobservationVectorInput=TRUE
#
#
#   if(is.na(observationVector[1])){
#     ValidobservationVectorInput=FALSE
#   }
#
#
#   if(ValidobservationVectorInput&is.na(observationIndpendentVariableVector[1])){
#     observationIndpendentVariableVector=independentVariableVector
#     warning("since observationIndpendentVariableVector is not provided replace it with independentVariableVector")
#
#   }
#
#   if(length(observationVector)!=length(observationIndpendentVariableVector)){
#     ValidobservationVectorInput=FALSE
#     warning("observationVector and observationIndpendentVariableVector need to be the same length hence observations will not be overlayed")
#   }
#
#   validDependentVariableTypeVectorInput=TRUE
#
#   if(is.na(dependentVariableTypeVector[1])){
#     validDependentVariableTypeVectorInput=FALSE
#   }else if(length(dependentVariableTypeVector)!=lengthSimulation){
#     dependentVariableTypeVector=NA
#     validDependentVariableTypeVectorInput=FALSE
#
#     warning("dependentVariableTypeVector was not the same length as the output of the simulationFunction so replaced with NA")
#   }
#
#
#
#   plot_df=data.frame()
#
#
#   for(i in seq(1,dim(parameter_matrix)[1] )){
#     plot_df=rbind(plot_df, data.frame(simulation=simulationFunction(as.numeric(parameter_matrix[i,])), independentVariable=independentVariableVector, dependentVariableType=dependentVariableTypeVector))
#   }
#
#   median_vec=c()
#   lower_percentile_vec=c()
#   upper_percentile_vec=c()
#   kind_df=unique(plot_df[,c("independentVariable", "dependentVariableType")])
#
#   for(i in seq(1,dim(kind_df)[1])){
#     kind=kind_df[i,]
#     if(validDependentVariableTypeVectorInput){
#       now_data_df=subset(plot_df, independentVariable==kind$independentVariable&dependentVariableType==kind$dependentVariableType )
#
#     }else{
#       now_data_df=subset(plot_df, independentVariable==kind$independentVariable )
#
#     }
#
#     nowQuantile=quantile(now_data_df$simulation, probs = c(confidenceLevels[1],0.5,confidenceLevels[2]), na.rm = TRUE)
#
#     lower_percentile_vec=c(lower_percentile_vec,nowQuantile[1])
#     median_vec=c(median_vec,nowQuantile[2])
#     upper_percentile_vec=c(upper_percentile_vec,nowQuantile[3])
#   }
#
#   plot_CI_df=kind_df
#   plot_CI_df$lower_percentile=as.numeric(lower_percentile_vec)
#   plot_CI_df$upper_percentile=as.numeric(upper_percentile_vec)
#   plot_CI_df$median=as.numeric(median_vec)
#
#
#   g=ggplot2::ggplot(plot_CI_df, aes(x=independentVariable , y=median ))+ggplot2::geom_line()+ggplot2::geom_ribbon(aes(ymin=lower_percentile, ymax=upper_percentile), alpha=0.2)+ ggplot2::labs(caption = paste0("solid line is the median of the model prediction and shaded area is its confidence interval of ",confidenceLevels[1]*100,"-",confidenceLevels[2]*100," percentile"))
#
#   if(validDependentVariableTypeVectorInput){
#     g=g+ggplot2::facet_wrap(.~dependentVariableType)
#   }
#   if(ValidobservationVectorInput){
#     g=g+ggplot2::geom_point(data=data.frame(independentVariable=observationIndpendentVariableVector, observation=observationVector, dependentVariableType=observationDependentVariableTypeVector), colour="red", aes(x=independentVariable,y=observation), inherit.aes = FALSE)
#   }
#
#   g=g+ggplot2::ylab("Dependent variable")
#
#   return(g)
# }





#' @title plot_simulationMatrixWithCI
#' @description
#' Plot simulation that are provided to plot confidence interval (or more like a confidence region).
#' @param simulationMatrix (required input) \emph{A matrix of numbers} where each row contains the simulated values that will be plotted.
#' @param independentVariableVector (default: NA) \emph{A vector of numbers} that represents the independent variables of each points of the simulation (e.g., observation time) where used for the values of x-axis when plotting. If set at NA then sequence of 1,2,3,... will be used.
#' @param dependentVariableTypeVector (default: NA) \emph{A vector of strings} specify the kind of variable the simulation values are. (i.e., if it simulate both PK and PD then indicate which simulation value is PK and which is PD).
#' @param confidenceLevels (default: c(25,75)) \emph{A vector of two numbers between 0 and 1} set the confidence interval that will be used for the plot.  Default is inter-quartile range.
#' @param observationVector (default: NA) \emph{A vector of numbers} used when wishing to overlay the plot of observations to the simulation.
#' @param observationIndpendentVariableVector (default: NA) \emph{A vector of numbers} used when wishing to overlay the plot of observations to the simulation.
#' @param observationDependentVariableTypeVector (default: NA) \emph{A vector of numbers} used when wishing to overlay the plot of observations to the simulation.
#' @param overLay (default: FALSE) \emph{TRUE or FALSE} if TRUE all variable types are overlayed in one plot, if FALSE the plot will be faceted by the variable type.
#' @return \emph{A ggplot object}.
#' @examples
#'\dontrun{
#'model_analytic_function=function(x){
#'
#'  observation_time=c(0.1,0.2,0.4,0.6,1,2,3,6,12)
#'  Dose=1000
#'  F=1
#'
#'  ka=x[1]
#'  V1=x[2]
#'  CL_2=x[3]
#'  t=observation_time
#'
#'  Cp=ka*F*Dose/(V1*(ka-CL_2/V1))*(exp(-CL_2/V1*t)-exp(-ka*t))
#'
#'  (Cp)
#'}
#'
#' observation=(c(4.91, 8.65, 12.4, 18.7, 24.3, 24.5, 18.4, 4.66, 0.238))
#'
#' CGNM_result=Cluster_Gauss_Newton_method(
#' nonlinearFunction=model_analytic_function,
#' targetVector = observation, num_iteration = 10, num_minimizersToFind = 100,
#' initial_lowerRange = c(0.1,0.1,0.1), initial_upperRange =  c(10,10,10),
#' lowerBound=rep(0,3), ParameterNames=c("Ka","V1","CL_2"), saveLog = FALSE)
#'
#' CGNM_bootstrap=Cluster_Gauss_Newton_Bootstrap_method(CGNM_result,
#'      nonlinearFunction=model_analytic_function, num_bootstrapSample=100)
#'
#'
#' plot_simulationMatrixWithCI(CGNM_result$bootstrapY,
#' independentVariableVector=observation_time, observationVector=observation)
#' }
#' @export
#' @import ggplot2

plot_simulationMatrixWithCI=function(simulationMatrix,  independentVariableVector=NA, dependentVariableTypeVector=NA, confidenceLevels=c(0.25,0.75), observationVector=NA, observationIndpendentVariableVector=NA, observationDependentVariableTypeVector=NA, overLay=FALSE){

  CGNM_result=NULL
  independentVariable=NULL
  dependentVariableType=NULL
  lower_percentile=NULL
  upper_percentile=NULL
  observation=NULL

  lengthSimulation=dim(simulationMatrix)[2]

  ValidIndependentVariableInput=TRUE

  if(is.na(independentVariableVector[1])){
    ValidIndependentVariableInput=FALSE
    warning("independentVariableVector was not provided so replaced with seq(1,length of simulation)")

  }else if(length(independentVariableVector)!=lengthSimulation){
    ValidIndependentVariableInput=FALSE
    warning("length of independentVariableVector was not the same length as the output of the simulationFunction so replaced with seq(1,length of output of simulation)")

  }

  if(!ValidIndependentVariableInput){
    independentVariableVector=seq(1,lengthSimulation)
  }

  ValidobservationVectorInput=TRUE


  if(is.na(observationVector[1])){
    ValidobservationVectorInput=FALSE
  }


  if(ValidobservationVectorInput&is.na(observationIndpendentVariableVector[1])){
    observationIndpendentVariableVector=independentVariableVector
    warning("since observationIndpendentVariableVector is not provided replace it with independentVariableVector")

  }

  if(length(observationVector)!=length(observationIndpendentVariableVector)){
    ValidobservationVectorInput=FALSE
    warning("observationVector and observationIndpendentVariableVector need to be the same length hence observations will not be overlayed")
  }

  validDependentVariableTypeVectorInput=TRUE

  if(is.na(dependentVariableTypeVector[1])){
    validDependentVariableTypeVectorInput=FALSE
  }else if(length(dependentVariableTypeVector)!=lengthSimulation){
    dependentVariableTypeVector=NA
    validDependentVariableTypeVectorInput=FALSE

    warning("dependentVariableTypeVector was not the same length as the output of the simulationFunction so replaced with NA")
  }



  plot_df=data.frame()


  for(i in seq(1,dim(simulationMatrix)[1] )){
    plot_df=rbind(plot_df, data.frame(simulation=simulationMatrix[i,], independentVariable=independentVariableVector, dependentVariableType=dependentVariableTypeVector))
  }

  median_vec=c()
  lower_percentile_vec=c()
  upper_percentile_vec=c()
  kind_df=unique(plot_df[,c("independentVariable", "dependentVariableType")])

  for(i in seq(1,dim(kind_df)[1])){
    kind=kind_df[i,]
    if(validDependentVariableTypeVectorInput){
      now_data_df=subset(plot_df, independentVariable==kind$independentVariable&dependentVariableType==kind$dependentVariableType )

    }else{
      now_data_df=subset(plot_df, independentVariable==kind$independentVariable )

    }

    nowQuantile=quantile(now_data_df$simulation, probs = c(confidenceLevels[1],0.5,confidenceLevels[2]), na.rm = TRUE)

    lower_percentile_vec=c(lower_percentile_vec,nowQuantile[1])
    median_vec=c(median_vec,nowQuantile[2])
    upper_percentile_vec=c(upper_percentile_vec,nowQuantile[3])
  }

  plot_CI_df=kind_df
  plot_CI_df$lower_percentile=as.numeric(lower_percentile_vec)
  plot_CI_df$upper_percentile=as.numeric(upper_percentile_vec)
  plot_CI_df$median=as.numeric(median_vec)


  if(validDependentVariableTypeVectorInput&overLay){
    g=ggplot2::ggplot(plot_CI_df, aes(x=independentVariable , y=median, colour=dependentVariableType ))+ggplot2::geom_line()+ggplot2::geom_ribbon(aes(ymin=lower_percentile, ymax=upper_percentile, fill=dependentVariableType), alpha=0.2)+ ggplot2::labs(caption = paste0("solid line is the median of the model prediction and shaded area is its confidence interval of ",confidenceLevels[1]*100,"-",confidenceLevels[2]*100," percentile"))

  }else{
    g=ggplot2::ggplot(plot_CI_df, aes(x=independentVariable , y=median ))+ggplot2::geom_line()+ggplot2::geom_ribbon(aes(ymin=lower_percentile, ymax=upper_percentile), alpha=0.2)+ ggplot2::labs(caption = paste0("solid line is the median of the model prediction and shaded area is its confidence interval of ",confidenceLevels[1]*100,"-",confidenceLevels[2]*100," percentile"))
    if(validDependentVariableTypeVectorInput){
      g=g+ggplot2::facet_wrap(.~dependentVariableType)
    }
  }

  if(ValidobservationVectorInput){
    if(overLay){
      g=g+ggplot2::geom_point(data=data.frame(independentVariable=observationIndpendentVariableVector, observation=observationVector, dependentVariableType=observationDependentVariableTypeVector), aes(x=independentVariable,y=observation, colour=dependentVariableType), inherit.aes = FALSE)

    }else{
      g=g+ggplot2::geom_point(data=data.frame(independentVariable=observationIndpendentVariableVector, observation=observationVector, dependentVariableType=observationDependentVariableTypeVector), colour="red", aes(x=independentVariable,y=observation), inherit.aes = FALSE)

    }
  }

  g=g+ggplot2::ylab("Dependent variable")+ggplot2::xlab("Independent variable")

  return(g)
}


