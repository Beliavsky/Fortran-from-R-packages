## TODO remove parallel statements from the documentation

repmat=function(subM,numrows,numcol){

  matrix_out=subM

  if(numrows>1){
    for(i in seq(1,numrows-1)){
      matrix_out=rbind(matrix_out,subM)
    }
  }

  tempM=matrix_out

  if(numcol>1){
    for(i in seq(1,numcol-1)){
      matrix_out=cbind(matrix_out,tempM)
    }
  }


  return(matrix_out)
}

col_normalize=function(data_in){
  out=data_in;
  for(i in seq(1,dim(data_in)[2])){
    if(sd(data_in[,i])!=0){
      out[,i]=(data_in[,i]-mean(data_in[,i]))/sd(data_in[,i])
    }
  }
  return(out)
}

col_sd=function(data_in){
  # out=c()
  # for(i in seq(1,dim(data_in)[2])){
  #   out=c(out,sd(data_in[,i],na.rm = TRUE))
  # }
  out=apply(data_in, 2, sd, na.rm=TRUE)

  return(out)
}

col_max=function(data_in){
  # out=c()
  # for(i in seq(1,dim(data_in)[2])){
  #   out=c(out,max(data_in[,i],na.rm = TRUE))
  # }
  out=apply(data_in, 2, max, na.rm=TRUE)

  return(out)
}

col_min=function(data_in){
  # out=c()
  # for(i in seq(1,dim(data_in)[2])){
  #   out=c(out,min(data_in[,i],na.rm = TRUE))
  # }
  out=apply(data_in, 2, min, na.rm=TRUE)

  return(out)
}

col_mean=function(data_in){
  # out=c()
  # for(i in seq(1,dim(data_in)[2])){
  #   out=c(out,mean(data_in[,i],na.rm = TRUE))
  # }

  out=apply(data_in, 2, mean, na.rm=TRUE)

  return(out)
}

col_median=function(data_in){
  # out=c()
  # for(i in seq(1,dim(data_in)[2])){
  #   out=c(out,median(data_in[,i],na.rm = TRUE))
  # }
  #
  out=apply(data_in, 2, median, na.rm=TRUE)
  return(out)
}


col_quantile=function(data_in, prob){
  # out=c()
  # for(i in seq(1,dim(data_in)[2])){
  #   out=c(out,as.numeric(quantile(data_in[,i],na.rm = TRUE, prob=prob)))
  # }
  out=apply(data_in, 2, quantile, prob=prob, na.rm=TRUE)

  return(out)
}

#' @title col_quantile
#' @description
#' Obtain columb wise quantile
#' @param data_in (required input) \emph{a matrix or a data.grame} where the column-wise quantile wishes to be determined.
#' @param prob (required input) \emph{a number} quantile expressed as in the probability.
#' @return \emph{a vector of number dim(data_in)[2]} containing: quantile of the each column where the probability is specified as "prob"
#' @examples
#' A=matrix(seq(1,100),nrow = 25)
#' col_quantile(A, 0.5)
#' @export
col_quantile=function(data_in, prob){
  out=c()
  for(i in seq(1,dim(data_in)[2])){
    out=c(out,as.numeric(quantile(data_in[,i],na.rm = TRUE, prob=prob)))
  }
  return(out)
}


elbow_index=function(vec_in){

  trapizoido_area=c()
  for(i in seq(1,length(vec_in)-1))
    trapizoido_area=c(trapizoido_area, (vec_in[1]+vec_in[i])*i+(vec_in[i+1]+vec_in[length(vec_in)])*(length(vec_in)-i))

  return(which(trapizoido_area==min(trapizoido_area)))
}

optimal_kmeans=function(data_in){
  #Assume data is stored in rows

  normalized_data=col_normalize(data_in)

  residual=c()

  numMaxCluster=min(round(dim(data_in)[1]/dim(data_in)[2]/2),10)

  if(numMaxCluster>1){

    for(i in seq(2,numMaxCluster)){
      kmeans_res=kmeans(normalized_data,i)
      residual=c(residual,kmeans_res$tot.withinss)
    }
    elbow_ind=(elbow_index(residual))

    if(length(elbow_ind)>0){
      kmeans(normalized_data,(elbow_ind[1]+1))
    }else{
      kmeans(normalized_data,1)
    }
  }else{
    kmeans(normalized_data,1)
  }

}

matlabSum <-function(matrix,direction){
  if (direction==1){
    return(t(as.matrix(colSums(matrix))))
  }else{
    return(as.matrix(rowSums(matrix)))
  }
}

matlabMedian<-function(matrix,direction){
  out=c()
  if (direction==1){
    for(i in seq(1,dim(matrix)[2])){
      out=c(out,median(matrix[,i], na.rm =TRUE))
    }
    return(t(as.matrix(out)))
  }else{
    for(i in seq(1,dim(matrix)[1])){
      out=c(out,median(matrix[i,], na.rm =TRUE))
    }
    return((as.matrix(out)))
  }
}

dot=function(vec1,vec2){
  return(as.numeric(t(as.matrix(vec1))%*%as.matrix(vec2)))
}

matlabRand<-function(numRow,numCol){
  outMat=runif(numCol)
  for(i in seq(2,numRow)){
    outMat=rbind(outMat,runif(numCol))
  }
  return(outMat)
}

tryCatch_nonlinearFunction=function(x,num_observations, nonlinearFunction, validTarget, ... ){

  if(is.null(dim(x))){
    numXrows=1
  }else{
    numXrows=dim(x)[1]
  }

  out <- tryCatch(
    {
      nonlinearFunction(x, ... )

    },
    error=function(cond) {
      if(numXrows==1){
        temoOut=rep(NaN,num_observations)
      }else{
        temoOut=matrix(rep(NaN,num_observations*numXrows),nrow = numXrows)
      }
      return(temoOut)
    },
    warning=function(cond) {
      nonlinearFunction(x, ... )
    },
    finally={

    }
  )
  if(numXrows==1){
    out[!validTarget]=0

    if(length(out)!=num_observations){
      out=rep(NaN,num_observations)
    }else if(is.na(sum(out))){
      out=rep(NaN,num_observations)
    }else if(is.infinite(sum(abs(out)))){
      out=rep(NaN,num_observations)
    }
  }else{

    if(dim(out)[2]!=num_observations){
      out=matrix(rep(NaN,num_observations*numXrows),nrow = numXrows)
    }else{
      out[is.na(rowSums(out)),]=NaN
      out[is.infinite(abs(rowSums(out))),]=NaN
      out[is.nan(rowSums(out)),]=NaN
    }

    out[,!validTarget]=0

  }


  return(out)
}


CGNM_input_test <- function(nonlinearFunction,  targetVector, initial_lowerRange, initial_upperRange, ... ){

  tempOut <- tryCatch(
    {
      (dim(nonlinearFunction(t(matrix(rep(initial_upperRange,2),ncol = 2)), ... )))

    },
    error=function(cond) {
      return(NULL)
    },
    warning=function(cond) {
      (dim(nonlinearFunction(t(matrix(rep(initial_upperRange,2),ncol = 2)), ... )))
    },
    finally={

    }
  )


  if(is.null(tempOut)){
    out=FALSE
  }else{
    out=(tempOut[1]==2)
  }


  if(length(initial_lowerRange)!=length(initial_upperRange)){
    stop("initial_lowerRange and initial_upperRange need to have the same length")

  }

  if(sum(is.na(as.numeric(initial_lowerRange)))+sum(is.na(as.numeric(initial_upperRange)))>0){
    stop("initial_lowerRange and initial_upperRange can only be the vector of non-numericl values")

  }

  if(out){

    message("nonlinearFunction is given as matrix to matrix function")
    testEval=nonlinearFunction(t(matrix(c(initial_lowerRange,(initial_upperRange+initial_lowerRange)/2,initial_upperRange),ncol = 3)), ... )
    RS_testEval=rowSums(testEval)

    if(is.na(RS_testEval[1])){
      message("WARNING: nonlinearFunction evaluation at initial_lowerRange NOT Successful.")
    }else{
      message("NonlinearFunction evaluation at initial_lowerRange Successful.")
    }

    if(is.na(RS_testEval[2])){
      message("WARNING: nonlinearFunction evaluation at (initial_upperRange+initial_lowerRange)/2 NOT Successful.")
    }else{
      message("NonlinearFunction evaluation at (initial_upperRange+initial_lowerRange)/2 Successful.")
    }

    if(is.na(RS_testEval[3])){
      message("WARNING: nonlinearFunction evaluation at initial_upperRange NOT Successful.")
    }else{
      message("NonlinearFunction evaluation at initial_upperRange Successful.")
    }


    if(dim(testEval)[2]!=length(targetVector)){
      stop("Length of output of the nonlinearFunction and targetVector are not the same. Double check the nonlinear function definition and targetVector input.")
    }

  }else{


    message("checking if the nonlinearFunction can be evaluated at the initial_lowerRange")
    testEval=nonlinearFunction(initial_lowerRange, ... )
    if(!is.na(sum(testEval))&&length(testEval)==length(targetVector)){
      message("Evaluation Successful")
    }else{
      message("WARNING: nonlinearFunction evaluation at initial_lowerRange NOT Successful.")
    }

    message("checking if the nonlinearFunction can be evaluated at the initial_upperRange")
    testEval=nonlinearFunction(initial_upperRange, ... )
    if(!is.na(sum(testEval))&&length(testEval)==length(targetVector)){
      message("Evaluation Successful")
    }else{
      message("WARNING: nonlinearFunction evaluation at initial_upperRange NOT Successful.")
    }

    message("checking if the nonlinearFunction can be evaluated at the (initial_upperRange+initial_lowerRange)/2")
    testEval=nonlinearFunction((initial_upperRange+initial_lowerRange)/2, ... )
    if(!is.na(sum(testEval))&&length(testEval)==length(targetVector)){
      message("Evaluation Successful")
    }else{
      message("WARNING: nonlinearFunction evaluation at (initial_upperRange+initial_lowerRange)/2 NOT Successful.")
    }

    if(length(testEval)!=length(targetVector)){
      stop("Length of output of the nonlinearFunction and targetVector are not the same. Double check the nonlinear function definition and targetVector input.")
    }
  }
  return(out)
}



#' @title Cluster_Gauss_Newton_method
#' @description Find multiple minimisers of the nonlinear least squares problem.
#' \deqn{argmin_x ||f(x)-y*||}
#' where
#' \enumerate{\item f: nonlinear function (e.g., mathematical model)
#' \item y*: target vector (e.g., observed data to fit the mathematical model)
#' \item x: variable of the nonlinear function that we aim to find the values that minimize (minimizers) the differences between the nonlinear function and target vector (e.g., model parameter)
#' }
#' Parameter estimation problems of mathematical models can often be formulated as nonlinear least squares problems.  In this context f can be thought at a model, x is the parameter, and y* is the observation.
#' CGNM iteratively estimates the minimizer of the nonlinear least squares problem from various initial estimates hence finds multiple minimizers.
#' Full detail of the algorithm and comparison with conventional method is available in the following publication, also please cite this publication when this algorithm is used in your research: Aoki et al. (2020) <doi.org/10.1007/s11081-020-09571-2>. Cluster Gauss–Newton method. Optimization and Engineering, 1-31.  As illustrated in this paper, CGNM is faster and more robust compared to repeatedly applying the conventional optimization/nonlinear least squares algorithm from various initial estimates. In addition, CGNM can realize this speed assuming the nonlinear function to be a black-box function (e.g. does not use things like adjoint equation of a system of ODE as the function does not have to be based on a system of ODEs.).
#' @param nonlinearFunction (required input) \emph{A function with input of a vector x of real number of length n and output a vector y of real number of length m.} In the context of model fitting the nonlinearFunction is \strong{the model}.  Given the CGNM does not assume the uniqueness of the minimizer, m can be less than n.  Also CGNM does not assume any particular form of the nonlinear function and also does not require the function to be continuously differentiable (see Appendix D of our publication for an example when this function is discontinuous). Also this function can be matrix to matrix equation.  This can be used for parallerization, see vignettes for examples.
#' @param targetVector (required input) \emph{A vector of real number of length m} where we minimize the Euclidean distance between the nonlinearFuncition and targetVector.  In the context of curve fitting targetVector can be though as \strong{the observational data}.
#' @param initial_lowerRange (required input) \emph{A vector of real number of length n} where each element represents  \strong{the lower range of the initial iterate}. Similarly to regular Gauss-Newton method, CGNM iteratively reduce the residual to find minimizers.  Essential differences is that CGNM start from the initial RANGE and not an initial point.
#' @param initial_upperRange (required input) \emph{A vector of real number of length n} where each element represents  \strong{the upper range of the initial iterate}.
#' @param lowerBound (default: NA) \emph{A vector of real number or NA of length n} where each element represents  \strong{the lower bound of the parameter search}.  If no lower bound set that element NA. Note that CGNM is an unconstraint optimization method so the final minimizer can be anywhere.  In the parameter estimation problem, there often is a constraints to the parameters (e.g., parameters cannot be negative). So when the upper or lower bound is set using this option, parameter transformation is conducted internally (e.g., if either the upper or lower bound is given parameters are log transformed, if the upper and lower bounds are given logit transform is used.)
#' @param upperBound (default: NA) \emph{A vector of real number or NA of length n} where each element represents  \strong{the upper bound of the parameter search}.  If no upper bound set that element NA.
#' @param ParameterNames (default: NA) \emph{A vector of string} of length n User can specify names of the parameters that will be used for the plots.
#' @param stayIn_initialRange (default: FALSE) \emph{TRUE or FALSE} if set TRUE, the parameter search will conducted strictly within the range specified by initial_lowerRange and initial_upperRange.
#' @param num_minimizersToFind  (default: 250) \emph{A positive integer} defining number of approximate minimizers CGNM will find. We usually \strong{use 250 when testing the model and 1000 for the final analysis}.  The computational cost increase proportionally to this number; however, larger number algorithm becomes more stable and increase the chance of finding more better minimizers. See Appendix C of our paper for detail.
#' @param num_iteration (default: 25)  \emph{A positive integer} defining maximum number of iterations. We usually \strong{set 25 while model building and 100 for final analysis}.  Given each point terminates the computation when the convergence criterion is met the computation cost does not grow proportionally to the number of iterations (hence safe to increase this without significant increase in the computational cost).
#' @param saveLog (default: TRUE) \emph{TRUE or FALSE} indicating either or not to save computation result from each iteration in CGNM_log folder. It requires disk write access right in the current working directory. \strong{Recommended to set TRUE if the computation is expected to take long time} as user can retrieve intrim computation result even if the computation is terminated prematurely (or even during the computation).
#' @param runName (default: "") \emph{string} that user can ue to identify the CGNM runs. The run history will be saved in the folder name CGNM_log_<runName>.  If this is set to "TIME" then runName is automatically set by the run start time.
#' @param textMemo (default: "") \emph{string} that user can write an arbitrary text (without influencing computation). This text is stored with the computation result so that can be used for example to describe model so that the user can recognize the computation result.
#' @param algorithmParameter_initialLambda (default: 1) \emph{A positive number} for initial value for the regularization coefficient lambda see Appendix B of of our paper for detail.
#' @param algorithmParameter_gamma (default: 2) \emph{A positive number} a positive scalar value for adjusting the strength of the weighting for the linear approximation see Appendix A of our paper for detail.
#' @param algorithmVersion (default: 3.0) \emph{A positive number} user can choose different version of CGNM algorithm currently 1.0 and 3.0 are available.  If number chosen other than 1.0 or 3.0 it will choose 1.0.
#' @param initialIterateMatrix (default: NA) \emph{A matrix} with dimension num_minimizersToFind x n.  User can provide initial iterate as a matrix  This input is used when the user wishes not to generate initial iterate randomly from the initial range.  The user is responsible for ensuring all function evaluation at each initial iterate does not produce NaN.
#' @param targetMatrix (default: NA) \emph{A matrix} with dimension num_minimizersToFind x m  User can define multiple target vectors in the matrix form.  This input is mainly used when running bootstrap method and not intended to be used for other purposes.
#' @param keepInitialDistribution (default: NA) \emph{A vector of TRUE or FALSE} of length n User can specify if the initial distribution of one of the input variable (e.g. parameter) to be kept as the initial iterate throughout CGNM iterations.
#' @param weightMatrix (default: NA) \emph{A matrix} with dimension num_minimizersToFind x m  User can define multiple weight vectors in the matrix form to weight the observations.  This input is mainly used when running case sampling bootstrap method and not intended to be used for other purposes.
#' @param MO_weights (default: NA) \emph{A numeric vector} where the weights for the middle out methods are specified.  The length of the vector should be the same as the number of parameters. MO can be used to incoperate prior knowledge of the parameter to be estimated, weight indicate an arbitrary confidence for the prior information. (MO method is still under methodological development.)
#' @param MO_values (default: NA) \emph{A numeric vector} where the values for the middle out methods are specified.  The length of the vector should be the same as the number of parameters. MO can be used to incoperate prior knowledge of the parameter to be estimated. (MO method is still under methodological development.)
#' @param ... Further arguments to be supplied to nonlinearFunction
#' @return list of a matrix X, Y,residual_history and initialX, as well as a list runSetting
#' \enumerate{\item X: \emph{a num_minimizersToFind by n matrix} which stores the approximate minimizers of the nonlinear least squares in each row. In the context of model fitting they are \strong{the estimated parameter sets}.
#' \item Y: \emph{a num_minimizersToFind by m matrix} which stores the nonlinearFunction evaluated at the corresponding approximate minimizers in matrix X above. In the context of model fitting each row corresponds to \strong{the model simulations}.
#' \item residual_history: \emph{a num_iteration by num_minimizersToFind matrix} storing sum of squares residual for all iterations.
#' \item initialX: \emph{a num_minimizersToFind by n matrix} which stores the set of initial iterates.
#' \item runSetting: a list containing all the input variables to Cluster_Gauss_Newton_method (i.e., nonlinearFunction, targetVector, initial_lowerRange, initial_upperRange ,algorithmParameter_initialLambda, algorithmParameter_gamma, num_minimizersToFind, num_iteration, saveLog, runName, textMemo).}
#' @examples
#' ##lip-flop kinetics (an example known to have two distinct solutions)
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
#' targetVector = observation, num_iteration = 10, num_minimizersToFind = 100,
#' initial_lowerRange = c(0.1,0.1,0.1), initial_upperRange =  c(10,10,10),
#' saveLog = FALSE)
#'
#' acceptedApproximateMinimizers(CGNM_result)
#'
#  ## flip-flop kinetics using RxODE (an example known to have two distinct solutions)
#' \dontrun{
#' library(RxODE)
#'
#' model_text="
#' d/dt(X_1)=-ka*X_1
#' d/dt(C_2)=(ka*X_1-CL_2*C_2)/V1"
#'
#' model=RxODE(model_text)
#' #define nonlinearFunction
#' model_function=function(x){
#'
#' observation_time=c(0.1,0.2,0.4,0.6,1,2,3,6,12)
#'
#' theta <- c(ka=x[1],V1=x[2],CL_2=x[3])
#' ev <- eventTable()
#' ev$add.dosing(dose = 1000, start.time =0)
#' ev$add.sampling(observation_time)
#' odeSol=model$solve(theta, ev)
#' log10(odeSol[,"C_2"])
#'
#' }
#'
#' observation=log10(c(4.91, 8.65, 12.4, 18.7, 24.3, 24.5, 18.4, 4.66, 0.238))
#'
#' CGNM_result=Cluster_Gauss_Newton_method(nonlinearFunction=model_function,
#' targetVector = observation, saveLog = FALSE,
#' initial_lowerRange = c(0.1,0.1,0.1),initial_upperRange =  c(10,10,10))}
#'
#' @export
#' @import stats MASS
Cluster_Gauss_Newton_method <- function(nonlinearFunction,targetVector, initial_lowerRange, initial_upperRange , lowerBound=NA, upperBound=NA, ParameterNames=NA, stayIn_initialRange=FALSE, num_minimizersToFind=250, num_iteration=25, saveLog=TRUE, runName="", textMemo="",algorithmParameter_initialLambda=1, algorithmParameter_gamma=2, algorithmVersion=3.0, initialIterateMatrix=NA, targetMatrix=NA, weightMatrix=NA, keepInitialDistribution=NA, MO_weights=NA, MO_values=NA,... ){
  MO_weights_in=MO_weights
  MO_values_in=MO_values

  saveFolderName="CGNM_log"
  if(runName=="TIME"){
    runName=as.character( Sys.time())
    runSetting$runName=runName
  }

  if(runName!=""){
    saveFolderName=paste0(runName,"_",saveFolderName)
  }
  #  saveFolderName=paste0(saveFolderName,ifelse(is.na(targetMatrix)[1],"",ifelse(all(targetMatrix[1,]==targetMatrix[2,]),"","_trap")))
  #  saveFolderName=paste0(saveFolderName,ifelse(is.na(weightMatrix)[1],"",ifelse(all(weightMatrix[1,]==weightMatrix[2,]),"","_trap")))

  #  saveFolderName=paste0(saveFolderName,ifelse(is.na(targetMatrix)[1],"",ifelse(all(targetMatrix[1,]==targetMatrix[2,])&all(weightMatrix[1,]==weightMatrix[2,]),"","_trap")))

  print(saveFolderName)
  if(saveLog){
    dir.create(saveFolderName)
  }


  if(length(lowerBound)==1&is.na(lowerBound)[1]){
    lowerBound=rep(NA, length(initial_lowerRange))
  }

  if(length(upperBound)==1&is.na(upperBound)[1]){
    upperBound=rep(NA, length(initial_upperRange))
  }

  if(length(MO_weights)==1&is.na(MO_weights)[1]){
    MO_weights=rep(0, length(initial_lowerRange))
    MO_values=(initial_lowerRange+initial_upperRange)/2
  }

  if(!( length(initial_lowerRange)==length(initial_upperRange))){
    stop("length of initial_lowerRange, initial_upperRange must be the same")
  }

  if(!(length(lowerBound)==length(upperBound)&length(initial_upperRange)==length(lowerBound))){
    stop("length of lowerBound, upperBound must be the same as initial_lowerRange")
  }

  if(!(length(initial_lowerRange)==length(MO_weights) & length(initial_lowerRange)==length(MO_values) )){
    stop("length of MO_weights, MO_values must be the same as initial_lowerRange")
  }

  length_paraVector=length(initial_lowerRange)


  if(!is.numeric(MO_weights)|sum(is.na(MO_weights))>0){
    stop("MO_weights needs to be numeric.")
  }

  MO_index=(MO_weights!=0)

  if(!is.numeric(MO_values[MO_index])|sum(is.na(MO_values[MO_index]))>0){
    stop("MO_values where the weights are non-zero must be numeric.")

  }

  validIndex=MO_index
  if((sum(MO_values[validIndex]>initial_lowerRange[validIndex])!=sum(validIndex))){
    stop("MO_values must be STRICTRY larger than initial lower range")
  }

  if((sum(MO_values[validIndex]<initial_upperRange[validIndex])!=sum(validIndex))){
    stop("MO_values must be STRICTRY smaller than initial upper range")
  }

  MO_values[!MO_index]=(initial_upperRange[!MO_index]+initial_lowerRange[!MO_index])/2

  validIndex=!is.na(lowerBound)
  if((sum(initial_lowerRange[validIndex]>lowerBound[validIndex])!=sum(validIndex))){
    stop("Initial lower range must be STRICTRY larger than lower BOUND")
  }

  validIndex=!is.na(upperBound)
  if((sum(initial_upperRange[validIndex]<upperBound[validIndex])!=sum(validIndex))){
    stop("Initial upper range must be STRICTRY less than upper BOUND")
  }


  if(is.na(ParameterNames)[1]){
    ParameterNames=paste0("theta_",seq(1,length_paraVector))
  }

  ReparameterizationDef=paste0("x",seq(1,length_paraVector))

  LB_index=(!is.na(as.numeric(lowerBound)))&is.na(upperBound)
  UB_index=(!is.na(as.numeric(upperBound)))&is.na(lowerBound)
  BB_index=(!is.na(as.numeric(lowerBound)))&(!is.na(as.numeric(upperBound)))

  if(sum(LB_index)+sum(UB_index)+sum(BB_index)>0){
    message(paste(paste(c(ParameterNames[LB_index],ParameterNames[UB_index],ParameterNames[BB_index]), collapse = ", "),"will be transformed internally to impose boundaries. See CGNM_result$runSetting$ReparameterizationDef for exact transformation. Transformed parameters are denoted as x and untransformed parameters are denoted as theta."))
  }

  LB_paraReDef=paste0("exp(x",seq(1,length(lowerBound)),")+",lowerBound)
  UB_paraReDef=paste0(upperBound,"-exp(x",seq(1,length(lowerBound)),")")
  BB_paraReDef=paste0("(",upperBound,"-",lowerBound,")*(exp(x",seq(1,length(lowerBound)),")/(exp(x",seq(1,length(lowerBound)),")+1))+",lowerBound)

  ReparameterizationDef[LB_index]=LB_paraReDef[LB_index]
  ReparameterizationDef[UB_index]=UB_paraReDef[UB_index]
  ReparameterizationDef[BB_index]=BB_paraReDef[BB_index]

  X_LR=initial_lowerRange
  X_LR[LB_index]=log(initial_lowerRange[LB_index]-lowerBound[LB_index])
  X_LR[UB_index]=log(-initial_lowerRange[UB_index]+upperBound[UB_index])
  X_LR[BB_index]=log((initial_lowerRange[BB_index]-lowerBound[BB_index])/(upperBound[BB_index]-initial_lowerRange[BB_index]))

  ParameterTransFormationDef=rep("", length(ParameterNames))
  ParameterTransFormationDef[LB_index]=paste0("log(",ParameterNames[LB_index],"-",lowerBound[LB_index],")")
  ParameterTransFormationDef[UB_index]=paste0("log(-",ParameterNames[UB_index],"+",upperBound[UB_index],")")
  ParameterTransFormationDef[BB_index]=paste0("log((",ParameterNames[BB_index],"-",lowerBound[BB_index],")/(",upperBound[BB_index],"-",ParameterNames[BB_index],"))")


  X_UR=initial_upperRange
  X_UR[LB_index]=log(initial_upperRange[LB_index]-lowerBound[LB_index])
  X_UR[UB_index]=log(-initial_upperRange[UB_index]+upperBound[UB_index])
  X_UR[BB_index]=log((initial_upperRange[BB_index]-lowerBound[BB_index])/(upperBound[BB_index]-initial_upperRange[BB_index]))


  MO_val_converted=MO_values
  MO_val_converted[LB_index]=log(MO_val_converted[LB_index]-lowerBound[LB_index])
  MO_val_converted[UB_index]=log(-MO_val_converted[UB_index]+upperBound[UB_index])
  MO_val_converted[BB_index]=log((MO_val_converted[BB_index]-lowerBound[BB_index])/(upperBound[BB_index]-MO_val_converted[BB_index]))


  if(!is.na(initialIterateMatrix)[1]){
    originalInit=initialIterateMatrix
    rowNumInitMatrix=dim(initialIterateMatrix)[1]
    initialIterateMatrix[,LB_index]=log(initialIterateMatrix[,LB_index]-repmat(lowerBound[LB_index],numrows = rowNumInitMatrix ,numcol = 1))
    initialIterateMatrix[,UB_index]=log(-initialIterateMatrix[,UB_index]+repmat(upperBound[UB_index],numrows = rowNumInitMatrix ,numcol = 1))
    initialIterateMatrix[,BB_index]=log((initialIterateMatrix[,BB_index]-repmat(lowerBound[BB_index],numrows = rowNumInitMatrix ,numcol = 1))/(repmat(upperBound[BB_index],numrows = rowNumInitMatrix ,numcol = 1)-initialIterateMatrix[,BB_index]))

  }

  lowerBoundMatrix=repmat(lowerBound[BB_index],numrows = num_minimizersToFind ,numcol = 1)
  upperBoundMatrix=repmat(upperBound[BB_index],numrows = num_minimizersToFind ,numcol = 1)

  if(sum(abs(MO_weights))>0){
    nonlinearFunction_varTrans=function(x, ... ){

      if(is.vector(x)){
        theta=x
        theta[LB_index]=exp(x[LB_index])-lowerBound[LB_index]
        theta[UB_index]=-exp(x[UB_index])+upperBound[UB_index]
        theta[BB_index]=(upperBound[BB_index]-lowerBound[BB_index])*(exp(x[BB_index])/(exp(x[BB_index])+1))+lowerBound[BB_index]
        out=nonlinearFunction(theta, ... )

      }else{
        theta=x

        theta[,LB_index]=exp(x[,LB_index])-repmat(lowerBound[LB_index],numrows = dim(x)[1] ,numcol = 1)
        theta[,UB_index]=-exp(x[,UB_index])+repmat(upperBound[UB_index],numrows = dim(x)[1] ,numcol = 1)
        theta[,BB_index]=(upperBoundMatrix-lowerBoundMatrix)[seq(1,dim(x)[1]),]*(exp(x[,BB_index])/(exp(x[,BB_index])+1))+lowerBoundMatrix[seq(1,dim(x)[1]),]

        out=nonlinearFunction(theta, ... )
      }

      if(is.vector(out)){
        out=c(out,x[MO_index]*MO_weights[MO_index])
      }else{
        for(MO_ind_nu in seq(1,length(MO_index))[MO_index]){
          out=cbind(out,matrix(x[,MO_ind_nu]*MO_weights[MO_ind_nu],ncol = 1))
        }
      }

      return(out)

    }


    targetVector=c(targetVector, MO_val_converted[MO_index]*MO_weights[MO_index])

  }else{
    nonlinearFunction_varTrans=function(x, ... ){

      if(is.vector(x)){
        theta=x
        theta[LB_index]=exp(x[LB_index])-lowerBound[LB_index]
        theta[UB_index]=-exp(x[UB_index])+upperBound[UB_index]
        theta[BB_index]=(upperBound[BB_index]-lowerBound[BB_index])*(exp(x[BB_index])/(exp(x[BB_index])+1))+lowerBound[BB_index]

      }else{
        theta=x

        theta[,LB_index]=exp(x[,LB_index])-repmat(lowerBound[LB_index],numrows = dim(x)[1] ,numcol = 1)
        theta[,UB_index]=-exp(x[,UB_index])+repmat(upperBound[UB_index],numrows = dim(x)[1] ,numcol = 1)
        theta[,BB_index]=(upperBoundMatrix-lowerBoundMatrix)[seq(1,dim(x)[1]),]*(exp(x[,BB_index])/(exp(x[,BB_index])+1))+lowerBoundMatrix[seq(1,dim(x)[1]),]
      }

      nonlinearFunction(theta, ... )
    }
  }


  textMemo=paste0("CGNM algorithm version ",algorithmVersion,"\n\n",textMemo)
  runSetting_list=list(nonlinearFunction=nonlinearFunction,
                       targetVector=targetVector,
                       initial_lowerRange=initial_lowerRange,
                       initial_upperRange=initial_upperRange,
                       stayIn_initialRange=stayIn_initialRange,
                       algorithmParameter_initialLambda=algorithmParameter_initialLambda,
                       algorithmParameter_gamma=algorithmParameter_gamma,
                       num_minimizersToFind=num_minimizersToFind,
                       num_iteration=num_iteration,
                       saveLog=saveLog,
                       runName=runName,
                       textMemo=textMemo,
                       algorithmVersion=algorithmVersion,
                       initialIterateMatrix=initialIterateMatrix,
                       targetMatrix=targetMatrix,
                       ParameterNames=ParameterNames,
                       ReparameterizationDef=ReparameterizationDef,
                       ParameterTransFormationDef=ParameterTransFormationDef,
                       lowerBound=lowerBound,
                       upperBound=upperBound,
                       weightMatrix=weightMatrix,
                       keepInitialDistribution=keepInitialDistribution,
                       targetMatrix=targetMatrix,
                       MO_weights=MO_weights_in,
                       MO_values=MO_values_in)


  if(saveLog){
    save(file=paste0(saveFolderName,'/runSetting_list.RDATA'), runSetting_list)
  }


  out=Cluster_Gauss_Newton_method_core(nonlinearFunction_varTrans, targetVector=targetVector, initial_lowerRange=X_LR, initial_upperRange=X_UR ,lowerBound=lowerBound, upperBound=upperBound, stayIn_initialRange, num_minimizersToFind, num_iteration, saveLog, runName, textMemo,algorithmParameter_initialLambda, algorithmParameter_gamma, algorithmVersion, initialIterateMatrix, targetMatrix, keepInitialDistribution, weightMatrix,ParameterNames,ReparameterizationDef, ParameterTransFormationDef, ... )

  out_X=as.data.frame(out$X)
  names(out_X)=paste0("x",seq(1,dim(out_X)[2]))

  Theta=as.data.frame(matrix(NA,nrow = dim(out_X)[1], ncol = length(ParameterNames)))
  names(Theta)=ParameterNames

  for(i in seq(1,length(ParameterNames))){
    Theta[,ParameterNames[i]]=with(out_X, eval(parse(text=ReparameterizationDef[i])))
  }

  out$finalParameterCombinations=Theta
  out$runSetting$MO_values=MO_values
  out$runSetting$MO_weights=MO_weights
  out$runSetting$ParameterTransFormationDef=ParameterTransFormationDef

  CGNM_result=out

  if(saveLog){
    save(file=paste0(saveFolderName,'/iteration_final.RDATA'), CGNM_result)
  }


  return(out)
}


Cluster_Gauss_Newton_method_core <- function(nonlinearFunction, targetVector, initial_lowerRange, initial_upperRange , stayIn_initialRange=FALSE, num_minimizersToFind=250, num_iteration=25, saveLog=FALSE, runName="", textMemo="",algorithmParameter_initialLambda=1, algorithmParameter_gamma=2, algorithmVersion=3.0, initialIterateMatrix=NA, targetMatrix=NA, keepInitialDistribution=NA,weightMatrix=NA,ParameterNames=NA,ReparameterizationDef=NA, ParameterTransFormationDef=NA, lowerBound=NA, upperBound=NA, folderSuffix="", ... ){

  if(sum(is.na(targetVector))>0|sum(is.nan(targetVector))>0){
    stop("targetVector contains NA or NAN")
  }

  if(sum(initial_lowerRange==initial_upperRange)==length(initial_upperRange)){
    stop("all initial_upperRange and initial_lowerRange values are the same, hence no parameter to estimate.")
  }

  CGNM_start_time=Sys.time()
  targetMatrix_in=targetMatrix

  dimTargetMatrix=c(0,0)
  dimInitialIterateMatrix=c(0,0)

  if(is.matrix(initialIterateMatrix)){
    dimInitialIterateMatrix=dim(initialIterateMatrix)
  }

  if(is.matrix(targetMatrix_in)){
    dimTargetMatrix=dim(targetMatrix_in)
  }

  if(dimInitialIterateMatrix[1]>0&dimTargetMatrix[1]>0){
    num_minimizersToFind=min(dimInitialIterateMatrix[1], dimTargetMatrix[1])
    initialIterateMatrix=initialIterateMatrix[1:num_minimizersToFind,]
    targetMatrix_in=targetMatrix_in[1:num_minimizersToFind,]
  }else if(dimTargetMatrix[1]>0){
    num_minimizersToFind=dimTargetMatrix[1]
  }else if(dimInitialIterateMatrix[1]>0){
    num_minimizersToFind=dimInitialIterateMatrix[1]
  }

  weightMatrix_forRes=weightMatrix

  if(!is.matrix(weightMatrix_forRes)){
    weightMatrix_forRes=matrix(1, nrow = num_minimizersToFind, ncol = length(targetVector))
  }else if(dim(weightMatrix_forRes)[1]!=num_minimizersToFind|dim(weightMatrix_forRes)[2]!=length(targetVector)){
    message("Wrong size of weightMatrix is given so all weight will be set to one.")
    weightMatrix_forRes=matrix(1, nrow = num_minimizersToFind, ncol = length(targetVector))
  }


  textMemo=paste0("CGNM algorithm version ",algorithmVersion,"\n\n",textMemo)
  runSetting=list(nonlinearFunction=nonlinearFunction, targetVector=targetVector, initial_lowerRange=initial_lowerRange, initial_upperRange=initial_upperRange, stayIn_initialRange=stayIn_initialRange,algorithmParameter_initialLambda=algorithmParameter_initialLambda, algorithmParameter_gamma=algorithmParameter_gamma, num_minimizersToFind=num_minimizersToFind, num_iteration=num_iteration, saveLog=saveLog, runName=runName, textMemo=textMemo, algorithmVersion=algorithmVersion, initialIterateMatrix=initialIterateMatrix, targetMatrix=targetMatrix,  keepInitialDistribution=keepInitialDistribution,weightMatrix=weightMatrix, ParameterNames=ParameterNames,ReparameterizationDef=ReparameterizationDef, ParameterTransFormationDef=ParameterTransFormationDef,lowerBound=lowerBound, upperBound=upperBound)

  testCGNMinput=TRUE
  nonlinearFunctionCanTakeMatrixInput=FALSE

  if(testCGNMinput){

    testCGNMinput_time_start=Sys.time()
    nonlinearFunctionCanTakeMatrixInput=CGNM_input_test(nonlinearFunction , targetVector, initial_lowerRange, initial_upperRange , ... )
    testCGNMinput_time_end=Sys.time()

    #    message(paste("Initial estimation of required computation time:",round((testCGNMinput_time_end-testCGNMinput_time_start)/3*num_minimizersToFind*num_iteration/60),"min"))
   #message(paste("CGNM iteration should finish before:",Sys.time()+(testCGNMinput_time_end-testCGNMinput_time_start)/3*num_minimizersToFind*num_iteration))


  }

  showIntermetiateResults=FALSE

  targetVector=as.numeric(targetVector)
  validTarget_vec=!is.na(targetVector)

  targetVector[!validTarget_vec]=0

  if(dimTargetMatrix[1]==0){
    targetMatrix=repmat(matrix(targetVector,nrow=1), num_minimizersToFind,1)
  }else{
    targetMatrix=targetMatrix_in
    targetMatrix[,!validTarget_vec]=0
  }


  X_ul_in=t(matrix(c(initial_lowerRange,initial_upperRange), nrow=length(initial_upperRange)))

  saveFolderName="CGNM_log"
  if(runName=="TIME"){
    runName=as.character( Sys.time())
    runSetting$runName=runName
  }

  if(runName!=""){
    saveFolderName=paste0(runName,"_",saveFolderName)
  }
  #  saveFolderName=paste0(saveFolderName,ifelse(is.na(targetMatrix)[1],"",ifelse(all(targetMatrix[1,]==targetMatrix[2,])&all(weightMatrix[1,]==weightMatrix[2,]),"","_bootstrap")))

  #saveFolderName=paste0(saveFolderName,ifelse(is.na(targetMatrix)[1],"",ifelse(all(targetMatrix[1,]==targetMatrix[2,]),"","_bootstrap")))
  #saveFolderName=paste0(saveFolderName,ifelse(is.na(weightMatrix)[1],"",ifelse(all(weightMatrix[1,]==weightMatrix[2,]),"","_bootstrap")))
  saveFolderName=paste0(saveFolderName,ifelse(folderSuffix=="","",paste0("_",folderSuffix)))

  print(saveFolderName)

  if(saveLog){
    dir.create(saveFolderName)
  }

  X_history=c()
  Y_history=c()

  algorithmParameter_gamma=algorithmParameter_gamma

  method <- 'CGNM'# set algorithm name for the log files)
  descriptionText <- paste0(toString(num_minimizersToFind),'samples_initLambda',toString(algorithmParameter_initialLambda),'_distanceOrder',toString(algorithmParameter_gamma))

  # setting the name depending on if we do restart or not.  Restart is an
  # experimental feature so I suggest not to use at this point.

  num_parameters <- dim(X_ul_in)[2]
  num_observations = length(targetVector)

  timeOneParaTook <- c()

  #  tic

  ############################
  # 1) Pre-iteration process #
  ############################
  X <- matrix(0,num_minimizersToFind,num_parameters)# initialize X (We keep the same naming for the variables as in the manuscript)
  Y <- matrix(1,num_minimizersToFind,num_observations)# initialize Y
  Y_new <- matrix(1,num_minimizersToFind,num_observations)# initialize Y
  lambda_vec <- matrix(1,1,num_minimizersToFind)*algorithmParameter_initialLambda# initialise regularisation parameter lambda

  residual_history <- c()# initialise a matrix that stores the SSR for all iterations


  ## Generate initial cluster
  is_alive <- matrix(0,1,num_minimizersToFind)# initialise the vector to keep track of if the nonlinear function was able to be evaluated at the randomly generated x.


  userDefInitialIterateAvailable=(dimInitialIterateMatrix[1]>0)

  # repeat the randomsampling of x until num_minimizersToFind of 'valid' x are sampled. Where we consider x to be valid if f(x) can be evaluated.
  while(sum(is_alive)<num_minimizersToFind){

    # random sampling of X
    if(userDefInitialIterateAvailable){
      X_temp=initialIterateMatrix
      userDefInitialIterateAvailable=FALSE
    }else{
      X_temp <- matlabRand(num_minimizersToFind,length(X_ul_in[1,]))*(repmat(X_ul_in[2,],num_minimizersToFind,1)-repmat(X_ul_in[1,],num_minimizersToFind,1))+repmat(X_ul_in[1,],num_minimizersToFind,1)
    }

#    print(initialIterateMatrix)
#    print(matlabRand(num_minimizersToFind,length(X_ul_in[1,]))*(repmat(X_ul_in[2,],num_minimizersToFind,1)-repmat(X_ul_in[1,],num_minimizersToFind,1))+repmat(X_ul_in[1,],num_minimizersToFind,1))


    # if(algorithmVersion==4){
    #   minX=initial_lowerRange
    #   maxX=initial_upperRange
    #
    #   X_temp=round((X_temp)%*%diag(1/(maxX-minX))*numDiscretization)%*%diag((maxX-minX))/numDiscretization
    # }

    # replace the rows of X matrix with randomly sampled x if x was
    # determined to be not valid
    X[is_alive==0,]=X_temp[is_alive==0,]

    # if("parallel" %in% tolower((.packages()))){
    #
    #   Y_list=mclapply(split(X, rep(seq(1:nrow(X)),ncol(X))), tryCatch_nonlinearFunction, num_observations=num_observations, nonlinearFunction=nonlinearFunction, mc.cores = 4)
    #
    # }else{

    toCompute_X=X[is_alive==0,]


    if(sum(is_alive==0)==1){
      toCompute_X=matrix(toCompute_X,nrow=1,ncol = length(toCompute_X))
    }


    if(nonlinearFunctionCanTakeMatrixInput){
      toCompute_Y=tryCatch_nonlinearFunction(x = toCompute_X, num_observations=num_observations, nonlinearFunction=nonlinearFunction, validTarget=validTarget_vec, ... )
    }else{
      Y_list=lapply(split(toCompute_X, rep(seq(1:nrow(toCompute_X)),ncol(toCompute_X))), tryCatch_nonlinearFunction, num_observations=num_observations, nonlinearFunction=nonlinearFunction, validTarget=validTarget_vec, ... )
      toCompute_Y=t(matrix(unlist(Y_list),ncol=length(Y_list)))
    }

    Y[is_alive==0,]=toCompute_Y

    # if(nonlinearFunctionCanTakeMatrixInput){
    #   Y=tryCatch_nonlinearFunction(X, num_observations=num_observations, nonlinearFunction=nonlinearFunction, validTarget=validTarget_vec)
    # }else{
    #   Y_list=lapply(split(X, rep(seq(1:nrow(X)),ncol(X))), tryCatch_nonlinearFunction, num_observations=num_observations, nonlinearFunction=nonlinearFunction, validTarget=validTarget_vec)
    #   Y=t(matrix(unlist(Y_list),ncol=length(Y_list)))
    # }
    #}




    # compute SSR
    residual <- matlabSum((weightMatrix_forRes*(Y-targetMatrix))^2,2)

    if(showIntermetiateResults){
      message(residual)
    }

    # determine valid x (in this case if the residual is not a number
    # we consider x to be valid)
    is_alive=is.nan(residual)==0

    message(paste("Generating initial cluster.",sum(is_alive),"out of",num_minimizersToFind,"done"))

    if(sum(is_alive)==0){

      ## this error message is not tested yet
      message("It is most likely that either nonlinearfunction definition, intial range selection, or targetVector definition is incorrect. Stop the CGNM run and fix the issue.")
      message("an example of targetMatrix: ", targetMatrix[1,])
      message("an example of evaluated function: ", Y[1,])
      message("Any non-numerical values on these may give you a clue for this issue.")
      stop()
    }
  }

  # store the residual vector for the initial cluster into the
  # residual_history matrix
  residual_history <- rbind(residual_history,residual)

  # store the matrix X and Y to X_history and Y_history, respectively
  X_history[[1]] <- X
  Y_history[[1]] <- Y

  lambda_history <- t(as.matrix(lambda_vec))

  prev_residual <- residual

  for (k in seq(1,num_iteration)){
    iteration_time_start=Sys.time()

    if(algorithmParameter_gamma=="AUTO"){
      if(k<50){
        gamma_toUse=2^round(k/15)
      }else{
        gamma_toUse=2^round(50/15)
      }
    }else{
      gamma_toUse=algorithmParameter_gamma
    }

    if(algorithmVersion==3){

      out_temp <- main_iteration_version3_fast(X, Y, lambda_vec, X_ul_in[1,], X_ul_in[2,], targetMatrix, gamma_toUse, stayIn_initialRange=stayIn_initialRange, keepInitialDistribution = keepInitialDistribution, weightMatrix=weightMatrix)

      # }else if(algorithmVersion==3.1){
      #
      #   switchIndex=10
      #
      #   if(k==switchIndex){
      #     lambda_vec <- matrix(1,1,num_minimizersToFind)*algorithmParameter_initialLambda# initialise regularisation parameter lambda
      #   }
      #
      #   if(k>=switchIndex){
      #     best_index=which(residual==min(residual))
      #     bestFit=Y[best_index,]
      #     targetMatrix=repmat(bestFit,dim(X)[1],1)
      #
      #     targetMatrix[best_index,]=targetVector
      #   }
      #
      #   out_temp <- main_iteration_version3_fast(X, Y, lambda_vec, X_ul_in[1,], X_ul_in[2,], targetMatrix, gamma_toUse)
      #
      # }else if(algorithmVersion==4){
      #
      #   if(k==10){
      #     lambda_vec <- matrix(1,1,num_minimizersToFind)*algorithmParameter_initialLambda# initialise regularisation parameter lambda
      #   }
      #
      #   if(k<10){
      #     out_temp <- main_iteration_version4(X, Y, lambda_vec, X_ul_in[1,], X_ul_in[2,], targetMatrix, gamma_toUse, numDiscretization=numDiscretization)
      #   }else{
      #
      #     out_temp <- main_iteration_version3(X, Y, lambda_vec, X_ul_in[1,], X_ul_in[2,], targetMatrix, gamma_toUse)
      #   }
      #

    }else{
      out_temp <- main_iteration_version1(X, Y, lambda_vec, X_ul_in[1,], X_ul_in[2,], targetMatrix, gamma_toUse)
    }

    X_new=X+out_temp

    # if("parallel" %in% tolower((.packages()))){
    #   Y_list=mclapply(split(X_new, rep(seq(1:nrow(X_new)),ncol(X_new))), tryCatch_nonlinearFunction, num_observations=num_observations, nonlinearFunction=nonlinearFunction, mc.cores = 4)
    #
    #   }else{

    toUpdate=(lambda_vec<(algorithmParameter_initialLambda*10^10))

    X_toCompute=X_new[toUpdate,]

    if(sum(toUpdate)==1){
      X_toCompute=matrix(X_toCompute,nrow=1,ncol = dim(X_new)[2])
    }

    Y_new=Y
    if(sum(toUpdate)>0){
      if(nonlinearFunctionCanTakeMatrixInput){
        Y_computed=tryCatch_nonlinearFunction(x=X_toCompute, num_observations=num_observations, nonlinearFunction=nonlinearFunction, validTarget=validTarget_vec, ... )
      }else{
        Y_list=lapply(split(X_toCompute, rep(seq(1:nrow(X_toCompute)),ncol(X_toCompute))), tryCatch_nonlinearFunction, num_observations=num_observations, nonlinearFunction=nonlinearFunction, validTarget=validTarget_vec, ... )
        Y_computed=t(matrix(unlist(Y_list),ncol=length(Y_list)))

      }
      Y_new[toUpdate,]=Y_computed
    }else{
      break
    }


    #}

    if(showIntermetiateResults){
      message(Y_new)
    }

    residual_new <- t(matlabSum((weightMatrix_forRes*(Y_new-targetMatrix))^2,2))


    residual <- prev_residual
    for (i in seq(1,num_minimizersToFind)){
      if (!is.nan(residual_new[i])&&!is.nan(prev_residual[i])&&(prev_residual[i]>residual_new[i])){
        X[i,] <- X[i,]+out_temp[i,]
        residual[i] <- residual_new[i]
        Y[i,] <- Y_new[i,]
        lambda_vec[i] <- lambda_vec[i]/10
      } else {
        if(lambda_vec[i]<(algorithmParameter_initialLambda*10^10)){
          lambda_vec[i] <- lambda_vec[i]*10
        }
      }
    }


    X_history[[k+1]] <- X
    Y_history[[k+1]] <- Y
    lambda_history <- cbind(lambda_history,t(lambda_vec))

    residual_history <- cbind(residual_history,residual)

    prev_residual <- residual

    # message(matlabMedian(residual_history,1))
    message(paste0("Iteration:",toString(k),"  Median sum of squares residual=", toString(median(prev_residual))))


    if(saveLog){

      CGNM_result=list(X=X,Y=Y,residual_history=residual_history, initialX=X_history[[1]], initialY=Y_history[[1]],lambda_history=lambda_history, runSetting=runSetting )
      save(file=paste0(saveFolderName,'/iteration_',toString(k),'.RDATA'), CGNM_result)
    }
    iteration_time_end=Sys.time()

    if(k==1){
      message(paste("Rough estimation of remaining computation time:",round(as.numeric(iteration_time_end-iteration_time_start, units="mins")*(num_iteration-k)*10)/10,"min"))
    }

    if(k%%10==1){
      message(paste("CGNM iteration estimated to finish at:",(iteration_time_end-iteration_time_start)*(num_iteration-k)+Sys.time()))
    }

    if(median(prev_residual)==0){
      break
    }
  }
  CGNM_result=list(X=X,Y=Y,residual_history=residual_history, lambda_history=lambda_history, runSetting=runSetting, initialX=X_history[[1]], initialY=Y_history[[1]] )



  if(saveLog){
    save(file=paste0(saveFolderName,'/iteration_final.RDATA'), CGNM_result)
  }

  CGNM_end_time=Sys.time()

  message(paste("CGNM computation time: ",round(as.numeric(CGNM_end_time-CGNM_start_time, units="mins")*10)/10,"min"))

  return(CGNM_result)

}




main_iteration_version4 <-  function(X_in, Y_in, lambdaV, minX, maxX, targetMatrix, algorithmParameter_gamma, numDiscretization=100){

  #// aliveIndex: index of the parameter set in the cluster whose SSR have decreased in the previous iteration.
  #// X_in: all set of parameters in the cluster from the previous itaration  X_in[i][j]: jth parameter value in the ith parameter set in the cluster
  #// Y_in: all set of the function value (solution of the forward problem, f(parameter)) from the previous iteration
  # global algorithmParameter_gamma

  ##//STEP 1: Linear Approximation

  YisNotNaN <- matrix(1,dim(X_in)[1],1)

  for (k in 1:dim(Y_in)[1]){
    for (i in 1:dim(Y_in)[2]){
      if (is.nan(Y_in[k,i])|is.infinite(Y_in[k,i])){
        YisNotNaN[k] <- 0
      }
    }
  }

  #// Linear approximation will be A x + y_o \approx f(x) = y for all y in Y matrix (stored in each row of Y matrix)
  #// This can be written as A_with_y0 X_append_one_column^T \approx Y^T


  kmean_result=optimal_kmeans(X_in)

  notAcceptableClusters=as.numeric(which(table(kmean_result$cluster)==1,arr.ind = TRUE))

  relative_distance <- matrix(0,dim(X_in)[1],dim(X_in)[1])
  rec_relative_distance <- matrix(0,dim(X_in)[1],dim(X_in)[1])

  for (i in seq(1,dim(X_in)[1])){
    relative_distance[i,i] <- 0
    rec_relative_distance[i,i] <- 0

    if(i>1){
      for (j in seq(1,(i-1))){
        relative_distance[i,j] <- 0
        relative_distance[j,i] <- 0

        if(kmean_result$cluster[i]==kmean_result$cluster[j]){

          for (k in seq(1,dim(X_in)[2])){
            relative_distance[i,j] <- relative_distance[i,j]+((X_in[i,k]-X_in[j,k])/(sd(X_in[kmean_result$cluster==kmean_result$cluster[j],k])))^2
          }
        }else{
          relative_distance[i,j] =0
        }

        if (YisNotNaN[i]==1&YisNotNaN[j]==1&relative_distance[i,j]!=0){

          relative_distance[j,i] <- relative_distance[i,j]
          rec_relative_distance[i,j] <- 1/relative_distance[i,j]
          rec_relative_distance[j,i] <- rec_relative_distance[i,j]
        } else {
          relative_distance[i,j] <- 0
          rec_relative_distance[i,j] <- 0
          rec_relative_distance[j,i] <- 0
        }

      }

    }

  }




  # Linear approximation will be A_with_y0  X_append_one_column^T \approx Y^T now solve it for A_with_y0
  # We solve for \min || X_append_one_column  A_with_y0^T - Y ||_F
  # use CGNR and solve each row of A_with_y0 (i.e., each column of A_with_y0^T)
  A <- matrix(0, dim(Y_in)[2],dim(X_in)[2])
  deltaX <- matrix(0, dim(X_in)[1],dim(X_in)[2])
  delta_X <- matrix(0, dim(X_in)[1],dim(X_in)[2])
  tempOnes <- matrix(1, dim(Y_in)[1],1)

  for (k in seq(1,dim(Y_in)[1])){
    if(lambdaV[k]<10^10&(!(kmean_result$cluster[k] %in% notAcceptableClusters))){
      ## 2-1): Construct weighted linear approximation of the nonlinear function
      for (i in seq(1,dim(X_in)[1])){
        deltaX[i,] <- X_in[i,]-X_in[k,]
      }

      # use the following code if one wishes to use MASS's matrix
      # inverse function (ginv) instead of CGNR

      rec_relative_distance_tempvec=(rec_relative_distance[ ,k])^algorithmParameter_gamma
      rec_relative_distance_tempvec[is.infinite(rec_relative_distance_tempvec)]=max(rec_relative_distance_tempvec[!is.infinite(rec_relative_distance_tempvec)])

      AT=t(diag(rec_relative_distance_tempvec)%*%deltaX)

      ATAinv=MASS::ginv(AT%*%t(AT))


      delta_Y=Y_in

      for(i in seq(1, dim(Y_in)[2])){
        delta_Y[,i]=Y_in[,i]-Y_in[k,i]
      }

      A=t(ATAinv%*%(AT%*%diag((rec_relative_distance[ ,k])^algorithmParameter_gamma)%*%(delta_Y)))


      #  for (i in seq(1,dim(Y_in)[2])){
      #    A[i,] <- CGNR_ATAx_ATb_with_weight(deltaX,Y_in[ ,i]-Y_in[k,i]*tempOnes,rec_relative_distance[ ,k]^algorithmParameter_gamma)
      #  }


      ## 2-2): Solve for x that minimizes the resodual using the weighted linear approximation
      delta_X[k,] <- CGNR_ATAx_ATb_with_reg(A, t(targetMatrix[k,]-Y_in[k,]),lambdaV[k])
    }

  }


  delta_X=round((delta_X)%*%diag(1/(maxX-minX))*numDiscretization)%*%diag((maxX-minX))/numDiscretization

  return(delta_X)
}




main_iteration_version3 <-  function(X_in, Y_in, lambdaV, minX, maxX, targetMatrix, algorithmParameter_gamma){

  #// aliveIndex: index of the parameter set in the cluster whose SSR have decreased in the previous iteration.
  #// X_in: all set of parameters in the cluster from the previous itaration  X_in[i][j]: jth parameter value in the ith parameter set in the cluster
  #// Y_in: all set of the function value (solution of the forward problem, f(parameter)) from the previous iteration
  # global algorithmParameter_gamma


  ##//STEP 1: Linear Approximation

  YisNotNaN <- matrix(1,dim(X_in)[1],1)

  for (k in 1:dim(Y_in)[1]){
    for (i in 1:dim(Y_in)[2]){
      if (is.nan(Y_in[k,i])|is.infinite(Y_in[k,i])){
        YisNotNaN[k] <- 0
      }
    }
  }

  #// Linear approximation will be A x + y_o \approx f(x) = y for all y in Y matrix (stored in each row of Y matrix)
  #// This can be written as A_with_y0 X_append_one_column^T \approx Y^T


  X_in_forkmeans=X_in

  kmean_result=optimal_kmeans(X_in_forkmeans)

  while(min(table(kmean_result$cluster))<dim(X_in)[2]){
    kmean_result=kmeans(col_normalize(X_in_forkmeans),centers = length(unique(kmean_result$cluster))-1)
  }

  notAcceptableClusters=as.numeric(which(table(kmean_result$cluster)==1,arr.ind = TRUE))

  relative_distance <- matrix(0,dim(X_in)[1],dim(X_in)[1])
  rec_relative_distance <- matrix(0,dim(X_in)[1],dim(X_in)[1])

  for (i in seq(1,dim(X_in)[1])){
    relative_distance[i,i] <- 0
    rec_relative_distance[i,i] <- 0

    if(i>1){
      for (j in seq(1,(i-1))){
        relative_distance[i,j] <- 0
        relative_distance[j,i] <- 0

        if(kmean_result$cluster[i]==kmean_result$cluster[j]){

          for (k in seq(1,dim(X_in)[2])){
            relative_distance[i,j] <- relative_distance[i,j]+((X_in[i,k]-X_in[j,k])/(sd(X_in[kmean_result$cluster==kmean_result$cluster[j],k])))^2
          }
        }else{
          relative_distance[i,j] =0
        }

        if (YisNotNaN[i]==1&YisNotNaN[j]==1&relative_distance[i,j]!=0){

          relative_distance[j,i] <- relative_distance[i,j]
          rec_relative_distance[i,j] <- 1/relative_distance[i,j]
          rec_relative_distance[j,i] <- rec_relative_distance[i,j]
        } else {
          relative_distance[i,j] <- 0
          rec_relative_distance[i,j] <- 0
          rec_relative_distance[j,i] <- 0
        }

      }

    }

  }




  # Linear approximation will be A_with_y0  X_append_one_column^T \approx Y^T now solve it for A_with_y0
  # We solve for \min || X_append_one_column  A_with_y0^T - Y ||_F
  # use CGNR and solve each row of A_with_y0 (i.e., each column of A_with_y0^T)
  A <- matrix(0, dim(Y_in)[2],dim(X_in)[2])
  deltaX <- matrix(0, dim(X_in)[1],dim(X_in)[2])
  delta_X <- matrix(0, dim(X_in)[1],dim(X_in)[2])
  tempOnes <- matrix(1, dim(Y_in)[1],1)

  for (k in seq(1,dim(Y_in)[1])){
    if(lambdaV[k]<10^10&(!(kmean_result$cluster[k] %in% notAcceptableClusters))){
      ## 2-1): Construct weighted linear approximation of the nonlinear function
      for (i in seq(1,dim(X_in)[1])){
        deltaX[i,] <- X_in[i,]-X_in[k,]
      }

      # use the following code if one wishes to use MASS's matrix
      # inverse function (ginv) instead of CGNR



      rec_relative_distance_tempvec=(rec_relative_distance[ ,k])^algorithmParameter_gamma
      rec_relative_distance_tempvec[is.infinite(rec_relative_distance_tempvec)]=max(rec_relative_distance_tempvec[!is.infinite(rec_relative_distance_tempvec)])

      AT=t(diag(rec_relative_distance_tempvec)%*%deltaX)

      ATAinv=MASS::ginv(AT%*%t(AT))

      delta_Y=Y_in

      for(i in seq(1, dim(Y_in)[2])){
        delta_Y[,i]=Y_in[,i]-Y_in[k,i]
      }

      A=t(ATAinv%*%(AT%*%diag((rec_relative_distance[ ,k])^algorithmParameter_gamma)%*%(delta_Y)))

      #  for (i in seq(1,dim(Y_in)[2])){
      #    A[i,] <- CGNR_ATAx_ATb_with_weight(deltaX,Y_in[ ,i]-Y_in[k,i]*tempOnes,rec_relative_distance[ ,k]^algorithmParameter_gamma)
      #  }


      ## 2-2): Solve for x that minimizes the resodual using the weighted linear approximation
      delta_X[k,] <- CGNR_ATAx_ATb_with_reg(A, t(targetMatrix[k,]-Y_in[k,]),lambdaV[k])
    }

  }

  return(delta_X)
}


main_iteration_version1 <-  function(X_in, Y_in, lambdaV, minX, maxX, targetMatrix, algorithmParameter_gamma){

  #// aliveIndex: index of the parameter set in the cluster whose SSR have decreased in the previous iteration.
  #// X_in: all set of parameters in the cluster from the previous itaration  X_in[i][j]: jth parameter value in the ith parameter set in the cluster
  #// Y_in: all set of the function value (solution of the forward problem, f(parameter)) from the previous iteration
  # global algorithmParameter_gamma


  ##//STEP 1: Linear Approximation

  YisNotNaN <- matrix(1,dim(X_in)[1],1)

  for (k in 1:dim(Y_in)[1]){
    for (i in 1:dim(Y_in)[2]){
      if (is.nan(Y_in[k,i])|is.infinite(Y_in[k,i])){
        YisNotNaN[k] <- 0
      }
    }
  }

  #// Linear approximation will be A x + y_o \approx f(x) = y for all y in Y matrix (stored in each row of Y matrix)
  #// This can be written as A_with_y0 X_append_one_column^T \approx Y^T


  x_width <- maxX-minX

  relative_distance <- matrix(0,dim(X_in)[1],dim(X_in)[1])
  rec_relative_distance <- matrix(0,dim(X_in)[1],dim(X_in)[1])

  for (i in seq(1,dim(X_in)[1])){
    relative_distance[i,i] <- 0
    rec_relative_distance[i,i] <- 0

    if(i>1){
      for (j in seq(1,(i-1))){
        relative_distance[i,j] <- 0
        relative_distance[j,i] <- 0

        for (k in seq(1,dim(X_in)[2])){
          relative_distance[i,j] <- relative_distance[i,j]+((X_in[i,k]-X_in[j,k])/x_width[k])^2
        }

        if (YisNotNaN[i]==1&YisNotNaN[j]==1){

          relative_distance[i,j] <- (relative_distance[i,j])

          relative_distance[j,i] <- relative_distance[i,j]
          rec_relative_distance[i,j] <- 1/relative_distance[i,j]
          rec_relative_distance[j,i] <- rec_relative_distance[i,j]
        } else {
          relative_distance[i,j] <- 0
          rec_relative_distance[i,j] <- 0
          rec_relative_distance[j,i] <- 0
        }

      }

    }

  }




  # Linear approximation will be A_with_y0  X_append_one_column^T \approx Y^T now solve it for A_with_y0
  # We solve for \min || X_append_one_column  A_with_y0^T - Y ||_F
  # use CGNR and solve each row of A_with_y0 (i.e., each column of A_with_y0^T)
  A <- matrix(0, dim(Y_in)[2],dim(X_in)[2])
  deltaX <- matrix(0, dim(X_in)[1],dim(X_in)[2])
  delta_X <- matrix(0, dim(X_in)[1],dim(X_in)[2])
  tempOnes <- matrix(1, dim(Y_in)[1],1)

  for (k in seq(1,dim(Y_in)[1])){

    if(lambdaV[k]<10^10){
      ## 2-1): Construct weighted linear approximation of the nonlinear function
      for (i in seq(1,dim(X_in)[1])){
        deltaX[i,] <- X_in[i,]-X_in[k,]
      }

      # use the following code if one wishes to use matlab's matrix
      # inverse function instead of CGNR
      #        AT=(diag((rec_relative_distance[ ,k]).^algorithmParameter_gamma)*deltaX)';
      #        ATAinv=inv(AT*AT');
      #        A=(ATAinv*(AT*diag((rec_relative_distance[ ,k]).^algorithmParameter_gamma)*(Y_in-repmat(Y_in[k,],dim(Y_in)[1],1))))';

      for (i in seq(1,dim(Y_in)[2])){
        A[i,] <- CGNR_ATAx_ATb_with_weight(deltaX,Y_in[ ,i]-Y_in[k,i]*tempOnes,rec_relative_distance[ ,k]^algorithmParameter_gamma)
      }

      ## 2-2): Solve for x that minimizes the resodual using the weighted linear approximation
      delta_X[k,] <- CGNR_ATAx_ATb_with_reg(A, t(targetMatrix[k,]-Y_in[k,]),lambdaV[k])
    }

  }

  return(delta_X)
}


main_iteration_version3_fast <-  function(X_in_original, Y_in, lambdaV_original, minX_original, maxX_original, targetMatrix, algorithmParameter_gamma, stayIn_initialRange=FALSE, keepInitialDistribution=NA, weightMatrix=NA){


  validXcol=!(col_max(X_in_original)==col_min(X_in_original))

  X_in=X_in_original[,validXcol]
  lambdaV=lambdaV_original
  minX=minX_original[validXcol]
  maxX=maxX_original[validXcol]


  #// aliveIndex: index of the parameter set in the cluster whose SSR have decreased in the previous iteration.
  #// X_in: all set of parameters in the cluster from the previous itaration  X_in[i][j]: jth parameter value in the ith parameter set in the cluster
  #// Y_in: all set of the function value (solution of the forward problem, f(parameter)) from the previous iteration
  # global algorithmParameter_gamma


  ##//STEP 1: Linear Approximation

  YisNotNaN <- matrix(1,dim(X_in)[1],1)

  for (k in 1:dim(Y_in)[1]){
    for (i in 1:dim(Y_in)[2]){
      if (is.nan(Y_in[k,i])|is.infinite(Y_in[k,i])){
        YisNotNaN[k] <- 0
      }
    }
  }

  #// Linear approximation will be A x + y_o \approx f(x) = y for all y in Y matrix (stored in each row of Y matrix)
  #// This can be written as A_with_y0 X_append_one_column^T \approx Y^T

  kmean_result=optimal_kmeans(X_in)

  while(min(table(kmean_result$cluster))<dim(X_in)[2]|min(kmean_result$withinss)==0){
    kmean_result=kmeans(col_normalize(X_in),centers = length(unique(kmean_result$cluster))-1)
  }


  notAcceptableClusters=as.numeric(which(table(kmean_result$cluster)==1,arr.ind = TRUE))

  relative_distance <- matrix(0,dim(X_in)[1],dim(X_in)[1])
  rec_relative_distance <- matrix(0,dim(X_in)[1],dim(X_in)[1])


  numClusters=length(unique(kmean_result$cluster))
  sdMatrix=matrix(0,  dim(X_in)[2],dim(X_in)[1])

  for(i in seq(1, numClusters)){
    for(j in seq(1, dim(X_in)[2])){
      sdMatrix[j,kmean_result$cluster==i]=sd(X_in[kmean_result$cluster==i,j])
    }
  }


  tX_in=t(X_in)

  for (i in seq(1,dim(X_in)[1])){

    useIndex=(kmean_result$cluster==kmean_result$cluster[i])
    sdNum=sdMatrix[,useIndex]
    sdNum=ifelse(sdNum<10^-9,10^-9, sdNum)

    relative_distance[i,useIndex] = colSums(((tX_in[,useIndex]-tX_in[,i])/(sdNum))^2)

  }


  rec_relative_distance=1/relative_distance
  rec_relative_distance[relative_distance==0]=0



  # Linear approximation will be A_with_y0  X_append_one_column^T \approx Y^T now solve it for A_with_y0
  # We solve for \min || X_append_one_column  A_with_y0^T - Y ||_F
  # use CGNR and solve each row of A_with_y0 (i.e., each column of A_with_y0^T)
  A <- matrix(0, dim(Y_in)[2],dim(X_in)[2])
  deltaX <- matrix(0, dim(X_in)[1],dim(X_in)[2])
  delta_X <- matrix(0, dim(X_in)[1],dim(X_in)[2])
  tempOnes <- matrix(1, dim(Y_in)[1],1)

  for (k in seq(1,dim(Y_in)[1])){
    if(lambdaV[k]<10^10&(!(kmean_result$cluster[k] %in% notAcceptableClusters))){
      ## 2-1): Construct weighted linear approximation of the nonlinear function
      for (i in seq(1,dim(X_in)[1])){
        deltaX[i,] <- X_in[i,]-X_in[k,]
      }

      # use the following code if one wishes to use MASS's matrix
      # inverse function (ginv) instead of CGNR

      A <- tryCatch(
        {
          rec_relative_distance_tempvec=(rec_relative_distance[ ,k])^algorithmParameter_gamma
          rec_relative_distance_tempvec[is.infinite(rec_relative_distance_tempvec)]=max(c(10^10,rec_relative_distance_tempvec[!is.infinite(rec_relative_distance_tempvec)]))

          AT=t(diag(rec_relative_distance_tempvec)%*%deltaX)

          ATAinv=MASS::ginv(AT%*%t(AT))

          delta_Y=Y_in

          for(i in seq(1, dim(Y_in)[2])){
            delta_Y[,i]=Y_in[,i]-Y_in[k,i]
          }

          A=t(ATAinv%*%(AT%*%diag((rec_relative_distance[ ,k])^algorithmParameter_gamma)%*%(delta_Y)));

          if(!is.na(keepInitialDistribution[1])&&length(keepInitialDistribution)==dim(X_in)[2]){
            A[,keepInitialDistribution]=0
          }

          (A)

        },
        error=function(cond) {

          print("using CGNR")
          for (i in seq(1,dim(Y_in)[2])){
            A[i,] <- CGNR_ATAx_ATb_with_weight(deltaX,Y_in[ ,i]-Y_in[k,i]*tempOnes,rec_relative_distance[ ,k]^algorithmParameter_gamma)
          }

          (A)
        },
        warning=function(cond) {

          rec_relative_distance_tempvec=(rec_relative_distance[ ,k])^algorithmParameter_gamma
          rec_relative_distance_tempvec[is.infinite(rec_relative_distance_tempvec)]=max(c(10^10,rec_relative_distance_tempvec[!is.infinite(rec_relative_distance_tempvec)]))

          AT=t(diag(rec_relative_distance_tempvec)%*%deltaX)

          ATAinv=MASS::ginv(AT%*%t(AT))

          delta_Y=Y_in

          for(i in seq(1, dim(Y_in)[2])){
            delta_Y[,i]=Y_in[,i]-Y_in[k,i]
          }

          A=t(ATAinv%*%(AT%*%diag((rec_relative_distance[ ,k])^algorithmParameter_gamma)%*%(delta_Y)));

          if(!is.na(keepInitialDistribution)&&length(keepInitialDistribution)==dim(X_in)[2]){

            A[,keepInitialDistribution]=0
          }

          (A)
        },
        finally={

        }
      )


      #
      # <<<<<<< Updated upstream
      # =======
      #       AT=t(diag(rec_relative_distance_tempvec)%*%deltaX)
      #
      #       ATAinv=MASS::ginv(AT%*%t(AT))
      #
      #       delta_Y=Y_in
      #
      #       for(i in seq(1, dim(Y_in)[2])){
      #         delta_Y[,i]=Y_in[,i]-Y_in[k,i]
      #       }
      #
      #       A=t(ATAinv%*%(AT%*%diag((rec_relative_distance[ ,k])^algorithmParameter_gamma)%*%(delta_Y)));
      #
      #       if(!is.na(keepInitialDistribution)[1]&&length(keepInitialDistribution)==dim(X_in)[2]){
      #         A[,keepInitialDistribution]=0
      #       }
      #
      #       #  for (i in seq(1,dim(Y_in)[2])){
      #       #    A[i,] <- CGNR_ATAx_ATb_with_weight(deltaX,Y_in[ ,i]-Y_in[k,i]*tempOnes,rec_relative_distance[ ,k]^algorithmParameter_gamma)
      #       #  }
      # >>>>>>> Stashed changes


      ## 2-2): Solve for x that minimizes the resodual using the weighted linear approximation
      if(!is.matrix(weightMatrix)){
        delta_X[k,] <- CGNR_ATAx_ATb_with_reg(A, t(targetMatrix[k,]-Y_in[k,]),lambdaV[k])
      }else if(dim(weightMatrix)[1]==dim(Y_in)[1]&dim(weightMatrix)[2]==dim(Y_in)[2]){
        delta_X[k,] <- CGNR_ATAx_ATb_with_reg_weight(A, t(targetMatrix[k,]-Y_in[k,]),lambdaV[k],weightMatrix[k,])
      }else{
        message("Invalid weight matrix is provided hence we will now use uniform weight of one.")
        delta_X[k,] <- CGNR_ATAx_ATb_with_reg(A, t(targetMatrix[k,]-Y_in[k,]),lambdaV[k])

      }

      if(stayIn_initialRange){
        scaling_vector=rep(1,dim(delta_X)[2])
        scaling_vector[((delta_X[k,])/(minX-X_in[k,]))>1]=((delta_X[k,])/(minX-X_in[k,])*2)[((delta_X[k,])/(minX-X_in[k,]))>1]
        scaling_vector[((delta_X[k,])/(maxX-X_in[k,]))>1]=((delta_X[k,])/(maxX-X_in[k,])*2)[((delta_X[k,])/(maxX-X_in[k,]))>1]
        #scaling_vector[((delta_X[k,])/(minX-X_in[k,]))>1]=((delta_X[k,])/(minX-X_in[k,])*1.0001)[((delta_X[k,])/(minX-X_in[k,]))>1]
        #scaling_vector[((delta_X[k,])/(maxX-X_in[k,]))>1]=((delta_X[k,])/(maxX-X_in[k,])*1.0001)[((delta_X[k,])/(maxX-X_in[k,]))>1]
        delta_X[k,]=delta_X[k,]/scaling_vector
      }
    }

  }

  delta_X_out=matrix(0, ncol = dim(X_in_original)[2], nrow = dim(X_in_original)[1])
  delta_X_out[,validXcol]=delta_X

  return(delta_X_out)
}

CGNR_ATAx_ATb <- function(A, b){#solves A^TA x = A^Tb

  num_max_iter <- dim(A)[2]


  x <- matrix(0,dim(A)[2],1)

  if (dim(A)[1]==length(b)){

    r <- t(A)%*%b
    p <- r

    innerProd_rr_old <- dot(r,r)

    for (i in seq(1,num_max_iter)){
      tempVector <- A%*%p
      tempD <- dot(tempVector,tempVector)
      if (is.na(tempD)||tempD==0){
        break
      }
      alpha <- innerProd_rr_old/tempD
      x <- x+alpha*p
      r <- r-alpha*t(A)%*%tempVector
      innerProd_rr_new <- dot(r,r)
      beta <- innerProd_rr_new/innerProd_rr_old
      innerProd_rr_old <- innerProd_rr_new
      if (innerProd_rr_old==0){
        break
      }
      p <- r+beta*p
    }
  } else {
    message('CGNR error: matrix A and vetor b sizes do not make sense.')
    message(paste('Size of A: ', toString(dim(A)[1])))
    message('Size of b: ', toString(length(b)))
  }
  return(x)
}

CGNR_ATAx_ATb_with_weight=function(A, b, weight){
  for (i in seq(1,min(length(weight),dim(A)[1]))){
    if (weight[i]==0){
      for (j in 1:dim(A)[2]){
        A[i,j] <- 0
      }
      b[i] <- 0
    } else {
      A[i,] <- weight[i]*A[i,]
      b[i] <- weight[i]*b[i]
    }
  }
  return(CGNR_ATAx_ATb(A[weight!=0,],b[weight!=0]))
}


CGNR_ATAx_ATb_with_reg_weight=function(A, b, lambda, weight){

  for (i in seq(1,min(length(weight),dim(A)[1]))){
    if (weight[i]==0){
      for (j in 1:dim(A)[2]){
        A[i,j] <- 0
      }
      b[i] <- 0
    } else {
      A[i,] <- weight[i]*A[i,]
      b[i] <- weight[i]*b[i]
    }
  }
  return(CGNR_ATAx_ATb_with_reg(A[weight!=0,],t(b[weight!=0]),lambda))
}


CGNR_ATAx_ATb_with_reg=function(A, b, lambda){#   //solves (A^TA+lambda I) x = A^Tb

  num_max_iter <- dim(A)[1]

  x <- matrix(0,dim(A)[2],1)

  if (dim(A)[1]==length(b)){

    r <- t(A)%*%t(b)
    p <- r

    innerProd_rr_old <- dot(r,r)
    AtA <- t(A)%*%A
    AtA_lambda <- AtA

    for (i in seq(1,dim(AtA)[1])){
      AtA_lambda[i,i] <- AtA_lambda[i,i]+lambda
    }


    for (i in seq(1,num_max_iter)){
      tempD <- dot(p,(AtA_lambda%*%p))
      if (is.na(tempD)||tempD==0){
        break
      }
      alpha <- innerProd_rr_old/tempD
      x <- x+alpha*p
      r <- r-alpha*(AtA_lambda%*%p)
      innerProd_rr_new <- dot(r,r)
      beta <- innerProd_rr_new/innerProd_rr_old
      innerProd_rr_old <- innerProd_rr_new
      if (innerProd_rr_old==0){
        break
      }
      p <- r+beta*p
    }
  } else {
    message('CGNR with reg error: matrix A and vetor b sizes do not make sense.')
    message(paste('Size of A: ', toString(dim(A)[1])))
    message(paste('Size of b: ', toString(length(b))))
  }
  return(x)
}



CGNR_ATAx_ATb_with_regVec=function(A, b, lambda){#   //solves (A^TA+lambda_vec) x = A^Tb

  num_max_iter <- dim(A)[1]

  x <- matrix(0,dim(A)[2],1)

  if (dim(A)[1]==length(b)){

    r <- t(A)%*%t(b)
    p <- r

    innerProd_rr_old <- dot(r,r)
    AtA <- t(A)%*%A
    AtA_lambda <- AtA

    for (i in seq(1,dim(AtA)[1])){
      AtA_lambda[i,i] <- AtA_lambda[i,i]+lambda[i]
    }


    for (i in seq(1,num_max_iter)){
      tempD <- dot(p,(AtA_lambda%*%p))
      if (is.na(tempD)||tempD==0){
        break
      }
      alpha <- innerProd_rr_old/tempD
      x <- x+alpha*p
      r <- r-alpha*(AtA_lambda%*%p)
      innerProd_rr_new <- dot(r,r)
      beta <- innerProd_rr_new/innerProd_rr_old
      innerProd_rr_old <- innerProd_rr_new
      if (innerProd_rr_old==0){
        break
      }
      p <- r+beta*p
    }
  } else {
    message('CGNR with reg error: matrix A and vetor b sizes do not make sense.')
    message(paste('Size of A: ', toString(dim(A)[1])))
    message(paste('Size of b: ', toString(length(b))))
  }
  return(x)
}



#' @title Cluster_Gauss_Newton_Bootstrap_method
#' @description Conduct residual resampling bootstrap analyses using CGNM.
#' @param CGNM_result (required input) \emph{A list} stores the computational result from Cluster_Gauss_Newton_method() function in CGNM package.
#' @param nonlinearFunction (required input) \emph{A function with input of a vector x of real number of length n and output a vector y of real number of length m.} In the context of model fitting the nonlinearFunction is \strong{the model}.  Given the CGNM does not assume the uniqueness of the minimizer, m can be less than n.  Also CGNM does not assume any particular form of the nonlinear function and also does not require the function to be continuously differentiable (see Appendix D of our publication for an example when this function is discontinuous).
#' @param num_bootstrapSample (default: 200) \emph{A positive integer} number of bootstrap samples to generate.
#' @param indicesToUseAsInitialIterates (default: NA) \emph{A vector of integers} indices to use for initial iterate of the bootstrap analyses.  For CGNM bootstrap, we use the parameters found by CGNM as the initial iterates, here you can manually spccify which of the approximate minimizers that was found by CGNM (where the CGNM computation result is given as CGNM_result file) to use as initial iterates.  (if NA, use indices chosen by the acceptedIndices() function with default setting).
#' @param bootstrapType (default:1) \emph{1 or 2} 1: residual resampling bootstrap method, 2: case sampling bootstrap method
#' @param ... Further arguments to be supplied to nonlinearFunction
#' @return list of a matrix X, Y,residual_history, initialX, bootstrapX, bootstrapY as well as a list runSetting.
#' \enumerate{\item X, Y, residual_history, initialX: identical to what was given as CGNM_result.
#' \item X: \emph{a num_bootstrapSample by n matrix} which stores the the X values that was sampled using residual resampling bootstrap analyses (In terms of model fitting this is the parameter combinations with variabilities that represent \strong{parameter estimation uncertainties}.).
#' \item Y: \emph{a num_bootstrapSample by m matrix} which stores the nonlinearFunction evaluated at the corresponding bootstrap analyses results in matrix bootstrapX above. In the context of model fitting each row corresponds to \strong{the model simulations}.
#' \item runSetting: identical to what is given as CGNM_result but in addition including num_bootstrapSample and indicesToUseAsInitialIterates.}
#' @examples
#' ##lip-flop kinetics (an example known to have two distinct solutions)
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
#' targetVector = observation, num_iteration = 10, num_minimizersToFind = 100,
#' initial_lowerRange = c(0.1,0.1,0.1), initial_upperRange =  c(10,10,10),
#' lowerBound=rep(0,3), ParameterNames=c("Ka","V1","CL_2"), saveLog = FALSE)
#'
#' CGNM_bootstrap=Cluster_Gauss_Newton_Bootstrap_method(CGNM_result,
#'      nonlinearFunction=model_analytic_function, num_bootstrapSample=100)
#'
#' plot_paraDistribution_byHistogram(CGNM_bootstrap)
#'
#' @export
#' @import stats MASS

Cluster_Gauss_Newton_Bootstrap_method <- function(CGNM_result, nonlinearFunction, num_bootstrapSample=200, indicesToUseAsInitialIterates=NA, bootstrapType=1, ...){

  cutoff_pvalue=0.05
  numParametersIncluded=NA
  useAcceptedApproximateMinimizers=TRUE

  targetVector=CGNM_result$runSetting$targetVector
  initial_lowerRange=CGNM_result$runSetting$initial_lowerRange
  initial_upperRange=CGNM_result$runSetting$initial_upperRange
  num_minimizersToFind=CGNM_result$runSetting$num_minimizersToFind
  num_iteration=CGNM_result$runSetting$num_iteration
  saveLog=CGNM_result$runSetting$saveLog
  runName=CGNM_result$runSetting$runName
  textMemo=CGNM_result$runSetting$textMemo
  algorithmParameter_initialLambda=CGNM_result$runSetting$algorithmParameter_initialLambda
  algorithmParameter_gamma=CGNM_result$runSetting$algorithmParameter_gamma
  algorithmVersion=CGNM_result$runSetting$algorithmVersion
  ReparameterizationDef=CGNM_result$runSetting$ReparameterizationDef
  ParameterTransFormationDef=CGNM_result$runSetting$ParameterTransFormationDef
  ParameterNames=CGNM_result$runSetting$ParameterNames

  if(is.null(CGNM_result$runSetting$MO_weights)){
    MO_weights=rep(0, length(initial_lowerRange))
    MO_values=(initial_lowerRange+initial_upperRange)/2
  }else{
    MO_values=CGNM_result$runSetting$MO_values
    MO_weights=CGNM_result$runSetting$MO_weights
  }
  if(length(MO_weights)==1&is.na(MO_weights)[1]){
    MO_weights=rep(0, length(initial_lowerRange))
    MO_values=(initial_lowerRange+initial_upperRange)/2
  }

  MO_index=(MO_weights!=0)
  MO_values[!MO_index]=(initial_upperRange[!MO_index]+initial_lowerRange[!MO_index])/2




  #  CGNM_result=Cluster_Gauss_Newton_method(nonlinearFunction, targetVector, initial_lowerRange, initial_upperRange , num_minimizersToFind, num_iteration, saveLog, runName, textMemo,algorithmParameter_initialLambda, algorithmParameter_gamma, algorithmVersion, initialIterateMatrix)
  if(is.na(indicesToUseAsInitialIterates[1])){
    acceptParaIndex=acceptedIndices(CGNM_result, cutoff_pvalue, numParametersIncluded, useAcceptedApproximateMinimizers)
  }else{
    acceptParaIndex=indicesToUseAsInitialIterates[indicesToUseAsInitialIterates<=dim(CGNM_result$X)[1]]
  }

  if(length(acceptParaIndex)<100){
    message(paste("WARNING: only", length(acceptParaIndex) ,"acceptable parameters found, we will run bootstrap with top 50 parameters."))
    #if(length(acceptParaIndex)<10){
    #  stop("Too few accepted parameters.")
    #}
    acceptParaIndex=acceptedIndices(CGNM_result, cutoff_pvalue, 100, FALSE)

  }


  bootstrapTargetMatrix=matrix(NA, nrow=num_bootstrapSample, ncol = dim(CGNM_result$Y)[2])
  bootstrapInitialMatrix=matrix(NA, nrow=num_bootstrapSample, ncol = dim(CGNM_result$X)[2])

  residual_vec=targetVector-CGNM_result$Y[acceptParaIndex[1],]
  temp_residual_vec=targetVector-CGNM_result$Y[acceptParaIndex[1],]
  residual_vec=temp_residual_vec[!is.na(temp_residual_vec)]

  for(i in seq(1,num_bootstrapSample)){
    bootstrapTargetMatrix[i,!is.na(temp_residual_vec)]=targetVector[!is.na(temp_residual_vec)]+sample(residual_vec,length(residual_vec), replace = TRUE)
  }

  bootstrapInitialMatrix=CGNM_result$X[sample(acceptParaIndex,num_bootstrapSample,replace=TRUE),]

  stayIn_initRange=FALSE

  if(!is.null(CGNM_result$runSetting$stayIn_initialRange)){
    stayIn_initRange=CGNM_result$runSetting$stayIn_initialRange
  }

  lowerBound=NA
  if(!is.null(CGNM_result$runSetting$lowerBound)){
    lowerBound=CGNM_result$runSetting$lowerBound
  }


  upperBound=NA
  if(!is.null(CGNM_result$runSetting$upperBound)){
    upperBound=CGNM_result$runSetting$upperBound
  }

  ParameterNames=NA
  if(!is.null(CGNM_result$runSetting$ParameterNames)){
    ParameterNames=CGNM_result$runSetting$ParameterNames
  }

  LB_index=(!is.na(as.numeric(lowerBound)))&is.na(upperBound)
  UB_index=(!is.na(as.numeric(upperBound)))&is.na(lowerBound)
  BB_index=(!is.na(as.numeric(lowerBound)))&(!is.na(as.numeric(upperBound)))

  lowerBoundMatrix=repmat(lowerBound[BB_index],numrows = num_minimizersToFind ,numcol = 1)
  upperBoundMatrix=repmat(upperBound[BB_index],numrows = num_minimizersToFind ,numcol = 1)



  MO_val_converted=MO_values
  MO_val_converted[LB_index]=log(MO_val_converted[LB_index]-lowerBound[LB_index])
  MO_val_converted[UB_index]=log(-MO_val_converted[UB_index]+upperBound[UB_index])
  MO_val_converted[BB_index]=log((MO_val_converted[BB_index]-lowerBound[BB_index])/(upperBound[BB_index]-MO_val_converted[BB_index]))
  #
  #
  # nonlinearFunction_varTrans=function(x, ... ){
  #
  #   if(is.vector(x)){
  #     theta=x
  #     theta[LB_index]=exp(x[LB_index])-lowerBound[LB_index]
  #     theta[UB_index]=-exp(x[UB_index])+upperBound[UB_index]
  #     theta[BB_index]=(upperBound[BB_index]-lowerBound[BB_index])*(exp(x[BB_index])/(exp(x[BB_index])+1))+lowerBound[BB_index]
  #
  #   }else{
  #     theta=x
  #
  #     theta[,LB_index]=exp(x[,LB_index])-repmat(lowerBound[LB_index],numrows = dim(x)[1] ,numcol = 1)
  #     theta[,UB_index]=-exp(x[,UB_index])+repmat(upperBound[UB_index],numrows = dim(x)[1] ,numcol = 1)
  #     theta[,BB_index]=(upperBoundMatrix-lowerBoundMatrix)[seq(1,dim(x)[1]),]*(exp(x[,BB_index])/(exp(x[,BB_index])+1))+lowerBoundMatrix[seq(1,dim(x)[1]),]
  #   }
  #
  #   nonlinearFunction(theta, ... )
  # }


  if(sum(abs(MO_weights))>0){
    nonlinearFunction_varTrans=function(x, ... ){

      if(is.vector(x)){
        theta=x
        theta[LB_index]=exp(x[LB_index])-lowerBound[LB_index]
        theta[UB_index]=-exp(x[UB_index])+upperBound[UB_index]
        theta[BB_index]=(upperBound[BB_index]-lowerBound[BB_index])*(exp(x[BB_index])/(exp(x[BB_index])+1))+lowerBound[BB_index]
        out=nonlinearFunction(theta, ... )

      }else{
        theta=x

        theta[,LB_index]=exp(x[,LB_index])-repmat(lowerBound[LB_index],numrows = dim(x)[1] ,numcol = 1)
        theta[,UB_index]=-exp(x[,UB_index])+repmat(upperBound[UB_index],numrows = dim(x)[1] ,numcol = 1)
        theta[,BB_index]=(upperBoundMatrix-lowerBoundMatrix)[seq(1,dim(x)[1]),]*(exp(x[,BB_index])/(exp(x[,BB_index])+1))+lowerBoundMatrix[seq(1,dim(x)[1]),]

        out=nonlinearFunction(theta, ... )
      }

      if(is.vector(out)){
        out=c(out,x[MO_index]*MO_weights[MO_index])
      }else{
        for(MO_ind_nu in seq(1,length(MO_index))[MO_index]){
          out=cbind(out,matrix(x[,MO_ind_nu]*MO_weights[MO_ind_nu],ncol = 1))
        }
      }

      return(out)

    }

  }else{
    nonlinearFunction_varTrans=function(x, ... ){

      if(is.vector(x)){
        theta=x
        theta[LB_index]=exp(x[LB_index])-lowerBound[LB_index]
        theta[UB_index]=-exp(x[UB_index])+upperBound[UB_index]
        theta[BB_index]=(upperBound[BB_index]-lowerBound[BB_index])*(exp(x[BB_index])/(exp(x[BB_index])+1))+lowerBound[BB_index]

      }else{
        theta=x

        theta[,LB_index]=exp(x[,LB_index])-repmat(lowerBound[LB_index],numrows = dim(x)[1] ,numcol = 1)
        theta[,UB_index]=-exp(x[,UB_index])+repmat(upperBound[UB_index],numrows = dim(x)[1] ,numcol = 1)
        theta[,BB_index]=(upperBoundMatrix-lowerBoundMatrix)[seq(1,dim(x)[1]),]*(exp(x[,BB_index])/(exp(x[,BB_index])+1))+lowerBoundMatrix[seq(1,dim(x)[1]),]
      }

      nonlinearFunction(theta, ... )
    }
  }


  if(bootstrapType==1){
    CGNM_bootstrap_result=Cluster_Gauss_Newton_method_core(nonlinearFunction = nonlinearFunction_varTrans, targetVector = targetVector, initial_lowerRange = initial_lowerRange, initial_upperRange = initial_upperRange, stayIn_initialRange = stayIn_initRange, num_minimizersToFind = num_minimizersToFind, num_iteration = num_iteration, saveLog = saveLog, runName = runName, textMemo = textMemo,algorithmParameter_initialLambda = algorithmParameter_initialLambda, algorithmParameter_gamma = algorithmParameter_gamma, algorithmVersion = algorithmVersion, initialIterateMatrix=bootstrapInitialMatrix, targetMatrix = bootstrapTargetMatrix, ParameterNames=ParameterNames,ReparameterizationDef = ReparameterizationDef,ParameterTransFormationDef=ParameterTransFormationDef,keepInitialDistribution = FALSE, folderSuffix="bootstrap", ... )
    CGNM_result$runSetting$targetMatrix=bootstrapTargetMatrix

  }else if(bootstrapType==2){
    bootstrapWeightMatrix=matrix(0,nrow = num_bootstrapSample,ncol = length(targetVector))

    for(i in seq(1,num_bootstrapSample)){
      randSample=table(sample(seq(1,length(targetVector)),length(targetVector), replace = TRUE))

      bootstrapWeightMatrix[i,as.numeric(names(randSample))]=as.numeric(randSample)
    }

    CGNM_bootstrap_result=Cluster_Gauss_Newton_method_core(nonlinearFunction = nonlinearFunction_varTrans, targetVector = targetVector, initial_lowerRange = initial_lowerRange, initial_upperRange = initial_upperRange, stayIn_initialRange = stayIn_initRange, num_minimizersToFind = num_minimizersToFind, num_iteration = num_iteration, saveLog = saveLog, runName = runName, textMemo = textMemo,algorithmParameter_initialLambda = algorithmParameter_initialLambda, algorithmParameter_gamma = algorithmParameter_gamma, algorithmVersion = algorithmVersion, initialIterateMatrix=bootstrapInitialMatrix, weightMatrix = bootstrapWeightMatrix, ParameterNames=ParameterNames,ReparameterizationDef = ReparameterizationDef,ParameterTransFormationDef=ParameterTransFormationDef,keepInitialDistribution = FALSE, folderSuffix="bootstrap", ... )
    CGNM_result$runSetting$targetMatrix=NA
  }else if(bootstrapType==3){
    bootstrapWeightMatrix=matrix(runif(num_bootstrapSample*length(targetVector)),nrow = num_bootstrapSample,ncol = length(targetVector))

    CGNM_bootstrap_result=Cluster_Gauss_Newton_method_core(nonlinearFunction = nonlinearFunction_varTrans, targetVector = targetVector, initial_lowerRange = initial_lowerRange, initial_upperRange = initial_upperRange, stayIn_initialRange = stayIn_initRange, num_minimizersToFind = num_minimizersToFind, num_iteration = num_iteration, saveLog = saveLog, runName = runName, textMemo = textMemo,algorithmParameter_initialLambda = algorithmParameter_initialLambda, algorithmParameter_gamma = algorithmParameter_gamma, algorithmVersion = algorithmVersion, initialIterateMatrix=bootstrapInitialMatrix, weightMatrix = bootstrapWeightMatrix, ParameterNames=ParameterNames,ReparameterizationDef = ReparameterizationDef,ParameterTransFormationDef=ParameterTransFormationDef,keepInitialDistribution = FALSE, folderSuffix="bootstrap", ... )
    CGNM_result$runSetting$targetMatrix=NA
  }

  CGNM_result$bootstrapX=CGNM_bootstrap_result$X
  CGNM_result$bootstrapY=CGNM_bootstrap_result$Y
  CGNM_result$runSetting$num_bootstrapSample=num_bootstrapSample
  CGNM_result$runSetting$indicesToUseAsInitialIterates=indicesToUseAsInitialIterates
  CGNM_result$runSetting$keepInitialDistribution=CGNM_bootstrap_result$keepInitialDistribution
  CGNM_result$runSetting$weightMatrix=CGNM_bootstrap_result$weightMatrix

  ReparameterizationDef=CGNM_result$runSetting$ReparameterizationDef
  ParameterTransFormationDef=CGNM_result$runSetting$ReparameterizationDef
  ParameterNames=CGNM_result$runSetting$ParameterNames

  out=data.frame(CGNM_result$bootstrapX)

  names(out)=paste0("x",seq(1,dim(out)[2]))

  bootstrapTheta=data.frame(row.names = seq(1,dim(out)[1]))

  for(i in seq(1,length(ParameterNames))){
    bootstrapTheta[,ParameterNames[i]]=with(out, eval(parse(text=ReparameterizationDef[i])))
  }

  CGNM_result$bootstrapTheta=bootstrapTheta
  CGNM_result$bootstrapParameterCombinations=bootstrapTheta



  if(CGNM_result$runSetting$saveLog){
    if(runName==""){
      save(file =paste0(runName,"CGNM_log_bootstrap/CGNM_bootstrapResult.RDATA"),CGNM_result)

    }else{
      save(file =paste0(runName,"_CGNM_log_bootstrap/CGNM_bootstrapResult.RDATA"),CGNM_result)

    }
  }

  return(CGNM_result)
}





#' @title Cluster_Gauss_Newton_EBE_method
#' @description obtain ebe empirical bayes estimate (EBE) using CGNM.
#' @param CGNM_result (required input) \emph{A list} stores the computational result from Cluster_Gauss_Newton_method() function in CGNM package.
#' @param nonlinearFunction (required input) \emph{A function with input of a vector x of real number of length n and output a vector y of real number of length m.} In the context of model fitting the nonlinearFunction is \strong{the model}.  Given the CGNM does not assume the uniqueness of the minimizer, m can be less than n.  Also CGNM does not assume any particular form of the nonlinear function and also does not require the function to be continuously differentiable (see Appendix D of our publication for an example when this function is discontinuous).
#' @param individualIndices_vec (required input) \emph{A string vector} where the each element corresponds to each observation and the value of the element is the ID that will be used for EBE. Must be the same as the length of the target vector.
#' @param numRepeat  (default: 1) \emph{A vector of integers} number of EBE per individual to obtain.
#' @param keepInitialDistribution (default: NA) \emph{A vector of TRUE or FALSE} of length n User can specify if the initial distribution of one of the input variable (e.g. parameter) to be kept as the initial iterate throughout CGNM iterations.
#' @param algorithmParameter_EBEweight (default: 9) \emph{A number}
#' @param ... Further arguments to be supplied to nonlinearFunction
#' @return list of a matrix X, Y,residual_history, initialX, bootstrapX, bootstrapY as well as a list runSetting.
#' \enumerate{\item X, Y, residual_history, initialX: identical to what was given as CGNM_result.
#' \item X: \emph{a num_bootstrapSample by n matrix} which stores the the X values that was sampled using residual resampling bootstrap analyses (In terms of model fitting this is the parameter combinations with variabilities that represent \strong{parameter estimation uncertainties}.).
#' \item Y: \emph{a num_bootstrapSample by m matrix} which stores the nonlinearFunction evaluated at the corresponding bootstrap analyses results in matrix bootstrapX above. In the context of model fitting each row corresponds to \strong{the model simulations}.
#' \item runSetting: identical to what is given as CGNM_result but in addition including num_bootstrapSample and indicesToUseAsInitialIterates.}
#' @examples
#' ##lip-flop kinetics (an example known to have two distinct solutions)
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
#' targetVector = observation, num_iteration = 10, num_minimizersToFind = 100,
#' initial_lowerRange = c(0.1,0.1,0.1), initial_upperRange =  c(10,10,10),
#' lowerBound=rep(0,3), ParameterNames=c("Ka","V1","CL_2"), saveLog = FALSE)
#'
#' CGNM_EBE=Cluster_Gauss_Newton_EBE_method(CGNM_result,
#'      nonlinearFunction=model_analytic_function, individualIndices_vec=seq(1,9))
#'
#'
#' @export
#' @import stats MASS

Cluster_Gauss_Newton_EBE_method <- function(CGNM_result, nonlinearFunction, individualIndices_vec, numRepeat=1, keepInitialDistribution=NA, algorithmParameter_EBEweight=1, ...){

  ## check if targetVector and individualIndices_vec has the same length
  if(length(CGNM_result$runSetting$targetVector)!=length(individualIndices_vec)){

  }else if(length(CGNM_result$runSetting$targetVector)-length(individualIndices_vec)==sum(CGNM_result$runSetting$MO_weights!=0)){

  }else{
    stop("Length of individualIndicies_vec needs to be the same as the targetVector")
  }


  unique_individualIndices_vec=unique(individualIndices_vec)
  numIndividual=length(unique_individualIndices_vec)

  message(paste0("Number of individuals: ", numIndividual))


  if(is.null(CGNM_result$bootstrapX)){
    acceptedMinimizers=as.data.frame(CGNM_result$X[acceptedIndices(CGNM_result),])
  }else{
    acceptedMinimizers=as.data.frame(CGNM_result$bootstrapX)
  }


  #  acceptedMinimizers=acceptedApproximateMinimizers(CGNM_result, algorithm = 2)

  if(dim(acceptedMinimizers)[1]<numRepeat){
    stop("Number of accepted minimizers is less than numRepeat.")
  }

  numBootstrap=numIndividual*numRepeat

  acceptedMinimizers_df=data.frame()

  for(i in seq(1, ceiling(numBootstrap/dim(acceptedMinimizers)[1]))){
    acceptedMinimizers_df=rbind(acceptedMinimizers_df, acceptedMinimizers)
  }

  acceptedMinimizers_matrix=as.matrix(acceptedMinimizers_df[seq(1,numBootstrap),])

  numMO=sum(CGNM_result$runSetting$MO_weights!=0)
  numTarget=length(CGNM_result$runSetting$targetVector)

  weightMatrix=matrix(1,numIndividual*numRepeat,numTarget)


  for(j in seq(1,numRepeat)){
    for(i in seq(1,numIndividual)){
      ID_nu=unique_individualIndices_vec[i]
      weightMatrix[i+numIndividual*(j-1),]=weightMatrix[i+numIndividual*(j-1),]+(individualIndices_vec==ID_nu)*algorithmParameter_EBEweight*length(individualIndices_vec)/sum(individualIndices_vec==ID_nu)
      if(numMO>0){
        weightMatrix[i+numIndividual*(j-1),seq(1+numTarget-numMO,numTarget)]=algorithmParameter_EBEweight

      }
    }
  }




  cutoff_pvalue=0.05
  numParametersIncluded=NA
  useAcceptedApproximateMinimizers=TRUE

  targetVector=CGNM_result$runSetting$targetVector
  initial_lowerRange=CGNM_result$runSetting$initial_lowerRange
  initial_upperRange=CGNM_result$runSetting$initial_upperRange
  num_iteration=CGNM_result$runSetting$num_iteration
  saveLog=CGNM_result$runSetting$saveLog
  runName=CGNM_result$runSetting$runName
  textMemo=CGNM_result$runSetting$textMemo
  algorithmParameter_initialLambda=CGNM_result$runSetting$algorithmParameter_initialLambda
  algorithmParameter_gamma=CGNM_result$runSetting$algorithmParameter_gamma
  algorithmVersion=CGNM_result$runSetting$algorithmVersion
  ReparameterizationDef=CGNM_result$runSetting$ReparameterizationDef
  ParameterNames=CGNM_result$runSetting$ParameterNames
  num_minimizersToFind=dim(weightMatrix)[1]
  ParameterTransFormationDef=CGNM_result$runSetting$ParameterTransFormationDef
  MO_weights=CGNM_result$runSetting$MO_weights
  MO_values=CGNM_result$runSetting$MO_values

  if(is.null(CGNM_result$runSetting$MO_weights)){
    MO_weights=rep(0, length(initial_lowerRange))
    MO_values=(initial_lowerRange+initial_upperRange)/2
  }else{
    MO_values=CGNM_result$runSetting$MO_values
    MO_weights=CGNM_result$runSetting$MO_weights
  }
  if(length(MO_weights)==1&is.na(MO_weights)[1]){
    MO_weights=rep(0, length(initial_lowerRange))
    MO_values=(initial_lowerRange+initial_upperRange)/2
  }

  MO_index=(MO_weights!=0)
  MO_values[!MO_index]=(initial_upperRange[!MO_index]+initial_lowerRange[!MO_index])/2


  stayIn_initRange=FALSE

  if(!is.null(CGNM_result$runSetting$stayIn_initialRange)){
    stayIn_initRange=CGNM_result$runSetting$stayIn_initialRange
  }

  lowerBound=NA
  if(!is.null(CGNM_result$runSetting$lowerBound)){
    lowerBound=CGNM_result$runSetting$lowerBound
  }


  upperBound=NA
  if(!is.null(CGNM_result$runSetting$upperBound)){
    upperBound=CGNM_result$runSetting$upperBound
  }

  ParameterNames=NA
  if(!is.null(CGNM_result$runSetting$ParameterNames)){
    ParameterNames=CGNM_result$runSetting$ParameterNames
  }

  LB_index=(!is.na(as.numeric(lowerBound)))&is.na(upperBound)
  UB_index=(!is.na(as.numeric(upperBound)))&is.na(lowerBound)
  BB_index=(!is.na(as.numeric(lowerBound)))&(!is.na(as.numeric(upperBound)))


  X_LR=initial_lowerRange
  X_LR[LB_index]=log(initial_lowerRange[LB_index]-lowerBound[LB_index])
  X_LR[UB_index]=log(-initial_lowerRange[UB_index]+upperBound[UB_index])
  X_LR[BB_index]=log((initial_lowerRange[BB_index]-lowerBound[BB_index])/(upperBound[BB_index]-initial_lowerRange[BB_index]))

  ParameterTransFormationDef=rep("", length(ParameterNames))
  ParameterTransFormationDef[LB_index]=paste0("log(",ParameterNames[LB_index],"-",lowerBound[LB_index],")")
  ParameterTransFormationDef[UB_index]=paste0("log(-",ParameterNames[UB_index],"+",upperBound[UB_index],")")
  ParameterTransFormationDef[BB_index]=paste0("log((",ParameterNames[BB_index],"-",lowerBound[BB_index],")/(",upperBound[BB_index],"-",ParameterNames[BB_index],"))")


  X_UR=initial_upperRange
  X_UR[LB_index]=log(initial_upperRange[LB_index]-lowerBound[LB_index])
  X_UR[UB_index]=log(-initial_upperRange[UB_index]+upperBound[UB_index])
  X_UR[BB_index]=log((initial_upperRange[BB_index]-lowerBound[BB_index])/(upperBound[BB_index]-initial_upperRange[BB_index]))


  MO_val_converted=MO_values
  MO_val_converted[LB_index]=log(MO_val_converted[LB_index]-lowerBound[LB_index])
  MO_val_converted[UB_index]=log(-MO_val_converted[UB_index]+upperBound[UB_index])
  MO_val_converted[BB_index]=log((MO_val_converted[BB_index]-lowerBound[BB_index])/(upperBound[BB_index]-MO_val_converted[BB_index]))


  lowerBoundMatrix=repmat(lowerBound[BB_index],numrows = num_minimizersToFind ,numcol = 1)
  upperBoundMatrix=repmat(upperBound[BB_index],numrows = num_minimizersToFind ,numcol = 1)

  if(sum(abs(MO_weights))>0){
    nonlinearFunction_varTrans=function(x, ... ){

      if(is.vector(x)){
        theta=x
        theta[LB_index]=exp(x[LB_index])-lowerBound[LB_index]
        theta[UB_index]=-exp(x[UB_index])+upperBound[UB_index]
        theta[BB_index]=(upperBound[BB_index]-lowerBound[BB_index])*(exp(x[BB_index])/(exp(x[BB_index])+1))+lowerBound[BB_index]
        out=nonlinearFunction(theta, ... )

      }else{
        theta=x

        theta[,LB_index]=exp(x[,LB_index])-repmat(lowerBound[LB_index],numrows = dim(x)[1] ,numcol = 1)
        theta[,UB_index]=-exp(x[,UB_index])+repmat(upperBound[UB_index],numrows = dim(x)[1] ,numcol = 1)
        theta[,BB_index]=(upperBoundMatrix-lowerBoundMatrix)[seq(1,dim(x)[1]),]*(exp(x[,BB_index])/(exp(x[,BB_index])+1))+lowerBoundMatrix[seq(1,dim(x)[1]),]

        out=nonlinearFunction(theta, ... )
      }

      if(is.vector(out)){
        out=c(out,x[MO_index]*MO_weights[MO_index])
      }else{
        for(MO_ind_nu in seq(1,length(MO_index))[MO_index]){
          out=cbind(out,matrix(x[,MO_ind_nu]*MO_weights[MO_ind_nu],ncol = 1))
        }
      }

      return(out)

    }

#    targetVector=c(targetVector, MO_val_converted[MO_index]*MO_weights[MO_index])

  }else{
    nonlinearFunction_varTrans=function(x, ... ){

      if(is.vector(x)){
        theta=x
        theta[LB_index]=exp(x[LB_index])-lowerBound[LB_index]
        theta[UB_index]=-exp(x[UB_index])+upperBound[UB_index]
        theta[BB_index]=(upperBound[BB_index]-lowerBound[BB_index])*(exp(x[BB_index])/(exp(x[BB_index])+1))+lowerBound[BB_index]

      }else{
        theta=x

        theta[,LB_index]=exp(x[,LB_index])-repmat(lowerBound[LB_index],numrows = dim(x)[1] ,numcol = 1)
        theta[,UB_index]=-exp(x[,UB_index])+repmat(upperBound[UB_index],numrows = dim(x)[1] ,numcol = 1)
        theta[,BB_index]=(upperBoundMatrix-lowerBoundMatrix)[seq(1,dim(x)[1]),]*(exp(x[,BB_index])/(exp(x[,BB_index])+1))+lowerBoundMatrix[seq(1,dim(x)[1]),]
      }

      nonlinearFunction(theta, ... )
    }
  }
  #
  #
  # nonlinearFunction_varTrans=function(x, ... ){
  #
  #   if(is.vector(x)){
  #     theta=x
  #     theta[LB_index]=exp(x[LB_index])-lowerBound[LB_index]
  #     theta[UB_index]=-exp(x[UB_index])+upperBound[UB_index]
  #     theta[BB_index]=(upperBound[BB_index]-lowerBound[BB_index])*(exp(x[BB_index])/(exp(x[BB_index])+1))+lowerBound[BB_index]
  #
  #   }else{
  #     theta=x
  #
  #     theta[,LB_index]=exp(x[,LB_index])-repmat(lowerBound[LB_index],numrows = dim(x)[1] ,numcol = 1)
  #     theta[,UB_index]=-exp(x[,UB_index])+repmat(upperBound[UB_index],numrows = dim(x)[1] ,numcol = 1)
  #     theta[,BB_index]=(upperBoundMatrix-lowerBoundMatrix)[seq(1,dim(x)[1]),]*(exp(x[,BB_index])/(exp(x[,BB_index])+1))+lowerBoundMatrix[seq(1,dim(x)[1]),]
  #   }
  #
  #   nonlinearFunction(theta, ... )
  # }

  CGNM_EBE_result=Cluster_Gauss_Newton_method_core(nonlinearFunction = nonlinearFunction_varTrans, targetVector = targetVector, initial_lowerRange = initial_lowerRange, initial_upperRange = initial_upperRange, stayIn_initialRange = stayIn_initRange, num_minimizersToFind = num_minimizersToFind, num_iteration = num_iteration, saveLog = saveLog, runName = runName, textMemo = textMemo,algorithmParameter_initialLambda = algorithmParameter_initialLambda, algorithmParameter_gamma = algorithmParameter_gamma, algorithmVersion = algorithmVersion, initialIterateMatrix=acceptedMinimizers_matrix, weightMatrix = weightMatrix, ParameterNames=ParameterNames,ReparameterizationDef = ReparameterizationDef,ParameterTransFormationDef=ParameterTransFormationDef,keepInitialDistribution = keepInitialDistribution, folderSuffix="EBE", ... )

  CGNM_result=CGNM_EBE_result

  CGNM_result$runSetting$targetMatrix=NA


  CGNM_result$EBE_X=CGNM_EBE_result$X
  CGNM_result$EBE_Y=CGNM_EBE_result$Y
  CGNM_result$runSetting$keepInitialDistribution=CGNM_EBE_result$runSetting$keepInitialDistribution
  CGNM_result$runSetting$weightMatrix=CGNM_EBE_result$runSetting$weightMatrix

  ReparameterizationDef=CGNM_result$runSetting$ReparameterizationDef
  ParameterNames=CGNM_result$runSetting$ParameterNames

  out=data.frame(CGNM_result$EBE_X)

  names(out)=paste0("x",seq(1,dim(out)[2]))

  EBETheta=data.frame(row.names = seq(1,dim(out)[1]))

  for(i in seq(1,length(ParameterNames))){
    EBETheta[,ParameterNames[i]]=with(out, eval(parse(text=ReparameterizationDef[i])))
  }

  CGNM_result$EBE_ParameterCombinations=EBETheta

  CGNM_result$runSetting$individualIndices_vec=individualIndices_vec
  CGNM_result$runSetting$numRepeat=numRepeat
  CGNM_result$runSetting$keepInitialDistribution=keepInitialDistribution
  CGNM_result$runSetting$algorithmParameter_EBEweight=algorithmParameter_EBEweight

  if(numRepeat>1){
    CGNM_result$EBE_ParameterCombinations_list=list()
    for(individual_nu in unique(individualIndices_vec)){
      CGNM_result$EBE_ParameterCombinations_list[[individual_nu]]=CGNM_result$EBE_ParameterCombinations[which(unique_individualIndices_vec==individual_nu)+seq(0,numRepeat-1)*length(unique_individualIndices_vec),]
    }
    CGNM_result$EBE_ParameterCombinations_list[[individual_nu]]=CGNM_result$EBE_ParameterCombinations_list[[individual_nu]][order(CGNM_result$residual_history[which(unique_individualIndices_vec==individual_nu)+seq(0,numRepeat-1)*length(unique_individualIndices_vec),dim(CGNM_result$residual_history)[2]]),]
  }

  if(CGNM_result$runSetting$saveLog){
    if(runName==""){
      save(file =paste0(runName,"CGNM_log_EBE/CGNM_EBE_Result.RDATA"),CGNM_result)
    }else{
      save(file =paste0(runName,"_CGNM_log_EBE/CGNM_EBE_Result.RDATA"),CGNM_result)
    }
  }

  return(CGNM_result)
}




#' @title shinyCGNM
#' @description Start Shiny app that assist the user to make R-script to conduct PBPK model fitting using CGNM.
#' @return NULL and start graphifcal user interface.
#' @examples
#'\dontrun{
#' shinyCGNM()
#'}
#' @export
#' @import shiny
shinyCGNM=function(){
  list.of.packages <- c("shiny", "shinyjs","DT","dplyr","shinydashboard", "rxode2", "mlr3misc", "shinybusy", "spsComps")
  new.packages <- list.of.packages[!(list.of.packages %in% .packages(all.available = TRUE))]
  if(length(new.packages)>0) {
    stop("The following packages are missing. Please install them:\n", paste(new.packages, collapse = "\n"),"\n\nYou can use the following command to install them:\ninstall.packages(c(\"",paste(new.packages, collapse = "\",\""),"\"))")
  }

  appDir <- system.file("shinyCGNM", package = "CGNM")
  if (appDir == "") {
    stop("Could not find shiny app directory. Try re-installing `CGNM`.", call. = FALSE)
  }

  shiny::runApp(appDir, display.mode = "normal")

}
