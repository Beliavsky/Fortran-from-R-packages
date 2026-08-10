
#' @title make_ShinyCGNM_observationData
#' @description
#' A helper function to write out the csv file that can be read in as the observation file in shinyCGNM
#' @param ID (required input) \emph{string or vector}
#' @param time (required input) \emph{number or numeric vector}
#' @param Observation_expression (required input) \emph{string or string vector}
#' @param Observed_value (required input) \emph{number or numeric vector}
#' @param ResidualError_model (required input) \emph{0 or 1} 0: additive residual model, 1: relative residual model
#' @param Memo (default: NA)  \emph{NA, string, or string vector} If TRUE plot absolute values of the residual.
#' @param fileName (default: NA) \emph{NA or string}
#' @return \emph{data.frame} if fileName is NA \emph{Null} if fileName is not NA but instead write out the csv file
#' @examples
#'
#' make_ShinyCGNM_observationData(
#' ID=1,
#' time=c(1,2,3,6,12,24),
#' Observation_expression="C_central",
#' Observed_value=c(0.1, 0.3, 0.6, 0.1, 0.05, 0.01),
#' ResidualError_model=1
#' )
#'
#' @export
#' @import utils

make_ShinyCGNM_observationData=function(ID, time, Observation_expression, Observed_value, ResidualError_model, Memo=NA, fileName=NA){

  out_df=data.frame(ID=ID,
                    time=time,
                    Observation_expression=Observation_expression,
                    Observed_value=Observed_value,
                    ResidualError_model=ResidualError_model,
                    Memo=Memo)

  if (is.na(fileName)){
    return(out_df)
  }else{
    write.csv(file = fileName, out_df, row.names = FALSE)
  }
}


#' @title make_ShinyCGNM_doseData
#' @description
#' A helper function to write out the csv file that can be read in as the dose file in shinyCGNM
#' @param ID (required input) \emph{string or vector}
#' @param dose (required input) \emph{number or numeric vector}
#' @param dosing.to (required input) \emph{string or string vector}
#' @param start.time (required input) \emph{number or numeric vector}
#' @param rate (required input) \emph{NA, number or numeric vector} set ti BA uf bolus
#' @param fileName (default: NA) \emph{NA or string}
#' @return \emph{data.frame} if fileName is NA \emph{Null} if fileName is not NA but instead write out the csv file
#' @examples
#'
#' make_ShinyCGNM_doseData(
#' ID=seq(1,5),
#' dose=10,
#' dosing.to="A_admin",
#' start.time=0,
#' rate=NA
#' )
#'
#' @export
#' @import utils

make_ShinyCGNM_doseData=function(ID, dose, dosing.to, start.time, rate, fileName=NA){

  out_df=data.frame(ID=ID,
                    dose=dose,
                    dosing.to=dosing.to,
                    start.time=start.time,
                    rate=rate)

  if (is.na(fileName)){
    return(out_df)
  }else{
    write.csv(file = fileName, out_df, row.names = FALSE)
  }
}
