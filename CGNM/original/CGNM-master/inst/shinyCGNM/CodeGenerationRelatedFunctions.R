
repeatDetect=function(vector_in){
  vector_nu=vector_in

  first_index=1

  firstIndex_vec=c()
  lastIndex_vec=c()
  value_vec=c()

  while(length(vector_nu)>0){

    indecies=which(vector_nu==vector_nu[1])

    last_index=which(indecies-seq(1,length(indecies))!=0)[1]-1

    if(is.na(last_index)){
      last_index=length(indecies)
    }

    firstIndex_vec=c(firstIndex_vec,first_index)
    lastIndex_vec=c(lastIndex_vec,last_index+first_index-1)
    value_vec=c(value_vec, vector_nu[1])

    if(is.na(last_index)|last_index==length(vector_nu)){
      vector_nu=c()
    }else{
      vector_nu=vector_nu[seq(last_index+1,length(vector_nu))]
      first_index=last_index+first_index
    }
  }

  return(data.frame(firstIndex=firstIndex_vec, lastIndex=lastIndex_vec, value=value_vec))

}

condenseVec_toText=function(vector_in){
  repeat_df=repeatDetect(vector_in)
  out_text="c("

  for(i in 1:nrow(repeat_df)){
    if(repeat_df$firstIndex[i]==repeat_df$lastIndex[i]){
      out_text=paste0(out_text, repeat_df$value[i],",")
    }else{
      out_text=paste0(out_text, "rep(", repeat_df$value[i], ",", repeat_df$lastIndex[i]-repeat_df$firstIndex[i]+1, "),")
    }
  }

  out_text=paste0(substr(out_text,1,nchar(out_text)-1),")")
  return(out_text)
}


condenseVec_toVec=function(vector_in){
  repeat_df=repeatDetect(vector_in)
  out_vec=c()

  for(i in 1:nrow(repeat_df)){
    if(repeat_df$firstIndex[i]==repeat_df$lastIndex[i]){
      out_vec=c(out_vec, repeat_df$value[i])
    }else{
      out_vec=c(out_vec,paste0("rep(", repeat_df$value[i], ",", repeat_df$lastIndex[i]-repeat_df$firstIndex[i]+1, ")"))
    }
  }

  return(out_vec)
}





pasteWithCollapse_WithApprox80charLimits=function(paste_vec, collapse){
 # paste_vec=condenseVec_toVec(paste_vec)
  numTotalChar=nchar(paste(paste_vec, collapse = collapse))
  numToDivide=ceiling(numTotalChar/80)
  numPerDivide=ceiling(length(paste_vec)/numToDivide)

  splitted_list=split(paste_vec, ceiling(seq_along(paste_vec)/numPerDivide))
  outText_vec=c()
  for(splitNameNu in names(splitted_list)){
    outText_vec=c(outText_vec,  paste(splitted_list[[splitNameNu]], collapse = collapse))
  }

  return(paste(outText_vec, collapse = gsub(",",",\n",collapse)))
}


puttogether_model_code = function(input, ll, rv, for_simulation=FALSE) {
  ParameterName_vec = ll$parameterInfo_dat$ParameterName
  Initial_upper_range_vec = ll$parameterInfo_dat$Initial_upper_range
  Initial_lower_range_vec = ll$parameterInfo_dat$Initial_lower_range
  Upper_bound_vec = ll$parameterInfo_dat$Upper_bound
  Lower_bound_vec = ll$parameterInfo_dat$Lower_bound
  MO_weight_vec = ll$parameterInfo_dat$MO_weight
  MO_value_vec = ll$parameterInfo_dat$MO_value
  runName=input$runNameText
  runName=gsub("[[:punct:]]","", runName)
  runName=gsub(" ","_", runName)

  unique_IDs = sort(unique(c(ll$ObservedData_dat$ID, ll$Dose_dat$ID)))

  UseIndVar = FALSE
  IndividualParameter_text = ""
  if (!is.null(ll$parameterInfo_dat$VaryByID)) {
    if (sum(ll$parameterInfo_dat$VaryByID != 0) > 0) {
      UseIndVar = TRUE

      ParameterName_vec = subset(ll$parameterInfo_dat, VaryByID == 0)$ParameterName
      Initial_upper_range_vec = subset(ll$parameterInfo_dat, VaryByID ==
                                         0)$Initial_upper_range
      Initial_lower_range_vec = subset(ll$parameterInfo_dat, VaryByID ==
                                         0)$Initial_lower_range
      Upper_bound_vec = subset(ll$parameterInfo_dat, VaryByID == 0)$Upper_bound
      Lower_bound_vec = subset(ll$parameterInfo_dat, VaryByID == 0)$Lower_bound

      MO_weight_vec = subset(ll$parameterInfo_dat, VaryByID == 0)$MO_weight
      MO_value_vec = subset(ll$parameterInfo_dat, VaryByID == 0)$MO_value


      temp_paraInfo = subset(ll$parameterInfo_dat, VaryByID != 0)
      indParameterNames = temp_paraInfo$ParameterName


      Parameter_text = paste(paste0(
        subset(ll$parameterInfo_dat, VaryByID == 0)$ParameterName,
        "=x[",
        seq(1, dim(
          subset(ll$parameterInfo_dat, VaryByID == 0)
        )[1]),
        "]"
      ),
      collapse = "\n")

      indParaNames = c()
      for (i in seq(1, dim(temp_paraInfo)[1])) {
        indParaNames = c(indParaNames,
                         paste0(temp_paraInfo$ParameterName[i], "_ID", unique_IDs))
        ParameterName_vec = c(
          ParameterName_vec,
          paste0(temp_paraInfo$ParameterName[i], "_ID", unique_IDs)
        )
        Initial_upper_range_vec = c(
          Initial_upper_range_vec,
          rep(
            temp_paraInfo$Initial_upper_range[i],
            length(unique_IDs)
          )
        )
        Initial_lower_range_vec = c(
          Initial_lower_range_vec,
          rep(
            temp_paraInfo$Initial_lower_range[i],
            length(unique_IDs)
          )
        )
        Upper_bound_vec = c(Upper_bound_vec,
                            rep(temp_paraInfo$Upper_bound[i], length(unique_IDs)))
        Lower_bound_vec = c(Lower_bound_vec,
                            rep(temp_paraInfo$Lower_bound[i], length(unique_IDs)))
        MO_weight_vec = c(MO_weight_vec,
                          rep(temp_paraInfo$MO_weight[i], length(unique_IDs)))
        MO_value_vec = c(MO_value_vec, rep(temp_paraInfo$MO_value[i], length(unique_IDs)))
      }


      Parameter_text = paste0(Parameter_text, "\n", paste0(paste0(
        indParaNames,
        "=x[",
        seq(
          dim(subset(ll$parameterInfo_dat, VaryByID == 0))[1] + 1,
          dim(subset(ll$parameterInfo_dat, VaryByID == 0))[1] + length(unique_IDs) *
            dim(temp_paraInfo)[1]
        ),
        "]"
      ),
      collapse = "\n"))

    } else{
      Parameter_text = paste(paste0(
        ll$parameterInfo_dat$ParameterName,
        "=x[",
        seq(1, dim(ll$parameterInfo_dat)[1]),
        "]"
      ),
      collapse = "\n")
    }
  } else{
    Parameter_text = paste(paste0(
      ll$parameterInfo_dat$ParameterName,
      "=x[",
      seq(1, dim(ll$parameterInfo_dat)[1]),
      "]"
    ),
    collapse = "\n")
  }

  odeCodeText = paste0("ODE_text=\"\n",
                       input$ODE_text,
                       "\"\n compiledModel=RxODE(ODE_text)")


  parameterText = paste0(ifelse(for_simulation,paste0("\nsimulation_function_",runName,"=function(x){\n"),paste0("\nmodel_function_",runName,"=function(x){\n")), Parameter_text)

  if (!UseIndVar) {
    parameterText = paste0(
      "\n\n",
      parameterText,
      "\n\nmodelingParameter=c(",
      pasteWithCollapse_WithApprox80charLimits(paste0(
        rv$ODEparameter, "=", rv$ODEparameter
      ), collapse = ","),
      ")\n\n"
    )
  }

  dataSet_text = paste0(
    #"\n\ndataSet=data.frame(seq=seq(1,",
    "\n\n",ifelse(for_simulation,"simulationDataSkelton", "dataSet"),"<<-subset(data.frame(seq=seq(1,",

    dim(ll$ObservedData_dat)[1],
    "),
      time=c(",
    pasteWithCollapse_WithApprox80charLimits(ll$ObservedData_dat$time, collapse = ","),
    "),
      ID=c(\"",
    pasteWithCollapse_WithApprox80charLimits(ll$ObservedData_dat$ID, collapse = "\",\""),
    "\"),
    ",ifelse(for_simulation,"",
      paste0("Observed_value=c(",
             pasteWithCollapse_WithApprox80charLimits(ll$ObservedData_dat$Observed_value, collapse = ","),
    "),")),"
      Observation_expression=c(\"",
    pasteWithCollapse_WithApprox80charLimits(ll$ObservedData_dat$Observation_expression, collapse = "\",\""),
    "\"))",ifelse(!for_simulation,
",!is.na(Observed_value))",")\n\n
uniqueTime=unique(simulationDataSkelton$time)\n"
  ))

  dose_obs_Text = paste0("\n\nsimResult_df=data.frame()\n")

  for (ID_nu in sort(unique_IDs)) {
    dose_df = subset(ll$Dose_dat, ID == ID_nu)
    obs_df = subset(ll$ObservedData_dat, ID == ID_nu)

    if (UseIndVar) {
      dose_obs_Text = paste0(dose_obs_Text, "\n\n## ID: ", ID_nu, "     ", paste( unique(obs_df$Memo), collapse = ", "),
                             "\n\n")



      dose_obs_Text = paste0(
        dose_obs_Text,
        paste(
          paste0(indParameterNames, "=", indParameterNames, "_ID", ID_nu),
          collapse = "\n"
        ),
        "\n\nmodelingParameter=c(",
        pasteWithCollapse_WithApprox80charLimits(
          paste0(rv$ODEparameter, "=", rv$ODEparameter),
          collapse = ","
        ),
        ")\n\n ev <- eventTable()"
      )

    } else{
      dose_obs_Text = paste0(dose_obs_Text, "\n\n ## ID: ", ID_nu, "     ", paste( unique(obs_df$Memo), collapse = ", "),"
      ev <- eventTable()")

    }

    if (dim(dose_df)[1] > 0) {
      for (i in seq(1, dim(dose_df)[1])) {
        if (is.na(dose_df$rate[i])) {
          dose_obs_Text = paste0(
            dose_obs_Text,
            "
      ev$add.dosing(dose = ",
            dose_df$dose[i],
            ", start.time=",
            dose_df$start.time[i],
            ", dosing.to=\"",
            dose_df$dosing.to[i],
            "\", nbr.doses=",
            dose_df$nbr.doses[i],
            ", dosing.interval=",
            dose_df$dosing.interval[i],
            ")"
          )
        } else{
          dose_obs_Text = paste0(
            dose_obs_Text,
            "
      ev$add.dosing(dose = ",
            dose_df$dose[i],
            ", rate = ",
            dose_df$rate[i],
            ", start.time=",
            dose_df$start.time[i],
            ", dosing.to=\"",
            dose_df$dosing.to[i],
            "\", nbr.doses=",
            dose_df$nbr.doses[i],
            ", dosing.interval=",
            dose_df$dosing.interval[i],
            ")"


          )
        }
      }
    }

   # ifelse(for_simulation,"uniqueTime",)


    if (dim(obs_df)[1] > 0) {
      dose_obs_Text =paste0(
        dose_obs_Text,
        "\n      ev$add.sampling(",ifelse(for_simulation,"uniqueTime",paste0("c(",
        paste(sort(as.numeric(
          unique(obs_df$time)
        )), collapse = ", "),
        ")")),
        ")\n      odeSol=compiledModel$solve(modelingParameter, ev)"
      )

      uniqueObsVariables = unique(obs_df$Observation_expression)

      for (obsVariable_nu in uniqueObsVariables) {
        dose_obs_Text = paste0(
          dose_obs_Text,
          "\n      simResult_df=rbind(simResult_df,data.frame(value=with(data.frame(odeSol), ",
          obsVariable_nu ,
          "), Observation_expression=\"",
          obsVariable_nu,
          "\", time=odeSol[,\"time\"], ID=\"",
          ID_nu,
          "\"))"
        )
      }

      if(for_simulation){
        dose_obs_Text=paste0(dose_obs_Text, paste0("\n\n      for(i in seq(2,dim(odeSol)[2])){
        simResult_df=rbind(simResult_df,data.frame(value=odeSol[,i], Observation_expression = colnames(odeSol)[i], time=odeSol[,1], ID=\"",ID_nu,"\"))
      }")
        )
      }
    }
  }
  dose_obs_Text = paste0(
    dose_obs_Text,
    "

      mergedData=merge(",ifelse(for_simulation,"simulationDataSkelton", "dataSet"),",simResult_df, all.x = TRUE)
      mergedData=mergedData[order(mergedData$seq),]
      return(mergedData$value)
}"
  )


  CGNM_run_text = paste0(
    paste0(
      "
ParaNames=c(\"",
      pasteWithCollapse_WithApprox80charLimits(ParameterName_vec, collapse = "\",\""),
      "\")
"
    ),

    paste0("UR=c(",
           paste(Initial_upper_range_vec, collapse = ","),
           ")
"),
    paste0("LR=c(",
           paste(Initial_lower_range_vec, collapse = ","),
           ")
"),
    paste0("U_bound=c(",
           paste(
             as.numeric(Upper_bound_vec), collapse = ","
           ),
           ")
"),
    paste0("L_bound=c(",
           paste(
             as.numeric(Lower_bound_vec), collapse = ","
           ),
           ")
"),
    paste0("observation=dataSet$Observed_value
")
  )

  if (ll$UseMO) {
    CGNM_run_text = paste0(
      CGNM_run_text,
      "MO_weights_vec=c(",
      paste(MO_weight_vec, collapse = ",") ,
      ")
MO_values_vec=c(",
      paste(MO_value_vec, collapse = ",") ,
      ")"
    )

  }

  testCode_text <-
    paste0(
      "library(CGNM)\nlibrary(rxode2)\n\n",ifelse(for_simulation,"","#If you wish to rerun CGNM regardless of if there is a log file available set it to TRUE\nForceReEstimation=FALSE\n\n"),
      odeCodeText,
      dataSet_text,
      parameterText,
      dose_obs_Text,
      ifelse(for_simulation,"",CGNM_run_text),
      "\n\n"
    )

  return(testCode_text)


}




makeCGNM_runCode = function(parallel = "none",
                            ll,
                            numIter = 25,
                            numMinimizersTofind = 250,
                            bootstrap = TRUE,
                            runName="") {

  runName=gsub("[[:punct:]]","", runName)
  runName=gsub(" ","_", runName)

  useResidualFunction = sum(ll$ObservedData_dat$ResidualError_model != 0) >
    0

  if (useResidualFunction) {
    CGNM_runOptions = "targetVector = rep(0,length(observation)), "
  } else{
    CGNM_runOptions = "targetVector = observation, "
  }

  if (ll$UseMO) {
    CGNM_runOptions = paste0(
      CGNM_runOptions,
      "initial_lowerRange = LR, initial_upperRange = UR, lowerBound = L_bound, upperBound = U_bound,ParameterNames = ParaNames, MO_weights=MO_weights_vec, MO_values=MO_values_vec, num_minimizersToFind = ",
      numMinimizersTofind,
      ", num_iteration=",
      numIter
    )

  } else{
    CGNM_runOptions = paste0(
      CGNM_runOptions,
      "initial_lowerRange = LR, initial_upperRange = UR, lowerBound = L_bound, upperBound = U_bound,ParameterNames = ParaNames, num_minimizersToFind = ",
      numMinimizersTofind,
      ", num_iteration=",
      numIter
    )
  }

  CGNM_runOptions = paste0(
    CGNM_runOptions,
    ",runName=\"",runName,"\""
  )


  out = ""
  if (useResidualFunction) {
    if (sum(ll$ObservedData_dat$ResidualError_model != 1) == 0) {
      out = paste0(
        out,
        "
residual_model=function(y_sim){
  out=(y_sim-dataSet$Observed_value)/y_sim

  return(out)
}"
      )
    } else{
      out = paste0(
        out,
        "
residual_model=function(y_sim){
  out=(y_sim-dataSet$Observed_value)/y_sim
  out[c(",
        paste(
          which(ll$ObservedData_dat$ResidualError_model == 0),
          collapse = ","
        )

        ,
        ")]=y_sim[c(",
        paste(
          which(ll$ObservedData_dat$ResidualError_model == 0),
          collapse = ","
        )

        ,
        ")]

  return(out)
}"
      )
    }
  }

  if (ll$UseMO) {
    out = paste0(
      out,
      "\n\n## CGNM R package above or equal to version 0.8.1 is necessary to run middle out method using MO_weights, MO_values options as implemented below."
    )
  }
  if(bootstrap){
    out=paste0(out, "\n\n", "if(!file.exists(\"",runName,"_CGNM_log_bootstrap/CGNM_bootstrapResult.RDATA\")|ForceReEstimation){"
    )
  }else{
    out=paste0(out, "\n\n", "if(!file.exists(\"",runName,"_CGNM_log/CGNM_result.RDATA\")|ForceReEstimation){"
    )
  }

  if (parallel == "none") {
    if (useResidualFunction) {
      out = paste0(
        out,
        "

model_function_withResidualmodel=function(x){
  return(residual_model(model_function_",runName,"(x)))
}

CGNM_result=Cluster_Gauss_Newton_method(model_function_withResidualmodel, ", CGNM_runOptions, ")
", ifelse(
  bootstrap,
  "CGNM_result=Cluster_Gauss_Newton_Bootstrap_method(CGNM_result,model_function_withResidualmodel)",
  ""
))

    } else{
      out = paste0(
        out,
        "

CGNM_result=Cluster_Gauss_Newton_method(model_function_",runName,", ",
        CGNM_runOptions,
        ")
",
        ifelse(
          bootstrap,
          paste0("CGNM_result=Cluster_Gauss_Newton_Bootstrap_method(CGNM_result,model_function_",runName,")"),
          ""
        )
      )
    }

  } else if (parallel == "win") {
    out = paste0(
      out,
      "

library(foreach)
library(doParallel)

numCoretoUse=detectCores()-1
registerDoParallel(numCoretoUse)
cluster=makeCluster(numCoretoUse, type =\"PSOCK\")
registerDoParallel(cl=cluster)

obsLength=length(observation)

# Given CGNM searches through wide range of parameter combination, it can encounter
# parameter combinations that is not feasible to evaluate. This try catch function
# is implemented within CGNM for regular functions but for the matrix functions
# user needs to implement outside of CGNM

modelFunction_tryCatch=function(x){
 out=tryCatch({",
      ifelse(
        useResidualFunction,
        paste0("residual_model(model_function_",runName,"(x))"),
        paste0("model_function_",runName,"(x)")
      ),
      "},
              error=function(cond) {rep(NA, obsLength)}
 )
 return(out)
}

model_matrix_function=function(x){
  X=as.matrix(x)

  if(is.matrix(X)){
    Y_list=foreach(i=1:dim(X)[1], .export = c(\"model_function_",runName,"\",",
ifelse(useResidualFunction, "\"residual_model\",", ""),
" \"modelFunction_tryCatch\", \"dataSet\", \"obsLength\", \"compiledModel\"), .packages = c(\"rxode2\"))%dopar%{
      modelFunction_tryCatch(as.numeric(X[i,]))
    }

    Y=t(matrix(unlist(Y_list),ncol=length(Y_list)))

  }else{

   Y= modelFunction_tryCatch(X)
  }

  return(Y)

}

CGNM_result=Cluster_Gauss_Newton_method(model_matrix_function, ",
CGNM_runOptions,
")
"
,
ifelse(
  bootstrap,
  "CGNM_result=Cluster_Gauss_Newton_Bootstrap_method(CGNM_result,model_matrix_function)",
  ""
),
"

stopCluster(cl=cluster)")

  } else if (parallel == "mac") {
    out = paste0(
      out,
      "

library(parallel)

    obsLength=length(observation)

    ## Given CGNM searches through wide range of parameter combination, it can encounter
    ## parameter combinations that is not feasible to evaluate. This try catch function
    ## is implemented within CGNM for regular functions but for the matrix functions
    ## user needs to implement outside of CGNM

    modelFunction_tryCatch=function(x){
      out=tryCatch({",
      ifelse(
        useResidualFunction,
        paste0("residual_model(model_function_",runName,"(x))"),
        paste0("model_function_",runName,"(x)")
      ),
      "},
                   error=function(cond) {rep(NA, obsLength)}
      )
      return(out)
    }

    model_matrix_function=function(x){
      Y_list=mclapply(split(x, rep(seq(1:nrow(x)),ncol(x))), modelFunction_tryCatch,mc.cores = (parallel::detectCores()-1), mc.preschedule = FALSE)

      Y=t(matrix(unlist(Y_list),ncol=length(Y_list)))

      return(Y)
    }

CGNM_result=Cluster_Gauss_Newton_method(model_matrix_function, ",
      CGNM_runOptions,
      ")
"
      ,
      ifelse(
        bootstrap,
        "CGNM_result=Cluster_Gauss_Newton_Bootstrap_method(CGNM_result,model_matrix_function)",
        ""
      )
    )

  }

if(bootstrap){
  out=paste0(out, "\n\n}else{\nload( file=\"",runName,"_CGNM_log_bootstrap/CGNM_bootstrapResult.RDATA\")\n}\n")
}else{
  out=paste0(out, "\n\n}else{\nload( file=\"",runName,"_CGNM_log/iteration_final.RDATA\")\n}\n")
}


if (ll$UseMO) {
  out = paste0(
    out,
    "


MO_para_names=CGNM_result$runSetting$ParameterNames[CGNM_result$runSetting$MO_weights!=0]
MO_values=CGNM_result$runSetting$MO_values[CGNM_result$runSetting$MO_weights!=0]

plot_goodnessOfFit(CGNM_result, independentVariableVector = c(dataSet$time, rep(0,length(MO_para_names))) ,dependentVariableTypeVector = c(paste(\"ID:\",dataSet$ID, dataSet$Observation_expression), MO_para_names) )+ggplot2::geom_point(colour=\"blue\")+ggplot2::labs(caption = \"Note the middleout values shown here are after transformation\")
plot_profileLikelihood(CGNM_result)+scale_x_continuous(trans=\"log10\")+ggplot2::geom_vline(data=data.frame(value=MO_values, parameterName=MO_para_names), aes(xintercept=value), colour=\"darkgrey\")
"
  )
} else{
  out = paste0(
    out,
    "\n\nplot_goodnessOfFit(CGNM_result, independentVariableVector = dataSet$time ,dependentVariableTypeVector = paste(\"ID:\",dataSet$ID, dataSet$Observation_expression))
plot_profileLikelihood(CGNM_result)+scale_x_continuous(trans=\"log10\")
"
  )
}


if (sum(ll$ObservedData_dat$ResidualError_model != 0) > 0) {
  out = paste0(
    out,
    "
plot_simulationWithCI(model_function_",runName,",parameter_matrix = CGNM_result$bootstrapParameterCombinations, independentVariableVector = dataSet$time, dependentVariableTypeVector =  paste( dataSet$Observation_expression,\"ID:\",dataSet$ID),
                      observationVector = dataSet$Observed_value, observationIndpendentVariableVector = dataSet$time, observationDependentVariableTypeVector =  paste( dataSet$Observation_expression,\"ID:\",dataSet$ID))+scale_y_continuous(trans=\"log10\")
"
  )

}
return(out)
}

simulation_code_text=function(runName){
  runName=gsub("[[:punct:]]","", runName)
  runName=gsub(" ","_", runName)

return( paste0("
load(\"",runName,"_CGNM_log_bootstrap/CGNM_bootstrapResult.RDATA\")
plot_simulationWithCI(simulation_function",ifelse(runName=="","","_"),runName,",parameter_matrix = CGNM_result$bootstrapParameterCombinations, independentVariableVector = simulationDataSkelton$time, dependentVariableTypeVector =  paste(simulationDataSkelton$Observation_expression, \"ID:\",simulationDataSkelton$ID))+scale_y_continuous(trans=\"log10\")
"))
}

