#' @title PSO Based Weighted Ensemble Algorithm for Volatility Modelling
#' @param Data Univariate Time Series Data
#' @param SplitR Split Ratio
#' @import stats xts rumidas rugarch WeightedEnsemble Metrics zoo
#' @return
#' \itemize{
#'   \item TrainFitted: Fitted values for the train series
#'   \item TestPred: Prediction for the test series
#'   \item Accuracy: Accuracy metric of the proposed model
#'   \item Weights: Weights of the ensemble
#' }
#' @export
#'
#' @examples
#' \donttest{
#' library("PWEV")
#' data<- as.ts(rnorm(200,100,50))
#' Result <- PWEV(data, 0.9)
#' }
#' @references
#' \itemize{
#'\item Paul, R.K., Das, T. and Yeasin, M., 2023. Ensemble of time series and machine learning model for forecasting volatility in agricultural prices. National Academy Science Letters, 46(3), pp.185-188.
#'\item Yeasin, M. and Paul, R.K., 2024. OptiSembleForecasting: optimization-based ensemble forecasting using MCS algorithm and PCA-based error index. The Journal of Supercomputing, 80(2), pp.1568-1597.
#' }

PWEV<-function(Data, SplitR){
  Accuracy_metric<-NULL
  data<-as.ts(Data)
  n <- length(data)

  data_train<-data[1:(n*SplitR)]
  data_test<-data[-c(1:(n*SplitR))]

  specs = ugarchspec(variance.model=list(model="sGARCH", garchOrder=c(1,1)))
  model_sgarch<-ugarchfit(spec=specs,data_train)
  train_sgarch<-model_sgarch@fit$fitted.values
  forecast_sgarch<-ugarchforecast(model_sgarch, n.ahead=length(data_test))
  test_sgarch<-forecast_sgarch@forecast$seriesFor


  specg = ugarchspec(variance.model=list(model="gjrGARCH", garchOrder=c(1,1)))
  model_gjrgarch<-ugarchfit(spec=specg,data_train)
  train_gjgarch<-model_gjrgarch@fit$fitted.values
  forecast_gjrgarch<-ugarchforecast(model_gjrgarch, n.ahead=length(data_test))
  test_gjgarch<-forecast_gjrgarch@forecast$seriesFor


  spece = ugarchspec(variance.model=list(model="iGARCH", garchOrder=c(1,1)))
  model_ewma<-ugarchfit(spec=spece,data_train)
  train_ewma<-model_ewma@fit$fitted.values
  forecast_ewma<-ugarchforecast(model_ewma, n.ahead=length(data_test))
  test_ewma<-forecast_ewma@forecast$seriesFor
  test_ewma1<-forecast_ewma@forecast$sigmaFor


  rv<-(data[-1]/data[-length(data)])^2
  RV_data<-as.xts(coredata(data),seq(as.Date("2010-01-01"),
                                     by = "month", length.out = length(data)))

  fit_mmidas<-umemfit("MEM","NO",x=RV_data,out_of_sample=length(data_test))
  train_mem<-as.vector(fit_mmidas$est_in_s)
  test_mem<-as.vector(fit_mmidas$est_oos)

  Train_Fitted<-cbind(data_train, train_sgarch, train_gjgarch,train_ewma, train_mem)
  Test_Predicted<-cbind(data_test,test_sgarch, test_gjgarch, test_ewma, test_mem)

  colnames(Train_Fitted)<-c("Actual", "SGARCH", "GJGARCH","EWMA", "MEM")
  colnames(Test_Predicted)<-c("Actual", "SGARCH", "GJGARCH","EWMA", "MEM")

  train<-Train_Fitted
  test<-Test_Predicted

  OptiSemble<-WeightedEnsemble(df=train)
  weight<-as.vector(OptiSemble$Weights)
  trainEnsemble<-as.vector(OptiSemble$Optimized_Result)

  testEnsemble<-rowSums(cbind(weight[1]*test[,2],weight[2]*test[,3],weight[3]*test[,4],weight[4]*test[,5]))

  Accuracy_metric<-function(Train_data, Test_data){
    compare_data<-Train_data
    AccuracyTable<-matrix(nrow = (ncol(compare_data)-1),ncol = 9)
    for (c in 2:ncol(compare_data)) {
      test_rmse<-rmse(compare_data[,1],compare_data[,c])
      test_mape<-mape(compare_data[,1],compare_data[,c])
      test_mae<-mae(compare_data[,1],compare_data[,c])
      test_rrse<-rrse(compare_data[,1],compare_data[,c])
      test_mdae <-mdae(compare_data[,1],compare_data[,c])
      test_rmsle<-rmsle(compare_data[,1],compare_data[,c])
      test_rae<-rae(compare_data[,1],compare_data[,c])
      test_smape<-smape(compare_data[,1],compare_data[,c])

      test_r_square<-as.numeric(((cor(compare_data[,1],compare_data[,c]))^2))


      AccuracyTable[(c-1),1]<-round(test_rmse,digits = 4)
      AccuracyTable[(c-1),2]<-round(test_mape,digits = 4)
      AccuracyTable[(c-1),3]<-round(test_mae,digits = 4)
      AccuracyTable[(c-1),4]<-round(test_rrse,digits = 4)
      AccuracyTable[(c-1),5]<-round(test_mdae,digits = 4)
      AccuracyTable[(c-1),6]<-round(test_rmsle,digits = 4)
      AccuracyTable[(c-1),7]<-round(test_rae,digits = 4)
      AccuracyTable[(c-1),8]<-round(test_smape,digits = 4)
      AccuracyTable[(c-1),9]<-round(test_r_square,digits = 4)
    }
    compare_data<-as.data.frame(Test_data)
    AccuracyTable1<-matrix(nrow = (ncol(compare_data)-1),ncol = 9)

    for (c in 2:ncol(compare_data)) {
      test_rmse<-rmse(compare_data[,1],compare_data[,c])
      test_mape<-mape(compare_data[,1],compare_data[,c])
      test_mae<-mae(compare_data[,1],compare_data[,c])
      test_rrse<-rrse(compare_data[,1],compare_data[,c])
      test_mdae <-mdae(compare_data[,1],compare_data[,c])
      test_rmsle<-rmsle(compare_data[,1],compare_data[,c])
      test_rae<-rae(compare_data[,1],compare_data[,c])
      test_smape<-smape(compare_data[,1],compare_data[,c])

      test_r_square<-as.numeric(((cor(compare_data[,1],compare_data[,c]))^2))


      AccuracyTable1[(c-1),1]<-round(test_rmse,digits = 4)
      AccuracyTable1[(c-1),2]<-round(test_mape,digits = 4)
      AccuracyTable1[(c-1),3]<-round(test_mae,digits = 4)
      AccuracyTable1[(c-1),4]<-round(test_rrse,digits = 4)
      AccuracyTable1[(c-1),5]<-round(test_mdae,digits = 4)
      AccuracyTable1[(c-1),6]<-round(test_rmsle,digits = 4)
      AccuracyTable1[(c-1),7]<-round(test_rae,digits = 4)
      AccuracyTable1[(c-1),8]<-round(test_smape,digits = 4)
      AccuracyTable1[(c-1),9]<-round(test_r_square,digits = 4)
    }
    ACC<-cbind(AccuracyTable,AccuracyTable1)
    name1<-paste0("Train_",c("RMSE","MAPE","MAE","RRSE","MADE","RMSLE","RAE","SMAPE", "R2"))
    name2<-paste0("Test_",c("RMSE","MAPE","MAE","RRSE","MADE","RMSLE","RAE","SMAPE", "R2"))
    name<-c(name1, name2)
    colnames(ACC)<-name
    rownames(ACC)<-colnames(Train_data[,-1])
    return(Accuracy=ACC)
  }


  TrainFinal<-cbind(train, trainEnsemble)
  TestFinal<-cbind(test, testEnsemble)
  colnames(TrainFinal)<-c(colnames(train), "Ensemble")
  colnames(TestFinal)<-colnames(TrainFinal)
  ACC<-Accuracy_metric(TrainFinal, Test_data = TestFinal)

  Results<-list(TrainFitted=TrainFinal, TestPred=TestFinal, Accuracy=ACC, Weights=weight)
  return(Results)
}
