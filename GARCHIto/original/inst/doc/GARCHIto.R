## -----------------------------------------------------------------------------
require(GARCHIto) 
data("sample_data")  
model_unified=UnifiedEst(sample_data$BPV, sample_data$return)
model_unified$coefficients # estimated model parameters 

## -----------------------------------------------------------------------------
# without the consideration of price jumps
model_realized_NJ=RealizedEst(sample_data$BPV)
model_realized_NJ$coefficients 
# with the consideration of price jumps 
model_realized=RealizedEst(sample_data$BPV, sample_data$JV)
model_realized$coefficients 

## ---- fig.align='center', fig.height=4, fig.width=7---------------------------
plot(model_unified$sigma, cex=0.5, type="o", ylim=c(0,0.00035),
     main="estimated conditional volatilities", ylab="", xlab="")
lines(model_realized_NJ$sigma,cex=0.5,type="o",col="blue",lty=2)
lines(model_realized$sigma,cex=0.5,type="o",col="red",lty=3, lwd=0.5)
legend("topleft", cex=0.8,
       legend=c("Unified GARCH-Ito", "Realized GARCH-Ito No Jump","Realized GARCH-Ito with Jump"),
       col = c("black", "blue", "red"),
       lty=c(1,2,3))

## -----------------------------------------------------------------------------
c(model_unified$pred, model_realized_NJ$pred, model_realized$pred) 

## -----------------------------------------------------------------------------
# conduct out of sample volatility forecasting and compute the mean squared prediction error 
error_unified=NULL
error_realized_NJ=NULL
error_realized=NULL
for (i in 560:603){
  sink("file")
  model1=UnifiedEst(sample_data$BPV[1:i], sample_data$return[1:i])
  error_unified=c(error_unified, (model1$pred-sample_data$BPV[i+1])^2)
  model2=RealizedEst(sample_data$BPV[1:i])
  error_realized_NJ=c(error_realized_NJ, (model2$pred-sample_data$BPV[i+1])^2)
  model3=RealizedEst(sample_data$BPV[1:i], sample_data$JV[1:i])
  error_realized=c(error_realized, (model3$pred-sample_data$BPV[i+1])^2)
  sink()
}

error=c(mean(error_unified), mean(error_realized_NJ), mean(error_realized))
names(error)=c("Unified GARCH-Ito", "Realized GARCH-Ito No Jump", "Realized GARCH-Ito with Jump")
error

## ---- eval=FALSE--------------------------------------------------------------
#  # without the consideration of price jumps
#  RealizedEst_Option(RV, NV) # homogeneous error
#  RealizedEst_Option(RV, NV, homogeneous=FALSE ) # heterogeneous error
#  # with the consideration of price jumps
#  RealizedEst_Option(RV, JV, NV) # homogeneous error
#  RealizedEst_Option(RV, JV, NV, homogeneous=FALSE) # heterogeneous error

