#' @export
rstandard.Arima <-
function(model,...)
{
  model$residuals/model$sigma2^.5
}