#' @export
`fitted.Arima` <-
function (object,...) 
{
fitted=eval(object$call$x)-object$residuals
fitted
}

