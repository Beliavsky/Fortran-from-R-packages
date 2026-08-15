rsbeta <-
function(n, delta=1)
{
if (is.null(n)) stop("sample size must be specified")
if (round(n) != n | n <= 0) stop("sample size must be a positive integer")
if(delta<=0) stop("delta should be positive")
rbeta(n, shape1=delta, shape2=delta)
}
