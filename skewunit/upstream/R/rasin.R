rasin <-
function(n) 
{
if (is.null(n)) stop("sample size must be specified")
if (round(n) != n | n <= 0) stop("sample size must be a positive integer")
rbeta(n, shape1=1/2, shape2=1/2)
}
