rtriang <-
function(n)
{
if (is.null(n)) stop("sample size must be specified")
if (round(n) != n | n <= 0) stop("sample size must be a positive integer")
u<-runif(n)
ifelse(u<=1/2, sqrt(u/2), 1-sqrt(0.5*(1-u)))
}
