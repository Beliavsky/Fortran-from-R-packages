

Rcpp::sourceCpp("src/arbitraryDimFLSSStests/arithmetic.cpp", verbose = T)


N = 1000
x = sample(2L, N, replace = T) - 1L
y = sample(2L, N, replace = T) - 1L
tmp = addHintTest(x, y)
sum(abs(tmp$correctRst - tmp$rst))




Rcpp::sourceCpp("src/arbitraryDimFLSSStests/arithmetic.cpp", verbose = T)

while(T){


N = sample(1000, 1)
x = sample(2L, N, replace = T) - 1L
y = sample(2L, N, replace = T) - 1L
if(paste0(rev(x), collapse = "") < paste0(rev(y), collapse = ""))
{
  tmp = x; x = y; y = tmp
}
tmp = subHintTest(x, y)
err = sum(abs(tmp$correctRst - tmp$rst))
if(err != 0) break


}










