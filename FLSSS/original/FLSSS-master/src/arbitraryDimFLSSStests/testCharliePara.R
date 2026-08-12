

# Rcpp::sourceCpp("src/arbitraryDimFLSSStests/testCharliePara.cpp", verbose = T)
N = 1e7
x = runif(N)
y = runif(N)


maxCore = 10
discount = 1
Nscaler = 3
Rcpp::sourceCpp("src/arbitraryDimFLSSStests/testCharliePara.cpp", verbose = T)
Rcpp::sourceCpp("src/arbitraryDimFLSSStests/testTbbPara.cpp", verbose = T)


a = x + 0.0
b = y + 0.0
system.time({
  tmp = testTbbPara(a, b, discount = discount, Nscaler = Nscaler, maxCore = maxCore)
}); tmp


a = x + 0.0
b = y + 0.0
system.time({
  tmp2 = testCharliePara(a, b, discount = discount, Nscaler = Nscaler, maxCore = maxCore)
})


tmp / tmp2 - 1




























