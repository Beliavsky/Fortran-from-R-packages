

Rcpp::sourceCpp("tests/hashtest/hashtest.cpp", verbose = T)
segSize = 10L; N = 5000000; r = 5
# while(T){
x = abs(matrix(as.integer(round(rnorm(N * segSize) * r)), nrow = segSize))
y = abs(matrix(as.integer(round(rnorm(N * segSize) * r)), nrow = segSize))
# y = x
tmp = test(x, y, verbose = F)
d = sum(tmp$isYin != tmp$isYinViaFlattened)
# cat(sum(tmp$isYin), " ")
# if(d != 0) break
# }











