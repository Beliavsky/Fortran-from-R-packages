


# RawVector ksumHash(NumericMatrix v, IntegerVector lb, IntegerVector ub,
#                    int maxCore, int upscale = 20, double memlimGB = 4)
set.seed(123)


Rcpp::sourceCpp("src/arbitraryDimFLSSStests/ksum.cpp", verbose = T)
d = 3
N = 70
# len = 5


# iter = 0L
# while (T) {
# iter = iter + 1L
# # if (iter %% 1000L == 0) cat(iter, "")
# cat(iter, "")


len = sample(3:6, 1)
lb = sort(sample(N, len) - 1L)
ub = sort(sample(N, len) - 1L)
tmp = list(pmin(lb, ub), pmax(lb, ub))
lb = tmp[[1]]; ub = tmp[[2]]
maxCore = 10
upscale = 30
v = matrix(runif(N * d, 0, 2^64*0.99), nrow = d)
v[, ncol(v)] = runif(d, 0, 2^63*0.99)


system.time({
  # sink("debug.txt")
  rst = ksumHash(v, lb, ub, maxCore, upscale, verbose = T)
  # sink()
  })


# if (diff(rst$ics) != 0) break
# }



reticulate::source_python("src/arbitraryDimFLSSStests/btsubsetsum.py")
NofSums(lb, ub)
tmp = btsubsetSum(as.numeric(v), lb, ub, uniqueOut = 0)























