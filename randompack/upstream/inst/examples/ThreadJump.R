library(parallel)

worker <- function(i) {              # R workers are separate processes
  library(randompack)                # a little inefficient, but ok
  rng <- randompack_rng("ranlux++")  # ...so cannot use duplicate
  rng$seed(123)
  rng$jump(32*i)
  rng$normal(1000)
}

cluster <- makeCluster(3)
results <- parLapply(cluster, 1:3, worker)  # must copy results on exit
stopCluster(cluster)                        # wait till workers finish
x <- rbind(results[[1]], results[[2]], results[[3]])
print(x)
