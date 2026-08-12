

cost = t(matrix(as.integer(c(8,18,22,5,11,11,22,11,17,22,11,20,13,13,7,22,15,22,24,8,8,24,18,8,24,14,11,15,24,8,10,15,19,25,6,13,10,25,19,24,13,12,5,18,10,24,8,5,22,22,21,22,13,16,21,5,25,13,12,9,24,6,22,24,11,21,11,14,12,10,20,6,13,8,19,12,19,18,10,21,5,9,11,9,22,8,12,13,9,25,19,24,22,6,19,14,25,16,13,5,11,8,7,8,25,20,24,20,11,6,10,10,6,22,10,10,13,21,5,19,19,19,5,11,22,24,18,11,6,13,24,24,22,6,22,5,14,6,16,11,6,8,18,10,24,10,9,10,6,15,7,13,20,8,7,9,24,9,21,9,11,19,10,5,23,20,5,21,6,9,9,5,12,10,16,15,19,18,20,18,16,21,11,12,22,16,21,25,7,14,16,10)), ncol = 8))
profit = t(matrix(as.integer(c(25,23,20,16,19,22,20,16,15,22,15,21,20,23,20,22,19,25,25,24,21,17,23,17,16,19,22,22,19,23,17,24,15,24,18,19,20,24,25,25,19,24,18,21,16,25,15,20,20,18,23,23,23,17,19,16,24,24,17,23,19,22,23,25,23,18,19,24,20,17,23,23,16,16,15,23,15,15,25,22,17,20,19,16,17,17,20,17,17,18,16,18,15,25,22,17,17,23,21,20,24,22,25,17,22,20,16,22,21,23,24,15,22,25,18,19,19,17,22,23,24,21,23,17,21,19,19,17,18,24,15,15,17,18,15,24,19,21,23,24,17,20,16,21,18,21,22,23,22,15,18,15,21,22,15,23,21,25,25,23,20,16,25,17,15,15,18,16,19,24,18,17,21,18,24,25,18,23,21,15,24,23,18,18,23,23,16,20,20,19,25,21)), ncol = 8))
budget = as.integer(c(36, 35, 38, 34, 32, 34, 31, 34))
sink("sink.txt")
sol = FLSSS::auxGAPbb(cost, profit, budget, maxCore = 7, tlimit = 3600, ub = "MT", greedyBranching = 1, optim = "max")
sink()











load("data/gap60instances.Rdata")
sol = list()
for(i in 1L : 16L)
{
  cost = gap60instances[[i]]$cost
  profit = gap60instances[[i]]$profit
  budget = gap60instances[[i]]$budget
  sink(paste0("sink", i, ".txt"))
  sol[[i]] = FLSSS::auxGAPbb(cost, profit, budget, maxCore = 7, tlimit = 3600, ub = "MT", greedyBranching = 1, optim = "max")
  sink()
}


sol2 = list()
for(i in 1L : length(gap60instances))
{
  cost = gap60instances[[i]]$cost
  profit = gap60instances[[i]]$profit
  budget = gap60instances[[i]]$budget
  sink(paste0("sink", i, ".txt"))
  sol2[[i]] = FLSSS::auxGAPbbDp(cost, profit, budget, maxCore = 7, tlimit = 3600, greedyBranching = 1, optim = "max")
  sink()
}








makeData <- function(v)
{
  K = v[1]
  v = v[-1]
  rst = list()
  for(i in 1L : K)
  {
    m = v[1]
    n = v[2]
    mn = m * n
    profit = t(matrix(v[3L : (mn + 2L)], nrow = n))
    cost = t(matrix(v[(mn + 2L + 1L) : (mn + 2L + mn)], nrow = n))
    budget = v[(mn + 2L + mn + 1L) : (mn + 2L + mn + 1L + m - 1L)]
    v = v[-(1L : (2L + mn + mn + m))]
    rst[[i]] = list(cost = cost, profit = profit, budget = budget)
  }
  rst
}
opath = paste0("http://people.brunel.ac.uk/~mastjjb/jeb/orlib/files/gap", 1L : 12L, ".txt")
tmp = lapply(opath, function(x) as.integer(scan(url(x))))
gap60instances = unlist(lapply(tmp, function(x) makeData(x)), recursive = FALSE)
# save(smallGapInstances, file = "tests/smallGapInstances.Rdata")


load("tests/smallGapInstances.Rdata")


sol = list()
for(i in 1L : length(smallGapInstances))
{
  cost = smallGapInstances[[i]]$cost
  profit = smallGapInstances[[i]]$profit
  budget = smallGapInstances[[i]]$budget
  sink(paste0("sink", i, ".txt"))
  sol[[i]] = FLSSS:::auxGAPbb(cost, profit, budget, maxCore = 7, tlimit = 3600, ub = "MT", greedyBranching = 1, optim = "max")
  sink()
}


sol2 = list()
for(i in 1L : length(smallGapInstances))
{
  cost = smallGapInstances[[i]]$cost
  profit = smallGapInstances[[i]]$profit
  budget = smallGapInstances[[i]]$budget
  sink(paste0("sink", i, ".txt"))
  sol2[[i]] = FLSSS:::auxGAPbb(cost, profit, budget, maxCore = 7, tlimit = 3600, ub = "MT", greedyBranching = 0, optim = "max")
  sink()
}


tmp = FLSSS::GAP(maxCore = 7L, smallGapInstances[[2]]$cost, smallGapInstances[[2]]$profit, smallGapInstances[[2]]$budget, heuristic = FALSE, tlimit = 3600 * 8, threadLoad = 8L, verbose = TRUE)
save(tmp, file = "tmp.Rdata")




















