

makeData <- function(v)
{
  m = v[1]
  n = v[2]
  mn = m * n
  loss = t(matrix(v[3L : (mn + 2L)], nrow = n))
  cost = t(matrix(v[(mn + 2L + 1L) : (mn + 2L + mn)], nrow = n))
  budget = v[(length(v) - m + 1L) : length(v)]
  list(cost = cost, loss = loss, budget = budget)
}


gapA = lapply(paste0("tests/gapInstances/gap_a/", list.files("tests/gapInstances/gap_a")), function(x) makeData(as.integer(scan(x))))
save(gapA, file = "tests/gapInstances/gapA.Rdata")
gapB = lapply(paste0("tests/gapInstances/gap_b/", list.files("tests/gapInstances/gap_b")), function(x) makeData(as.integer(scan(x))))
save(gapB, file = "tests/gapInstances/gapB.Rdata")
gapC = lapply(paste0("tests/gapInstances/gap_c/", list.files("tests/gapInstances/gap_c")), function(x) makeData(as.integer(scan(x))))
save(gapC, file = "tests/gapInstances/gapC.Rdata")
gapD = lapply(paste0("tests/gapInstances/gap_d/", list.files("tests/gapInstances/gap_d")), function(x) makeData(as.integer(scan(x))))
save(gapD, file = "tests/gapInstances/gapD.Rdata")
gapE = lapply(paste0("tests/gapInstances/gap_e/", list.files("tests/gapInstances/gap_e")), function(x) makeData(as.integer(scan(x))))
save(gapE, file = "tests/gapInstances/gapE.Rdata")




load("tests/gapInstances/gapA.Rdata"); load("tests/gapInstances/gapB.Rdata"); load("tests/gapInstances/gapC.Rdata"); load("tests/gapInstances/gapD.Rdata"); load("tests/gapInstances/gapE.Rdata")


i = 6L; gapType = gapA
cost = gapType[[i]]$cost
profitOrLoss = gapType[[i]]$loss
budget = gapType[[i]]$budget
all(apply(cost, 2, function(x) all(x <= budget)))


sink("sink.txt")
sol = FLSSS:::auxGAPbb(cost, profitOrLoss, budget, maxCore = 7, tlimit = 10, ub = "MT", greedyBranching = 1, optim = "min")
sink()
str(sol); all(sol$agentCost <= budget)


sink("sink.txt")
sol2 = FLSSS:::auxGAPbbDp(cost, profitOrLoss, budget, maxCore = 7, tlimit = 60, greedyBranching = 1, optim = "min")
str(sol2); all(sol2$agentCost <= budget)
sink()




# Budget exceedance =
#   24, 285,
# overloadedAgentTask =
#   7, 41, 51, 87,
# 6, 14, 19, 20, 22, 24, 29, 32, 33, 38, 45, 54, 61, 69, 72, 77, 81, 96, 98, 99,
# overloadedAgentWeight =
#   11, 11, 9, 5,
# 23, 19, 25, 6, 12, 19, 14, 15, 14, 21, 17, 17, 8, 24, 14, 15, 5, 19, 17, 22,
# next agent =
#   3, 1, 3, 3,
# 3, 3, 3, 3, 3, 3, 0, 3, 0, 3, 3, 0, 3, 3, 3, 3, 0, 3, 3, 3,
# overloadedAgentPenalty =
#   12, 8, 8, 4,
# 0, 0, 0, 0, 0, 0, 4, 0, 1, 0, 0, 19, 0, 0, 0, 0, 13, 0, 0, 0,


weight = c(11, 11, 9, 5)
value = c(12, 8, 8, 4)
cap = 24
cap = sum(weight) - cap
tmp = FLSSS::auxKnapsack01bb(weight, value, cap)


weight = c(20, 12, 14, 23, 19, 23, 7, 8, 23, 11, 15, 20, 21, 19, 10, 15, 17, 20, 24, 5, 19, 12, 22, 17, 11, 12)
value = abs(c(5, 2, 9, 5, 13, 1, 6, 1, 0, 2, 0, 2, 5, 3, 1, 5, 0, 5, 0, 3, 1, 2, 3, 1, 0, 2))
cap = 54
cap = sum(weight) - cap
sink("sink.txt")
tmp = FLSSS::auxKnapsack01bb(weight, value, cap)
sink()




# overloadedAgentWeight =
#   20, 12, 14, 23, 19, 23, 7, 8, 23, 11, 15, 20, 21, 19, 10, 15, 17, 20, 24, 5, 19, 12, 22, 17, 11, 12,
# 22, 15, 12, 13, 10, 11, 11, 11, 17, 6, 8, 20, 6, 21, 11, 17, 11, 21, 22, 15, 8, 25, 15, 14, 18, 9, 9,
# 20, 23, 17, 24, 20, 20, 19, 15, 18, 5, 21, 12, 10, 24, 11, 17, 8, 16, 22, 23, 21, 16, 18, 24, 18, 22,
# overloadedAgentPenalty =
#   5, 2, 9, 5, 13, 1, 6, 1, 0, 2, 0, 2, 5, 3, 1, 5, 0, 5, 0, 3, 1, 2, 3, 1, 0, 2,
# 0, 3, 0, 5, 0, 2, 5, 5, 3, 5, 3, 4, 1, 2, 9, 2, 3, 6, 5, 5, 1, 1, 0, 0, 6, 0, 2,
# 1, 4, 16, 6, 3, 0, 1, 5, 2, 5, 1, 6, 10, 6, 6, 0, 9, 3, 6, 14, 1, 7, 7, 4, 5, 2,









