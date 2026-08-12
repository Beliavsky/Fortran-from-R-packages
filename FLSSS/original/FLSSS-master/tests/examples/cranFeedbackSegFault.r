

for(i in 1L : 1000L){
# Test case P07 from
# https://people.sc.fsu.edu/~jburkardt/datasets/knapsack_01/knapsack_01.html
   # =====================================================================================
costs = matrix(c(70, 73, 77, 80, 82, 87, 90, 94, 98, 106, 110, 113, 115, 118, 120),
                 ncol = 1)
gains = c(135, 139, 149, 150, 156, 163, 173, 184, 192, 201, 210, 214, 221, 229, 240)
budgets = 750
# library(FLSSS)
rst2 = FLSSS::mmKnapsack(
  maxCore = 7L, len = 8L, itemsProfits = gains, itemsCosts = costs,
  capacities = budgets, heuristic = FALSE, tlimit = 3,
  threadLoad = 4L, verbose = TRUE)
}
























