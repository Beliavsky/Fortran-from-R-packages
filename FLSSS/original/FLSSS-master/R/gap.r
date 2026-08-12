gapfindB <-function(X, Min, Max, LB = NULL, UB = NULL)
{
  nagent = length(Min) - 1L
  ntask = ncol(X) / nagent
  v = list()
  k = 1L
  for(i in seq(1L, by = nagent, len = ntask))
  {
    v[[k]] = X[, i : (i + nagent - 1L)]
    v[[k]] = rbind(v[[k]], 0L : (nagent - 1L))
    k = k + 1L
  }
  if(is.null(LB) | is.null(UB))
  {
    LB = rep(1L, ntask)
    UB = rep(nagent, ntask)
  }


  # print("========================")
  # print(Min - rowSums(as.data.frame(mapply(function(x, y)
  # {
  #   x[, y]
  # }, v, UB, SIMPLIFY = F))))
  # print("========================")
  # print("========================")
  # print(Max - rowSums(as.data.frame(mapply(function(x, y)
  # {
  #   x[, y]
  # }, v, LB, SIMPLIFY = F))))
  # print("========================")


  first = T
  while(T)
  {
    LBresv = LB
    for(i in 1L : ntask)
    {
      S = Min - rowSums(as.data.frame(mapply(function(x, y)
      {
        x[, y]
      }, v[-i], UB[-i], SIMPLIFY = F)))
      names(S) = NULL
      currentTask = v[[i]]
      # cat("Min_rowSums =", S, "\n")
      for(k in LB[i] : UB[i])
      {
        # cat("k =", k, ", ")
        # print(currentTask[, k])
        # if(all(currentTask[, k] >= S - 1e-10))
        if(currentTask[nrow(currentTask), k] >= S[nrow(currentTask)] - 1e-10)
        {
          # print("all(currentTask[, k] >= S - 1e-10)")
          break
        }
      }
      # if(k >= UB[i] & !all(currentTask[, k] >= S - 1e-10)) return(list(LB, UB, F))
      if(k >= UB[i] & !(currentTask[nrow(currentTask), k] >= S[nrow(currentTask)] - 1e-10)) return(list(LB, UB, F))
      LB[i] = k
      # cat('LB[i] = ', LB[i], '\n')
    }


    # cat("LB =", LB - 1L, "\n")


    if(all(LBresv == LB) & !first) break
    first = F


    UBresv = UB
    for(i in 1L : ntask)
    {
      S = Max - rowSums(as.data.frame(mapply(function(x, y)
      {
        x[, y]
      }, v[-i], LB[-i], SIMPLIFY = F)))
      currentTask = v[[i]]
      # cat("Max_rowSums =", S, "\n")
      for(k in UB[i] : LB[i])
      {
        if(all(currentTask[, k] <= S + 1e-10))
        {
          # print("all(currentTask[, k] <= S + 1e-10)")
          break
        }
      }
      if(k <= LB[i] & !all(currentTask[, k] <= S + 1e-10)) return(list(LB, UB, F))
      UB[i] = k
    }


    # cat("UB =", UB - 1L, "\n")


    if(all(UBresv == UB)) break
  }
  return(list(LB - 1L, UB - 1L))
}




GAPpre <- function(maxCore = 7L, agentsCosts, agentsProfits, agentsBudgets, heuristic = FALSE, tlimit = 60, threadLoad = 8L, verbose = TRUE)
{
  # Quick fail
  {
    minCosts = apply(agentsCosts, 1, function(x) min(x))
    if(!all(minCosts <= agentsBudgets))
    {
      stop("Agent budgets supressed any assignment.")
    }
    if(ncol(agentsCosts) == 1L)
    {
      stop("There is only one task.")
    }
    if(length(agentsBudgets) == 1L)
    {
      stop("There is only one agent.")
    }
  }


  profits = agentsProfits
  agents = nrow(agentsCosts)
  if(agents < 2L) return("Need at least two agents.")
  tasks = ncol(agentsCosts)
  costColumns = mapply(function(i, cost)
  {
    column = numeric(agents * tasks)
    column[seq(i, len = tasks, by = agents)] = cost
    column
  }, 1L : agents, as.data.frame(t(agentsCosts)))


  indCol = rep(0L : (agents - 1L), tasks)
  V = costColumns
  agentsIndex = indCol
  i = 1L
  while(i <= nrow(V))
  {
    rg = i : (i + agents - 1L)
    tmpOrder = order(profits[rg])
    V[rg, ] = V[rg, ][tmpOrder, ]
    agentsIndex[rg] = agentsIndex[rg][tmpOrder]
    profits[rg] = profits[rg][tmpOrder]
    i = i + agents
  }


  Vresv = V
  scaleFactor = max(agentsCosts) * 1.1


  for(i in 1L : ncol(V))
  {
    V[, i] =  V[, i] + indCol * scaleFactor
  }


  keyTarget = ((agents - 1L) * tasks) : 0L


  MAXmat = as.matrix(as.data.frame(lapply(keyTarget, function(x)
  {
    scaleFactor * x + agentsBudgets
  })))


  V = V / scaleFactor
  MAXmat = MAXmat / scaleFactor


  MAXmat = rbind(MAXmat, keyTarget)
  dimnames(MAXmat) = NULL


  return(list(V = V, MAXmat = MAXmat))
}




GAP <- function(maxCore = 7L, agentsCosts, agentsProfits, agentsBudgets, heuristic = FALSE, tlimit = 60, threadLoad = 8L, verbose = TRUE)
{
  message("This function solves the generalized assignment problem via hypercube contraction and is most likely less efficient than functions like 'auxGAPbb()'.")
  # Quick fail
  {
    minCosts = apply(agentsCosts, 1, function(x) min(x))
    if(!all(minCosts <= agentsBudgets))
    {
      stop("Agent budgets supressed any assignment.")
    }
    if(ncol(agentsCosts) == 1L)
    {
      stop("There is only one task.")
    }
    if(length(agentsBudgets) == 1L)
    {
      stop("There is only one agent.")
    }
  }


  profits = agentsProfits
  agents = nrow(agentsCosts)
  if(agents < 2L) return("Need at least two agents.")
  tasks = ncol(agentsCosts)
  costColumns = mapply(function(i, cost)
  {
    column = numeric(agents * tasks)
    column[seq(i, len = tasks, by = agents)] = cost
    column
  }, 1L : agents, as.data.frame(t(agentsCosts)))


  indCol = rep(0L : (agents - 1L), tasks)
  V = costColumns
  agentsIndex = indCol
  i = 1L
  while(i <= nrow(V))
  {
    rg = i : (i + agents - 1L)
    tmpOrder = order(profits[rg])
    V[rg, ] = V[rg, ][tmpOrder, ]
    agentsIndex[rg] = agentsIndex[rg][tmpOrder]
    profits[rg] = profits[rg][tmpOrder]
    i = i + agents
  }


  Vresv = V
  scaleFactor = max(agentsCosts) * 1.1


  for(i in 1L : ncol(V))
  {
    V[, i] =  V[, i] + indCol * scaleFactor
  }


  keyTarget = ((agents - 1L) * tasks) : 0L


  MAXmat = as.matrix(as.data.frame(lapply(keyTarget, function(x)
  {
    scaleFactor * x + agentsBudgets
  })))


  V = V / scaleFactor
  MAXmat = MAXmat / scaleFactor


  MAXmat = rbind(MAXmat, keyTarget)
  dimnames(MAXmat) = NULL


  rst = z_GAP(maxCore, t(V), profits, MAXmat, rep(0L, tasks), rep(agents - 1L, tasks), tlimit, threadLoad, verbose = verbose, heuristic = heuristic)


  rst = rst + 1L


  if(rst[1] == 1L & rst[length(rst)] == 1L) rst = NA
  if(is.na(rst[1])) foundAgent = NA
  else foundAgent = agentsIndex[rst] + 1L


  names(agentsBudgets) = NULL
  list(assignedAgents = data.frame(task = 1L : tasks, agent = foundAgent),
       assignmentProfit = sum(profits[rst]),
       assignmentCosts = colSums(Vresv[rst, ]),
       agentsBudgets = agentsBudgets,
       unconstrainedMaxProfit = sum(apply(agentsProfits, 2, function(x) max(x))),
       FLSSSsolution = rst,
       FLSSSvec = V,
       MAXmat = MAXmat,
       foreShadowFLSSSvec = Vresv)
}























