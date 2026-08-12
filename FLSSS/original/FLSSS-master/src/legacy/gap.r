



GAP <- function(maxCore = 7L, agentsCosts, agentsProfits, agentsBudgets, heuristic = FALSE, tlimit = 60, useBiSrchInFB = FALSE, threadLoad = 8L, verbose = TRUE)
{
  # quick fail
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
  V = cbind(numeric(nrow(V)), V)
  scaleFactor = c(1, apply(V[, -1], 2, function(x) max(x)))


  for(i in 1L : ncol(V))
  {
    V[, i] =  V[, i] + indCol * scaleFactor[i]
  }


  keyTarget = ((agents - 1L) * tasks) : 0L
  # keyTarget = 0L : ((agents - 1L) * tasks)
  target = c(0, agentsBudgets / 2)
  ME = c(1 / 8, agentsBudgets / 2)


  targetMat = as.data.frame(lapply(keyTarget, function(x)
  {
    scaleFactor * x + target
  }))
  colnames(targetMat) = NULL


  dlst = 0L
  dl = 1L
  dust = 0L
  du = ncol(V)


  targetMat = as.matrix(targetMat)
  dimnames(targetMat) = NULL


  rst = z_GAP(maxCore, len = ncol(agentsCosts), t(V), numeric(0), dlst, dl, dust, du, targetMat, profits, ME, zeroBasedLB = seq(0L, by = agents, len = tasks), zeroBasedUB = seq(agents - 1L, by = agents, len = tasks), tlimit, useBiSrchInFB, threadLoad, verbose, heuristic) + 1L


  if(rst[1] == 1L & rst[length(rst)] == 1L) rst = NA
  if(is.na(rst[1])) foundAgent = NA
  else foundAgent = agentsIndex[rst] + 1L


  list(assignedAgents = data.frame(task = 1L : tasks, agent = foundAgent),
       assignmentProfit = sum(profits[rst]),
       assignmentCosts = colSums(Vresv[rst, ]),
       agentsBudgets = agentsBudgets,
       unconstrainedMaxProfit = sum(apply(agentsProfits, 2, function(x) max(x))),
       FLSSSsolution = rst,
       FLSSSvec = V,
       FLSSStargets = as.matrix(targetMat),
       FLSSSme = ME,
       foreShadowFLSSSvec = Vresv)
}




GAPintegerized <- function(maxCore = 7L, agentsCosts, agentsProfits, agentsBudgets, heuristic = FALSE, precisionLevel = integer(length(agentsBudgets)), returnBeforeMining = FALSE, tlimit = 60, useBiSrchInFB = FALSE, threadLoad = 8L, verbose = TRUE)
{
  # quick fail
  {
    if(.Machine$sizeof.pointer == 4L)
    {
      message("32-bit architecture unsupported")
      return()
    }
    minCosts = apply(agentsCosts, 1, function(x) min(x))
    if(!all(minCosts <= agentsBudgets))
    {
      stop("Agent budgets supress any assignment.")
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
  mVresv = costColumns


  # integerize
  if(all(precisionLevel == 0L))
  {
    precisionLevel = rep(ncol(agentsCosts) * 8L, length(agentsBudgets))
  }


  {
    tmp = z_integerize(tasks, costColumns, agentsBudgets / 2, agentsBudgets / 2, precisionLevel)
    mV = tmp$integerized
    mTarget = tmp$target
    mME = tmp$ME
  }


  INT = list(mV = mV, mTarget = mTarget, mME = mME)


  indCol = rep(0L : (agents - 1L), tasks)
  agentsIndex = indCol
  i = 1L
  while(i <= nrow(mV)) # order costs according to profits in each agent block
  {
    rg = i : (i + agents - 1L)
    tmpOrder = order(profits[rg])
    mV[rg, ] = mV[rg, ][tmpOrder, ]
    agentsIndex[rg] = agentsIndex[rg][tmpOrder]
    profits[rg] = profits[rg][tmpOrder]
    i = i + agents
  }


  mV = cbind(integer(nrow(mV)), mV)
  scaleFactor = c(1L, apply(mV[, -1L], 2L, function(x) abs(max(x))))


  for(i in 1L : ncol(mV))
  {
    mV[, i] =  mV[, i] + indCol * scaleFactor[i]
  }


  keyTarget = ((agents - 1L) * tasks) : 0L
  # keyTarget = c(rev(keyTarget[1L : 20L]), keyTarget[-(1L : 20L)])


  mTarget = c(0L, mTarget)
  mME = c(0L, mME)


  targetMat = as.data.frame(lapply(keyTarget, function(x)
  {
    scaleFactor * x + mTarget
  }))
  targetMat = as.matrix(targetMat)
  dimnames(targetMat) = NULL


  # find the largest subset sums
  {
    tmp = seq(1L, by = agents, len = tasks)
    lowestSS = apply(mV, 2, function(x)
    {
      S = 0L
      for(i in 1L : length(tmp))
      {
        S = S + x[tmp[i]]
      }; S
    })
    largestSS = apply(mV, 2, function(x)
    {
      S = 0L
      for(i in 1L : length(tmp))
      {
        S = S + x[tmp[i] + agents - 1L]
      }; S
    })
  }


  # adjust targetMat and find max magnitude
  targetMatAndMaxMag = z_filterTargetFindLargestMagnitude(tasks, mV, targetMat, mME, lowestSS, largestSS)
  targetMat = targetMatAndMaxMag$targetMat
  maxMag = targetMatAndMaxMag$maxMag


  # return(list(mV = mV, targetMat = targetMat, mME = mME, lowestSS = lowestSS, largestSS = largestSS, maxMag = maxMag))


  # crunch integers
  {
    tmp = z_crunchIntegers(tasks, mV, targetMat, mME, maxMag = maxMag)
    mV = tmp$mV
    targetMat = tmp$targetMat
    mME = tmp$mME
    maskV = tmp$maskV
  }


  dlst = 0L
  dl = 1L
  dust = 0L
  du = ncol(mV)
  if(returnBeforeMining) return(c(INT, list(compressedDim = ncol(mV))))
  if(verbose) cat("Dimensionality reduced from", agents + 1L, "to", du, "\n")


  rst = z_GAP(maxCore, len = ncol(agentsCosts), t(mV), maskV, dlst, dl, dust, du, targetMat, profits, mME, zeroBasedLB = seq(0L, by = agents, len = tasks), zeroBasedUB = seq(agents - 1L, by = agents, len = tasks), tlimit, useBiSrchInFB, threadLoad, verbose, heuristic) + 1L


  if(rst[1] == 1L & rst[length(rst)] == 1L) rst = NA
  if(is.na(rst[1])) foundAgent = NA
  else foundAgent = agentsIndex[rst] + 1L


  list(assignedAgents = data.frame(task = 1L : tasks, agent = foundAgent),
       assignmentProfit = sum(profits[rst]),
       assignmentCosts = colSums(mVresv[agentsIndex[rst] + seq(0L, by = agents, len = tasks) + 1L, ]),
       agentsBudgets = agentsBudgets,
       unconstrainedMaxProfit = sum(apply(agentsProfits, 2, function(x) max(x))),
       FLSSSsolution = rst,
       FLSSSvec = mV,
       FLSSStargets = as.matrix(targetMat),
       FLSSSme = mME,
       foreShadowFLSSSvec = mVresv)
}



