



auxGAPbb <- function(cost, profitOrLoss, budget, maxCore = 7, tlimit = 60, ub = "MT", greedyBranching = TRUE, optim = "max", multhreadOn = "nodes", threadLoad = 32)
{
  if(multhreadOn == "nodes")
    auxGAPbbMulthreadNodes(cost, profitOrLoss, budget, maxCore, threadLoad, tlimit, ub, greedyBranching, optim)
  else auxGAPbbMulthreadKPs(cost, profitOrLoss, budget, maxCore, tlimit, ub, greedyBranching, optim)
}


auxGAPbbDp <- function(cost, profitOrLoss, budget, maxCore = 7, tlimit = 60, greedyBranching = TRUE, optim = "max", multhreadOn = "nodes", threadLoad = 32)
{
  if(multhreadOn == "nodes")
    auxGAPbbDpMulthreadNodes(cost, profitOrLoss, budget, maxCore, threadLoad, tlimit, greedyBranching, optim)
  else auxGAPbbDpMulthreadKPs(cost, profitOrLoss, budget, maxCore, tlimit, greedyBranching, optim)
}


auxGAPga <- function(cost, profitOrLoss, budget, trials, populationSize, generations, randomSeed = NULL, maxCore = 7, optim = "max")
{
  if(!is.null(randomSeed)) set.seed(randomSeed)
  randomSeeds = as.integer(runif(trials) * (2 ^ 31 - 2))
  auxGAPgaGivenRandomSeeds(cost, profitOrLoss, budget, randomSeeds, populationSize, generations, optim, maxCore)
}
































