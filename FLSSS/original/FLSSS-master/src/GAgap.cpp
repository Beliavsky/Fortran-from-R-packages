// [[Rcpp::depends(RcppParallel)]]
# include <Rcpp.h>
# include <RcppParallel.h>
# include <random>
# include "header/macros.hpp"
# include "header/dnyTasking.hpp"
using namespace Rcpp;
# define unifInt std::uniform_int_distribution<indtype>
# define RNG std::mt19937_64


template<typename valtype>
struct WV
{
  valtype weight, value;
};


template<typename valtype, typename indtype>
struct sol
{
  valtype fitness, unfitness;
  vec<indtype> gene; // size = Ntask.
  vec<valtype> residualBudget;
  vec<indtype> auxContainer; // capacity = Ntask
  
  
  void initialize(indtype Nagent, indtype Ntask)
  {
    fitness = 0;
    unfitness = 0;
    gene.resize(Ntask);
    residualBudget.resize(Nagent);
    auxContainer.reserve(Ntask);
  }
  
  
  void swap(sol &x)
  {
    fitness = x.fitness;
    unfitness = x.unfitness;
    gene.swap(x.gene);
    residualBudget.swap(x.residualBudget);
  }
  
  
  void grow(WV<valtype> **info, valtype *budget)
  {
    residualBudget.assign(budget, budget + residualBudget.size());
    fitness = 0;
    for(indtype j = 0, jend = gene.size(); j < jend; ++j)
    {
      indtype a = gene[j];
      residualBudget[a] -= info[j][a].weight;
      fitness += info[j][a].value;
    }
    unfitness = 0;
    for(indtype i = 0, iend = residualBudget.size(); i < iend; ++i)
    {
      unfitness += std::max<valtype> (-residualBudget[i], 0);
    }
  }
  
  
  // update residualBudget, unfitness, fitness
  // pa, pt: previous agent and task.
  // a, t: new agent and task
  void updateAccumulativeInfo(
      indtype pa, indtype pt, indtype a, indtype t, WV<valtype> **info)
  {
    fitness += info[t][a].value - info[pt][pa].value;
    valtype tmp = residualBudget[pa] + info[pt][pa].weight;
    if(residualBudget[pa] < 0)
    {
      if(tmp >= 0) unfitness += residualBudget[pa];
      else unfitness -= info[pt][pa].weight;
    }
    residualBudget[pa] = tmp;
    tmp = residualBudget[a] - info[t][a].weight;
    if(tmp < 0)
    {
      if(residualBudget[a] < 0) unfitness += info[t][a].weight;
      else unfitness -= tmp;
    }
    residualBudget[a] = tmp;
  }
  
  
  void educateA(WV<valtype> **info, RNG &rng)
  {
    for(indtype i = 0, iend = residualBudget.size(); i < iend; ++i)
    {
      if(residualBudget[i] >= 0) continue;
      auxContainer.resize(0);
      for(indtype j = 0, jend = gene.size(); j < jend; ++j)
      {
        if(gene[j] == i) auxContainer.push_back(j);
      }
      unifInt U(0, auxContainer.size() - 1);
      indtype selectedTask = auxContainer[U(rng)];
      // current assignment is agent i, task t;
      valtype currentValue = info[selectedTask][i].value;
      WV<valtype> *col = info[selectedTask];
      indtype agentSwitchTo = i + 1;
      for(; agentSwitchTo < iend; ++agentSwitchTo)
      {
        if(col[agentSwitchTo].weight <= residualBudget[agentSwitchTo]
             and col[agentSwitchTo].value > currentValue) break;
      }
      if(agentSwitchTo == iend)
      {
        agentSwitchTo = 0;
        for(; agentSwitchTo < i; ++agentSwitchTo)
        {
          if(col[agentSwitchTo].weight <= residualBudget[agentSwitchTo]
               and col[agentSwitchTo].value > currentValue) break;
        }
      }
      if(agentSwitchTo == i) continue; // no swap
      updateAccumulativeInfo(i, selectedTask, agentSwitchTo, selectedTask, info);
      gene[selectedTask] = agentSwitchTo;
    }
  }
  
  
  void educateB(WV<valtype> **info)
  {
    indtype Ntask = gene.size(), Nagent = residualBudget.size();
    {
      valtype S = 0;
      for(int i = 0, iend = residualBudget.size(); i < iend; ++i)
      {
        if(residualBudget[i] < 0) S += -residualBudget[i];
      }
    }
    for(indtype j = 0; j < Ntask; ++j)
    {
      indtype &genej = gene[j];
      valtype currentVal = info[j][genej].value;
      for(indtype i = 0; i < Nagent; ++i)
      {
        if(i == genej) continue;
        if(info[j][i].weight > residualBudget[i] or
             currentVal >= info[i][j].value) continue;
        updateAccumulativeInfo(genej, j, i, j, info);
        genej = i;
        break;
      }
    }
  }
};


template<typename valtype, typename indtype>
struct gapGA
{
  RNG rng;
  indtype Nagent, Ntask, whichMostFitnessAndFeasible;
  valtype highestFeasibleFitness;
  vec<sol<valtype, indtype> > population;
  sol<valtype, indtype> kid;
  vec<indtype> geniusGene;
  valtype geniusFitness;
  // vec<WV<valtype> > infoContent;
  // vec<WV<valtype>*> info;
  WV<valtype> **info;
  valtype *budget; // can share.
  
  
  // return C
  void initialize(int rsd, NumericMatrix cost, NumericMatrix profit,
                  NumericVector Budget, int Npopu, String optim, WV<valtype> **INFO)
  {
    rng.seed(rsd);
    
    
    Nagent = cost.nrow();
    Ntask = cost.ncol();
    // infoContent.resize(INT(Nagent) * Ntask);
    // info.resize(Ntask);
    info = INFO;
    
    
    kid.initialize(Nagent, Ntask);
    geniusGene.resize(Ntask);
    geniusFitness = -std::numeric_limits<valtype>::max();
    
    
    budget = &Budget[0]; // implies valtype should always be double.
    
    
    population.resize(Npopu);
    unifInt U(0, Nagent - 1);
    for(indtype i = 0; i < Npopu; ++i)
    {
      population[i].initialize(Nagent, Ntask);
      for(indtype k = 0; k < Ntask; ++k)
      {
        population[i].gene[k] = U(rng);
      }
      population[i].grow(info, budget); // set parameters.
    }
    
    
    whichMostFitnessAndFeasible = -1;
    highestFeasibleFitness = -std::numeric_limits<valtype>::max();
    for(indtype i = 0, iend = population.size(); i < iend; ++i)
    {
      if(highestFeasibleFitness >= population[i].fitness or
           std::abs(population[i].unfitness) >= 1e-10) continue;
      highestFeasibleFitness = population[i].fitness;
      whichMostFitnessAndFeasible = i;
    }
  }
  
  
  indtype biTournament()
  {
    indtype psize = population.size();
    unifInt U(0, psize - 1);
    indtype i = U(rng), j = U(rng);
    if(i == j)
    {
      if(i < psize - 1) ++j;
      else
      {
        U = unifInt(0, psize - 2);
        j = U(rng);
      }
    }
    if(population[i].fitness < population[j].fitness) return j;
    else return i;
  }
  
  
  void selectParent(indtype &dadInd, indtype &mumInd)
  {
    dadInd = biTournament();
    mumInd = biTournament();
    while(dadInd == mumInd) mumInd = biTournament();
  }
  
  
  // Ntask >= 2
  void giveBirth(indtype dadInd, indtype mumInd)
  {
    unifInt U(0, Ntask - 1);
    indtype COpoint = U(rng);
    // crossover:
    sol<valtype, indtype>
      &dad = population[dadInd], &mum = population[mumInd];
      std::copy(dad.gene.begin(), dad.gene.begin() + COpoint, kid.gene.begin());
      std::copy(mum.gene.begin() + COpoint, mum.gene.end(),
                kid.gene.begin() + COpoint);
      // mutation:
      indtype geneA = U(rng);
      indtype geneB = U(rng);
      if(geneA == geneB)
      {
        if(geneA < Ntask - 1) ++geneB;
        else
        {
          U = unifInt(0, Ntask - 2);
          geneB = U(rng);
        }
      }
      std::swap(kid.gene[geneA], kid.gene[geneB]);
  }
  
  
  // Generate a kid. If it is accepted and its fitness is greater than
  // than the current best, return true, else false.
  bool generate()
  {
    indtype dadInd, mumInd;
    selectParent(dadInd, mumInd);
    giveBirth(dadInd, mumInd);
    kid.grow(&info[0], budget);
    kid.educateA(&info[0], rng);
    kid.educateB(&info[0]);
    // compare gene; loop within the population.
    indtype whichUnfitMost = std::numeric_limits<indtype>::min();
    valtype mostUnfitness = -std::numeric_limits<valtype>::max();
    indtype whichFitLeast = std::numeric_limits<indtype>::min();
    valtype leastFitness = std::numeric_limits<valtype>::max();
    indtype *kidgene = &kid.gene[0];
    for(indtype i = 0, iend = population.size(); i < iend; ++i)
    {
      if(mostUnfitness < population[i].unfitness)
      {
        mostUnfitness = population[i].unfitness;
        whichUnfitMost = i;
      }
      if(leastFitness > population[i].fitness)
      {
        leastFitness = population[i].fitness;
        whichFitLeast = i;
      }
      indtype *geneSeq = &population[i].gene[0];
      indtype j = 0;
      for(; j < Ntask; ++j)
      {
        if(geneSeq[j] != kidgene[j]) break;
      }
      if(j == Ntask) return false;
    }
    
    
    indtype swapPosition = 0;
    if(std::abs(mostUnfitness) >= 1e-10) // There are infeasible solutions in population.
    {
      swapPosition = whichUnfitMost;
    }
    else
    {
      swapPosition = whichFitLeast; // it is possible that this is the only feasible solution currently in population.
    }
    population[swapPosition].swap(kid);
    if(swapPosition == whichMostFitnessAndFeasible)
    {
      whichMostFitnessAndFeasible = -1;
      highestFeasibleFitness = -std::numeric_limits<valtype>::max();
      geniusGene.assign(population[swapPosition].gene.begin(), population[swapPosition].gene.end());
      geniusFitness = highestFeasibleFitness;
      return false;
    }
    sol<valtype, indtype> &newMember = population[swapPosition];
    if(std::abs(newMember.unfitness) < 1e-10 // feasible
         and newMember.fitness > highestFeasibleFitness)
    {
      whichMostFitnessAndFeasible = swapPosition;
      highestFeasibleFitness = newMember.fitness;
      return true; // best solution got updated.
    }
    return false;
  }
  
  
  valtype runNgenerations(INT generN, vec<indtype> &rst) // M is population size
  {
    INT Nstay = 0;
    while(true)
    {
      bool tmp = generate();
      if(tmp) Nstay = 0;
      else ++Nstay;
      if(Nstay > generN) break;
    }
    if(geniusFitness > highestFeasibleFitness)
    {
      rst.assign(geniusGene.begin(), geniusGene.end());
      return geniusFitness;
    }
    if(whichMostFitnessAndFeasible >= 0) rst.assign(
      population[whichMostFitnessAndFeasible].gene.begin(),
      population[whichMostFitnessAndFeasible].gene.end());
    return highestFeasibleFitness;
  }
};




template<typename valtype, typename indtype>
struct gapGApara: public RcppParallel::Worker
{
  INT generationN;
  gapGA<valtype, indtype> *Gvec;
  vec<vec<indtype> > &rst;
  valtype *rstVal;
  dynamicTasking *dT;
  void operator() (std::size_t st, std::size_t end)
  {
    for(;;)
    {
      INT objI = 0;
      if(!dT->nextTaskID(objI)) break;
      rstVal[objI] = Gvec[objI].runNgenerations(generationN, rst[objI]);
    }
  }
  gapGApara(INT generationN, vec<gapGA<valtype, indtype> > &G,
            vec<vec<indtype> > &rst, vec<valtype> &rstValue, std::size_t maxCore):
    generationN(generationN), rst(rst)
  {
    Gvec = &G[0];
    rst.resize(G.size(), vec<indtype>(G[0].Ntask));
    rstValue.resize(G.size());
    rstVal = &rstValue[0];
    dynamicTasking dt(maxCore, G.size()); dT = &dt;
    parallelFor(0, dT->NofCore, *this);
  }
};




// [[Rcpp::export]]
List auxGAPgaGivenRandomSeeds(
    NumericMatrix cost, NumericMatrix profitOrLoss,
    NumericVector budget, IntegerVector randomSeeds,
    int populationSize = 100, int generations = 1000,
    String optim = "max", int maxCore = 7)
{
  int trials = randomSeeds.size();
  maxCore = std::min<int> (trials, maxCore);
  vec<gapGA<double, int> > Gvec(trials);
  int Nagent = cost.nrow(), Ntask = cost.ncol();
  
  
  vec<WV<double> > infoContent(Ntask * Nagent);
  vec<WV<double>*> infoV(Ntask);
  for(int j = 0; j < Ntask; ++j)
  {
    infoV[j] = &infoContent[0] + j * Nagent;
  }
  WV<double> **I = &infoV[0];
  double C = -std::numeric_limits<double>::max();
  bool maxOpt = (optim == "max");
  if(!maxOpt) C = *std::max_element(profitOrLoss.begin(), profitOrLoss.end()) + 1;
  for(int j = 0; j < Ntask; ++j)
  {
    for(int i = 0; i < Nagent; ++i)
    {
      I[j][i].weight = cost[j * Nagent + i];
      if(maxOpt) I[j][i].value = profitOrLoss[j * Nagent + i];
      else I[j][i].value = C - profitOrLoss[j * Nagent + i];
    }
  }
  
  
  for(int i = 0, iend = Gvec.size(); i < iend; ++i)
  {
    Gvec[i].initialize(
        randomSeeds[i], cost, profitOrLoss, budget, populationSize, optim, I);
  }
  vec<vec<int> > rst;
  vec<double> rstVal;
  gapGApara<double, int> (generations, Gvec, rst, rstVal, maxCore);
  
  
  IntegerMatrix allGenes(Ntask, Gvec.size() * populationSize);
  NumericVector allProfitOrLoss(Gvec.size() * populationSize);
  NumericVector allBudgetExceedance(Gvec.size() * populationSize);
  int *v = &allGenes[0];
  int t = 0;
  for(int i = 0, iend = Gvec.size(); i < iend; ++i)
  {
    for(int j = 0, jend = Gvec[i].population.size(); j < jend; ++j)
    {
      allProfitOrLoss[t] = Gvec[i].population[j].fitness;
      if(!maxOpt) allProfitOrLoss[t] = Ntask * C - allProfitOrLoss[t];
      allBudgetExceedance[t] = Gvec[i].population[j].unfitness;
      ++t;
      for(int k = 0; k < Ntask; ++k)
      {
        v[k] = Gvec[i].population[j].gene[k] + 1;
      }
      v += Ntask;
    }
  }
  
  
  int whichMax = std::max_element(rstVal.begin(), rstVal.end()) - rstVal.begin();
  vec<int> &solution = rst[whichMax];
  double solutionRevenue = rstVal[whichMax];
  if(solutionRevenue == -std::numeric_limits<double>::max())
  {
    return List::create(Named("totalProfitOrLoss") = solutionRevenue,
                        Named("agentCost") = NumericVector(0),
                        Named("assignment") = IntegerVector(0),
                        Named("populationInfo") = List::create(
                          Named("allGenes") = allGenes,
                          Named("allBudgetExceedance") = allBudgetExceedance,
                          Named("allProfitOrLoss") = allProfitOrLoss));
  }
  
  
  if(!maxOpt)
  {
    solutionRevenue = C * Ntask - solutionRevenue;
  }
  
  
  NumericVector agentCost(Nagent);
  for(int j = 0; j < Ntask; ++j)
  {
    agentCost[solution[j]] += cost[j * Nagent + solution[j]];
  }
  
  
  IntegerVector assignment(Ntask);
  for(int i = 0; i < Ntask; ++i)
  {
    assignment[i] = solution[i] + 1;
  }
  
  
  return List::create(Named("totalProfitOrLoss") = solutionRevenue,
                      Named("agentCost") = agentCost,
                      Named("assignment") = assignment,
                      Named("populationInfo") = List::create(
                        Named("allGenes") = allGenes,
                        Named("allBudgetExceedance") = allBudgetExceedance,
                        Named("allProfitOrLoss") = allProfitOrLoss));
}








































