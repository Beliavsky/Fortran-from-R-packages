# include <Rcpp.h>
# include "header/BaB01kp.hpp"
using namespace Rcpp;


// List extra01knapsackBaB(
//     NumericVector weight, NumericVector value, NumericVector caps,
//     IntegerVector itemNcaps = IntegerVector(0), int maxCore = 7,
//     double tlimit = 60, String ub = "MT")


// Upon completion, best stores items that should stay in the original
// knapsack. The return value equals the revenue penalty. reassign stores
// items that should be reassigned.
// *weight is the array of penalties.
template<typename valtype, typename indtype, typename fmove>
valtype specialBiKpBaB(
    // vec<indtype> &originalTaskIndex,
    valtype *weight, valtype *value, indtype Nitem, valtype cap,
    vec<kpEle<valtype, indtype> > &Xcontain, // is of size Nitem + 2 at least.
    vec<indtype> &unitValOrder, // is of size Nitem at least.
    vec<valtype> &valuePerWeight, // is of size Nitem at least.
    vec<indtype> &stay, // is of size Nitem at least.
    vec<indtype> &reassign, // is of size Nitem at least.
    vec<indtype> &current // is of size Nitem at least.
  )
{
  indtype &Xsize = Nitem;
  indtype lenCap = Xsize + 2; // lenCap is a infinity actually. Dummy.
  // vec<kpEle<double, int> > Xcontain(Xsize + 2);
  Xcontain.resize(Xsize + 2); // Just for formality. Reallocation will never happen.
  kpEle<valtype, indtype> *X = &Xcontain[0] + 1;


  // vec<int> unitValOrder(Xsize);
  unitValOrder.resize(Xsize);
  valuePerWeight.resize(Xsize);
  {
    X[-1].accValue = 0;
    X[-1].accWeight = 0;
    for(indtype i = 0; i < Xsize; ++i)
    {
      unitValOrder[i] = i;
      valuePerWeight[i] = value[i] / weight[i];
    }
    cmp tmpCompare;
    tmpCompare.val = &valuePerWeight[0];
    std::sort(unitValOrder.begin(), unitValOrder.end(), tmpCompare);


    // Initialize X.accWeight, X.accValue, X.valuePerWeight.
    for(indtype i = 0; i < Xsize; ++i)
    {
      X[i].accWeight = X[i - 1].accWeight + weight[unitValOrder[i]];
      X[i].accValue = X[i - 1].accValue + value[unitValOrder[i]];
      X[i].valuePerWeight = valuePerWeight[unitValOrder[i]];
    }
    // ================================ transform to maximum problem
    cap = X[Xsize - 1].accWeight - cap;
    // ================================
    X[Xsize].accWeight = X[Xsize - 1].accWeight + (cap + 1);
    X[Xsize].accValue = X[Xsize - 1].accValue + 0;


    X[Xsize - 1].minWeightAfter = cap + 1;
    for(indtype i = Xsize - 2; i >= -1; --i)
    {
      X[i].minWeightAfter = std::min<valtype> (
        weight[unitValOrder[i + 1]], X[i + 1].minWeightAfter);
    }
  }


  // no time constraint.
  valtype penalty = X[Xsize - 1].accValue - bkp<valtype, indtype, fmove, false> (
    X, Xsize, cap, Xsize + 2, stay, current); // lenCap is infinite but is set to Xsize + 2.;
  reassign.resize(0);
  for(indtype k = 0, i = 0, iend = stay.size(); i < iend; )
  {
    if(k != stay[i])
    {
      // reassign.push_back(originalTaskIndex[unitValOrder[k]]);
      reassign.push_back(unitValOrder[k]);
      ++k;
    }
    else
    {
      // stay[i] = originalTaskIndex[unitValOrder[stay[i]]];
      stay[i] = unitValOrder[stay[i]];
      ++i;
    }
  }


  return penalty;
}




template<typename valtype, typename indtype>
struct WVA // weight, value, agent.
{
  // indtype nextEligibleAgent;
  valtype weight, value;
};


/*
template<typename valtype, typename indtype>
struct orderByVal
{
  bool operator() (
      const WVA<valtype, indtype> &i, const WVA<valtype, indtype> &j)
  {return i.value < j.value;}
};


// For each of the current tasks assigned to an agent, find out
// which agent would the task be assigned to if there will be
// a reassignment. There will a 0-1 task matrix.
// AT is a 2D array. AT[i][j] points to agent i, task j,
// fathom is a 2D 0-1 array. fathom[i][j] points to agent i, task j;
// fathom will become a fathomed solution.
template<typename valtype, typename indtype>
void updateFathomOneAgent(
    indtype fromWhichAgent, vec<indtype> &reassigned,
    WVA<valtype, indtype> **AT, unsigned char **fathom)
// reassigned stores tasks to be reassigned.
{
  for(indtype j = 0, jend = reassigned.size(); j < jend; ++j)
  {
    indtype task = reassigned[j];
    fathom[AT[task][fromWhichAgent].nextEligibleAgent][task] = true;
    fathom[fromWhichAgent][task] = false;
  }
}


// Return true if fathom is a valid solution.
// If it is a valid solution, the return value is the total revenue, else -1
template<typename valtype, typename indtype>
valtype examineFathomed(
    unsigned char **fathom, WVA<valtype, indtype> **AT,
    indtype Nagent, indtype Ntask, valtype *b // , vec<valtype> &agentTotalCost
                     ) // agentCosts has size Ntask
{
  bool valid = true;
  valtype totalVal = 0;
  for(indtype i = 0; i < Nagent; ++i)
  {
    valtype totalWeight = 0;
    for(indtype j = 0; j < Ntask; ++j)
    {
      if(fathom[i][j])
      {
        totalWeight += AT[i][j].weight;
        totalVal += AT[i][j].value;
      }
    }
    if(totalWeight > b[i]) valid = false;
  }
  if(valid) return totalVal;
  return -1;
}


template<typename indtype>
struct branchVariable
{
  bool isin;
  indtype agent, task;
};
*/


// 1. Make initial solution: for an unassigned task, assign the agent who has sufficient
//    budget and would generate the greatest revenue for performing this task.
// 2. Validate initial solution. If yes, a solution is found, branch ends.
// 3. Solve a series 0-1 knapsack problems. Fathom a solution.
// 4. Validate the fathomed solution. If valid, a solution is found, branch ends.
// 5. Choose the branching variable. Adjust the budget of related agents.
// 6. Go to Step 1.


template<typename valtype, typename indtype, typename fmove>
void gap(WVA<valtype, indtype> **AT, indtype Nagent, indtype Ntask,
         valtype *agentBudget)
{
  vec<vec<indtype> > agentTask(Nagent, vec<indtype> (Ntask));
  vec<vec<indtype> > agentTaskNextAgent(Nagent, vec<indtype> (Ntask));
  vec<vec<indtype> > fathomed(Nagent, vec<indtype> (Ntask));
  vec<vec<valtype> > penalty(Nagent, vec<valtype> (Ntask));
  for(indtype i = 0; i < Nagent; ++i)
  {
    agentTask[i].resize(0);
    agentTaskNextAgent[i].resize(0);
    penalty[i].resize(0);
    fathomed[i].resize(0);
  }


  vec<valtype> agentCost(Nagent, 0);


  for(indtype j = 0; j < Ntask; ++j)
  {
    indtype maxi = -1, max2i = -1;
    valtype maxval = std::numeric_limits<valtype>::min();
    valtype max2val = maxval;
    for(indtype i = 0; i < Nagent; ++i)
    {
      if(AT[i][j].weight <= agentBudget[i] and AT[i][j].value > maxval)
      {
        max2i = maxi; max2val = maxval;
        maxi = i; maxval = AT[i][j].value;
      }
    }
    agentTask[maxi].push_back(j);
    agentTaskNextAgent[maxi].push_back(max2i); // which agent j will be assigned to, if j will be assigned.
    penalty[maxi].push_back(maxval - max2val);
    agentCost[maxi] += maxval;
  }


  vec<indtype> overloadedAgent(Nagent);
  overloadedAgent.resize(0);
  for(indtype i = 0; i < Nagent; ++i)
  {
    if(agentCost[i] <= agentBudget[i]) continue;
    overloadedAgent.push_back(i);
  }


  vec<kpEle<valtype, indtype> > Xcontain(Ntask + 2); // is of size Nitem + 2 at least.
  vec<indtype> unitValOrder(Ntask); // is of size Nitem at least.
  vec<valtype> valuePerWeight(Ntask); // is of size Nitem at least.
  vec<indtype> stay(Ntask); // is of size Nitem at least.
  vec<indtype> reassign(Ntask); // is of size Nitem at least.
  vec<indtype> current(Ntask);
  vec<valtype> atweight(Ntask); // Store item weights in the knapsack.
  valtype totalPenalty = 0;


  for(indtype i = 0, iend = overloadedAgent.size(); i < iend; ++i)
  {
    indtype a = overloadedAgent[i];
    atweight.resize(agentTask[a].size());
    for(indtype k = 0, kend = atweight.size(); k < kend; ++k)
    {
      atweight[k] = AT[a][agentTask[a][k]].weight;
    }
    totalPenalty += specialBiKpBaB(
      // &agentTask[a][0],
      &atweight[0], &penalty[a][0], penalty[a].size(),
      agentCost[a] - agentBudget[a],
      Xcontain, // is of size Nitem + 2 at least.
      unitValOrder, // is of size Nitem at least.
      valuePerWeight, // is of size Nitem at least.
      stay, // is of size Nitem at least.
      reassign, // is of size Nitem at least.
      current // is of size Nitem at least.
    );


    indtype *oti = &agentTask[a][0]; // ot stands for original task index.
    for(indtype k = 0, kend = reassign.size(); k < kend; ++k)
    {
      indtype t = reassign[k];
      indtype anext = agentTaskNextAgent[a][t];
      B[anext][oti[t]] = 1;
    }
  }



































}
























