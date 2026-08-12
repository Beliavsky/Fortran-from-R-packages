// [[Rcpp::depends(RcppParallel)]]
# include <Rcpp.h>
# include <RcppParallel.h>
# include "header/BaB01kp.hpp"
// # include "header/DP01kp.hpp"
using namespace Rcpp;


template<typename valtype>
struct WV
{
  valtype weight, value;
};


template<typename valtype, typename indtype>
struct stackEle
{
  indtype agent, task;
  valtype desirability;
  stackEle(){}
  stackEle(indtype agent, indtype task): agent(agent), task(task){};
  stackEle(indtype agent, indtype task, valtype desirability):
    agent(agent), task(task), desirability(desirability){};
  void assign(indtype x, indtype y, valtype z){agent = x; task = y; desirability = z;}
};


// Given that none of B[i] equals 2
// Will reset unassigned variables to -1.
template<typename valtype, typename indtype>
bool best2agent(WV<valtype> *t, signed char *Bt, indtype Nagent,
                indtype &best, indtype &best2, valtype *residualBudget)
{
  valtype bestval = -std::numeric_limits<valtype>::max(), best2val = bestval;
  best = -1; best2 = -1;
  for(indtype i = 0; i < Nagent; ++i)
  {
    if(Bt[i] == 0) continue;
    Bt[i] = -1;
    if(t[i].weight > residualBudget[i] or best2val >= t[i].value) continue;
    if(t[i].value > bestval)
    {
      best2val = bestval; best2 = best;
      bestval = t[i].value; best = i;
    }
    else
    {
      best2val = t[i].value; best2 = i;
    }
  }
  if(best == -1) return false;
  if(best2 == -1) best2 = best;
  return true;
}


template<typename valtype, typename indtype>
void popAssignedLockedStackByK(
    vec<stackEle<valtype, indtype> > &T, signed char **B,
    WV<valtype> **info, indtype Nagent, valtype *residualBudget, indtype K)
{
  for(stackEle<valtype, indtype> *x = &T.back(),
      *xend = x - K; x > xend; --x)
  {
    B[x->task][x->agent] = -1;
    residualBudget[x->agent] += info[x->task][x->agent].weight;
    B[x->task][Nagent] = 0;
  }
  T.resize(T.size() - K);
}


// Returning true means no need for backtrack.
template<typename valtype, typename indtype>
bool findOverloadedAgentsPenaltyWeightNextAgent(
    valtype &totalRevenue,
    vec<indtype> &overloaded, WV<valtype> **info, signed char **B,
    indtype Nagent, indtype Ntask,
    valtype *residualBudget, valtype *budgetExceedance,
    valtype *agentCost,
    vec<vec<indtype> > &agentTask,
    vec<vec<valtype> > &agentWeight, // of size Nagent
    vec<vec<valtype> > &agentPenalty, // will be used as values in knapsacking
    vec<vec<indtype> > &nextAgent, vec<stackEle<valtype, indtype> > &T)
{
  indtype globalStackSizeIncrease = 0;


  while(true)
  {
    indtype stackSizeIncrease = 0;
    for(indtype i = 0; i < Nagent; ++i)
    {
      nextAgent[i].resize(0);
      agentWeight[i].resize(0);
      agentPenalty[i].resize(0);
      agentTask[i].resize(0);
    }
    std::fill(agentCost, agentCost + Nagent, 0);
    totalRevenue = 0;
    for(indtype j = 0; j < Ntask; ++j)
    {
      if(B[j][Nagent])
      {
        if(stackSizeIncrease == 0)
        {
          indtype i = 0;
          while(B[j][i] != 2) ++i;
          totalRevenue += info[j][i].value;
        }
        continue;
      }
      indtype best, best2;
      bool needNoBacktrack = best2agent( // Will also set nonzero, nontwo variables to -1.
        info[j], B[j], Nagent, best, best2, residualBudget);
      if(!needNoBacktrack)
      {
        popAssignedLockedStackByK(T, B, info, Nagent, residualBudget, globalStackSizeIncrease);
        return false; // need back track.
      }
      if(best == best2) // Loop over all tasks becomes necessary because variables will be locked.
      {
        T.push_back(stackEle<valtype, indtype> (best, j, 0));
        ++stackSizeIncrease;
        ++globalStackSizeIncrease;
        B[j][best] = 2; // assign and lock.
        residualBudget[best] -= info[j][best].weight;
        B[j][Nagent] = 1; // indicate lock.
        continue;
      }
      if(stackSizeIncrease != 0) continue; // the continuing loop-over the remaining tasks is just to
                                           // check if any more best == best2 situations, and lock, update
                                           // residualBudget consequently.
      B[j][best] = 1;
      totalRevenue += info[j][best].value;
      valtype w = info[j][best].weight;
      agentWeight[best].push_back(w);
      agentCost[best] += w;
      nextAgent[best].push_back(best2);
      agentPenalty[best].push_back(info[j][best].value - info[j][best2].value);
      agentTask[best].push_back(j);
    }
    if(stackSizeIncrease == 0) break;
  }


  overloaded.resize(0);
  for(indtype i = 0; i < Nagent; ++i)
  {
    budgetExceedance[i] = agentCost[i] - residualBudget[i];
    if(budgetExceedance[i] <= 0) continue;
    if(agentWeight[i].size() == 0)
    {
      popAssignedLockedStackByK(T, B, info, Nagent, residualBudget, globalStackSizeIncrease);
      return false;
    }
    overloaded.push_back(i);
  }
  return true;
}




template<typename valtype, typename indtype>
bool thereIsOverlodedAgent(WV<valtype> **info,
                           signed char **B, indtype Nagent, indtype Ntask,
                           valtype *residualBudget, valtype *agentCost,
                           valtype &totalRevenue)
{
  totalRevenue = 0;
  std::fill(agentCost, agentCost + Nagent, 0);
  for(indtype j = 0; j < Ntask; ++j)
  {
    for(indtype i = 0; i < Nagent; ++i)
    {
      if(B[j][i] == 1)
      {
        agentCost[i] += info[j][i].weight;
        totalRevenue += info[j][i].value;
      }
      else if(B[j][i] == 2) totalRevenue += info[j][i].value;
    }
  }
  for(indtype i = 0; i < Nagent; ++i)
  {
    if(agentCost[i] > residualBudget[i]) return true;
  }
  return false;
}




template<typename valtype, typename indtype, typename fmove>
valtype specialBiKpBaB(
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
  // indtype lenCap = Xsize + 2; // lenCap is a infinity actually. Dummy.
  // vec<kpEle<double, int> > Xcontain(Xsize + 2);
  Xcontain.resize(Xsize + 3); // Just for formality. Reallocation will never happen.
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
    // cmp<valtype, indtype> tmpCompare;
    // tmpCompare.val = &valuePerWeight[0];
    // std::sort(unitValOrder.begin(), unitValOrder.end(), tmpCompare);
    std::sort(unitValOrder.begin(), unitValOrder.end(),
              cmp<valtype, indtype> (&valuePerWeight[0]));


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
    // X[Xsize].accWeight = X[Xsize - 1].accWeight + (cap + 1);
    // X[Xsize].accValue = X[Xsize - 1].accValue + 0;


    X[Xsize].accWeight = X[Xsize - 1].accWeight + (cap + 1);
    X[Xsize].valuePerWeight = X[Xsize - 1].valuePerWeight / 2;
    X[Xsize].accValue = X[Xsize - 1].accValue + X[Xsize].valuePerWeight * (cap + 1);
    X[Xsize].minWeightAfter = cap + 2;


    X[Xsize + 1].accWeight = X[Xsize].accWeight + (cap + 2);
    X[Xsize + 1].valuePerWeight = X[Xsize].valuePerWeight / 2;
    X[Xsize + 1].accValue = X[Xsize].accValue + X[Xsize + 1].valuePerWeight * (cap + 2);
    X[Xsize + 1].minWeightAfter = cap + 3;


    X[Xsize - 1].minWeightAfter = cap + 1;
    for(indtype i = Xsize - 2; i >= -1; --i)
    {
      X[i].minWeightAfter = std::min<valtype> (
        weight[unitValOrder[i + 1]], X[i + 1].minWeightAfter);
    }
  }


  // no time constraint.
  valtype penalty = X[Xsize - 1].accValue - bkp<valtype, indtype, fmove, false> (
    X, Xsize, cap, Xsize, stay, current, std::numeric_limits<double>::max());
  // lenCap is infinite but is set to Xsize + 2.;


  // Rcout << "Inside knapsacking, stay = ";
  reassign.resize(0);
  indtype k = 0;
  for(indtype i = 0, iend = stay.size(); i < iend; ++k)
  {
    if(k != stay[i])
    {
      reassign.push_back(unitValOrder[k]);
    }
    else
    {
      stay[i] = unitValOrder[k];
      // Rcout << stay[i] << ", ";
      ++i;
    }
  }
  // Rcout << "\n";
  for(; k < Nitem; ++k) reassign.push_back(unitValOrder[k]);
  return penalty;
}




// return totalPenalty.
template<typename valtype, typename indtype, typename fmove>
valtype specialBiKpBaBmulti(
    vec<indtype> &overloadedAgent,
    vec<vec<valtype> > &weight, vec<vec<valtype> > &penalty,
    vec<vec<indtype> > &stay, vec<vec<indtype> > &reassign,
    vec<valtype> &budgetExceedance,
    vec<valtype> &valuePerWeight, // aux
    vec<kpEle<valtype, indtype> > &Xcontain, // aux
    vec<indtype> &unitValOrder, // aux
    vec<indtype> &current)
{
  valtype totalPenalty = 0;
  for(int i = 0, iend = overloadedAgent.size(); i < iend; ++i)
  {
    indtype a = overloadedAgent[i];
    totalPenalty += specialBiKpBaB<valtype, indtype, fmove> (
      &weight[a][0], &penalty[a][0], penalty[a].size(),
      budgetExceedance[a], Xcontain, unitValOrder,
      valuePerWeight, stay[i], reassign[i], current);
  }
  return totalPenalty;
}




// Update B and find the branching variable
template<typename indtype> // Update B by reassigned
void updateBafterKnapsacking(
    vec<indtype> &overloaded, signed char **B,
    vec<vec<indtype> > &nextAgent,
    vec<vec<indtype> > &targetTask,
    vec<vec<indtype> > &reassigned) // reassigned corresponds to overloaded
{
  for(indtype i = 0, iend = overloaded.size(); i < iend; ++i)
  {
    indtype a = overloaded[i];
    vec<indtype> &agentSwapTo = nextAgent[a], &taskToReassign = targetTask[a];
    for(indtype k = 0, kend = reassigned[i].size(); k < kend; ++k)
    {
      indtype tmp = reassigned[i][k];
      indtype tmp2 = taskToReassign[tmp];
      B[tmp2][a] = -1;
      B[tmp2][agentSwapTo[tmp]] = 1;
    }
  }
}




template<typename valtype, typename indtype, bool greedyBranch> // Update by stay.
void pushAllBranchingVariableIntoStack(
    vec<stackEle<valtype, indtype> > &T,
    signed char **B, indtype Nagent,
    vec<indtype> &overloaded, vec<vec<indtype> > &stay,
    vec<vec<indtype> > &targetTask,
    vec<vec<valtype> > &weight, // overloaded agent, task weights
    vec<vec<valtype> > &penalty, // overloaded agent, task values
    valtype *residualBudget, WV<valtype> **info)
{
  // Push all elements in stay into stack.
  stackEle<valtype, indtype> *Tst = &*T.end();
  for(indtype i = 0, iend = overloaded.size(); i < iend; ++i)
  {
    indtype a = overloaded[i];
    for(indtype k = 0, kend = stay[i].size(); k < kend; ++k)
    {
      indtype tmp = stay[i][k];
      valtype desirability = penalty[a][tmp] / weight[a][tmp] * residualBudget[a];
      indtype tmpTask = targetTask[a][tmp];
      T.push_back(stackEle<valtype, indtype> (a, tmpTask, desirability));
      if(greedyBranch)
      {
        residualBudget[a] -= weight[a][tmp];
        B[tmpTask][a] = 2;
        B[tmpTask][Nagent] = 1;
      }
      for(stackEle<valtype, indtype> *t = &T.back() - 1; t >= Tst; --t)
      {
        if(t->desirability >= (t + 1)->desirability) break;
        std::swap(t[0], t[1]);
      }
    }
  }


  if(greedyBranch) return;


  // Check if all newly pushed elements in T are appropriate, and if not, pop.
  // This check step is not necessary, but the branch tree may be shaped better.
  {
    indtype a = Tst->agent, t = Tst->task;
    residualBudget[a] -= info[t][a].weight;
    B[t][a] = 2;
    B[t][Nagent] = 1;
  }
  indtype i = 1;
  for(indtype iend = &*T.end() - Tst; i < iend; ++i)
  {
    indtype a = Tst[i].agent, t = Tst[i].task;
    valtype tmpResidualBudget = residualBudget[a] - info[t][a].weight;
    valtype *w = &weight[a][0];
    indtype *ts = &targetTask[a][0];
    bool stackStop = false;
    for(indtype j = 0, jend = weight[a].size(); j < jend; ++j)
    {
      if(ts[j] == t or w[j] <= tmpResidualBudget) continue;
      stackStop = true;
      break;
    }
    if(stackStop) break;
    residualBudget[a] = tmpResidualBudget;
    B[t][a] = 2;
    B[t][Nagent] = 1;
  }
  T.resize(Tst - &*T.begin() + i);
  /*
  Rcout << "stack increased by = " << T.size() - (Tst - &*T.begin()) << "\n";
  Rcout << "After stack push:\n";
  for(indtype t = 1, tend = T.size(); t < tend; ++t)
  {
    Rcout << T[t].agent << ", " << T[t].task << ", " << int(B[T[t].task][T[t].agent]) << ",     ";
  }
  Rcout << "\n";
  */
}




// // Backtrack
// returning false means backtrach fails.
template<typename valtype, typename indtype>
bool backtrack(vec<stackEle<valtype, indtype> > &T, signed char **B,
               indtype Nagent, WV<valtype> **info,
               valtype *residualBudget)
{
  while(T.size() > 0 and B[T.back().task][T.back().agent] == 0)
  {
    B[T.back().task][T.back().agent] = -1;
    T.pop_back();
  }
  if(T.size() == 0) return false;
  indtype t = T.back().task, a = T.back().agent;
  B[t][a] = 0;
  B[t][Nagent] = 0;
  residualBudget[a] += info[t][a].weight;
  return true;
}


// B has been initialized. Elements are either -1 or 1.
// After this function, if stackDepthLim == -1, it means the function returns because
// T.size() - firstZeroN > stackDepthLim
template<typename valtype, typename indtype, typename fmove, bool greedyBranch>
valtype gapBabStageI(vec<signed char> &currentSolution, vec<signed char> &Bcontainer,
               indtype Nagent, indtype Ntask, WV<valtype> **info,
               valtype *residualBudget, std::time_t timer, double tlimit,
               int &Nnode, int &Nkp, indtype &stackDepthLim,
               vec<stackEle<valtype, indtype> > &T)
{
  vec<signed char*> Bv(Ntask);
  for(indtype j = 0; j < Ntask; ++j)
  {
    INT tmp = j * (INT(Nagent) + 1);
    Bv[j] = &Bcontainer[0] + tmp;
    Bv[j][Nagent] = 0;
  }
  signed char **B = &Bv[0];


  vec<indtype> overloadedAgent(Nagent);
  vec<vec<indtype> > overloadedAgentTask(Nagent, vec<indtype>(Ntask));
  vec<vec<valtype> > overloadedAgentWeight(Nagent, vec<valtype>(Ntask));
  vec<vec<valtype> > overloadedAgentPenalty(Nagent, vec<valtype>(Ntask)); // will be used as values in knapsacking
  vec<vec<indtype> > nextAgent(Nagent, vec<indtype>(Ntask));
  vec<vec<indtype> > reassign(Nagent, vec<indtype>(Ntask));
  vec<vec<indtype> > stay(Nagent, vec<indtype>(Ntask));


  vec<valtype> budgetExceedance(Nagent);
  T.reserve(INT(Nagent) * Ntask);
  currentSolution.resize(INT(Nagent) * Ntask);
  valtype currentSolutionRevenue = -std::numeric_limits<valtype>::max();


  // Auxiliary containers.
  // maxCore = std::min<int> (maxCore, Nagent);
  vec<kpEle<valtype, indtype> > Xcontain(Ntask + 3); // aux
  vec<valtype> valuePerWeight(Ntask);
  vec<indtype> unitValOrder(Ntask); // aux
  vec<indtype> current(Ntask); // aux
  vec<valtype> agentCosts(Nagent);


  // bool postKnapsack = false;
  while(true)
  {
    std::time_t now; std::time(&now);
    if(std::difftime(now, timer) > tlimit) break;
    while(true) // Repeat backtracking until knapsacking becomes necessary.
    {
      valtype revenueUB = 0;
      bool needNoBacktrack = findOverloadedAgentsPenaltyWeightNextAgent(
        revenueUB, overloadedAgent, info, B,
        Nagent, Ntask, residualBudget, &budgetExceedance[0], &agentCosts[0],
        overloadedAgentTask, overloadedAgentWeight, // of size Nagent
        overloadedAgentPenalty, // will be used as values in knapsacking
        nextAgent, T);
      // Rcout << "needNoBacktrack = " << needNoBacktrack << "\n";
      if(needNoBacktrack and revenueUB > currentSolutionRevenue)
      {
        if(overloadedAgent.size() > 0) break;
        currentSolution.assign(Bcontainer.begin(), Bcontainer.end());
        currentSolutionRevenue = revenueUB;
      }
      /*
      Rcout << "revenueUB = " << revenueUB << "\n";
      Rcout << "before backtrack, T.size() = " << T.size() << "\n";
      for(indtype i = 1, iend = T.size(); i < iend; ++i)
      {
      Rcout << T[i].agent << ", " << T[i].task << ", " <<
        int(B[T[i].task][T[i].agent]) << ",    ";
      }
      Rcout << "\n";
      for(indtype i = 0; i < Nagent; ++i)
      {
      for(indtype j = 0; j < Ntask; ++j)
      {
      if(B[j][i] >= 0) Rcout << " " << int(B[j][i]) << ", ";
      else Rcout << int(B[j][i]) << ", ";
      }
      Rcout << "\n";
      }
      */


      bool bt = backtrack(T, B, Nagent, info, residualBudget);
      /*
      Rcout << "after backtrack, T.size() = " << T.size() << "\n";
      for(indtype i = 1, iend = T.size(); i < iend; ++i)
      {
      Rcout << T[i].agent << ", " << T[i].task << ", " <<
        int(B[T[i].task][T[i].agent]) << ",    ";
      }
      Rcout << "\n";
      for(indtype i = 0; i < Nagent; ++i)
      {
      for(indtype j = 0; j < Ntask; ++j)
      {
      if(B[j][i] >= 0) Rcout << " " << int(B[j][i]) << ", ";
      else Rcout << int(B[j][i]) << ", ";
      }
      Rcout << "\n";
      }*/
      if(!bt) return currentSolutionRevenue;
    }


    /*
    {
    Rcout << "Current max profit = " << currentSolutionRevenue << "\n";
    Rcout << "After initialization, B = \n";
    for(indtype i = 0; i < Nagent; ++i)
    {
    for(indtype j = 0; j < Ntask; ++j)
    {
    if(B[j][i] >= 0) Rcout << " " << int(B[j][i]) << ", ";
    else Rcout << int(B[j][i]) << ", ";
    }
    Rcout << "\n";
    }
    Rcout << "Overloaded agent = \n";
    for(indtype i = 0, iend = overloadedAgent.size(); i < iend; ++i)
    {
    Rcout << overloadedAgent[i] << ", ";
    }
    Rcout << "\n";
    Rcout << "Budget exceedance = \n";
    for(indtype i = 0, iend = overloadedAgent.size(); i < iend; ++i)
    {
    Rcout << budgetExceedance[overloadedAgent[i]] << ", ";
    }
    Rcout << "\n";
    Rcout << "overloadedAgentTask = \n";
    for(indtype i = 0, iend = overloadedAgent.size(); i < iend; ++i)
    {
    indtype a = overloadedAgent[i];
    for(indtype j = 0, jend = overloadedAgentTask[a].size(); j < jend; ++j)
    {
    Rcout << overloadedAgentTask[a][j] << ", ";
    }
    Rcout << "\n";
    }
    Rcout << "overloadedAgentWeight = \n";
    for(indtype i = 0, iend = overloadedAgent.size(); i < iend; ++i)
    {
    indtype a = overloadedAgent[i];
    for(indtype j = 0, jend = overloadedAgentWeight[a].size(); j < jend; ++j)
    {
    Rcout << overloadedAgentWeight[a][j] << ", ";
    }
    Rcout << "\n";
    }
    Rcout << "overloadedAgentPenalty = \n";
    for(indtype i = 0, iend = overloadedAgent.size(); i < iend; ++i)
    {
    indtype a = overloadedAgent[i];
    for(indtype j = 0, jend = overloadedAgentPenalty[a].size(); j < jend; ++j)
    {
    Rcout << overloadedAgentPenalty[a][j] << ", ";
    }
    Rcout << "\n";
    }
    Rcout << "next agent = \n";
    for(indtype i = 0, iend = overloadedAgent.size(); i < iend; ++i)
    {
    indtype a = overloadedAgent[i];
    for(indtype j = 0, jend = nextAgent[a].size(); j < jend; ++j)
    {
    Rcout << nextAgent[a][j] << ", ";
    }
    Rcout << "\n";
    }
    }
    */


    ++Nnode;
    Nkp += overloadedAgent.size();
    // valtype totalPenalty =
    specialBiKpBaBmulti<valtype, indtype, fmove> (
        overloadedAgent, overloadedAgentWeight, overloadedAgentPenalty,
        stay, reassign, budgetExceedance, valuePerWeight, Xcontain, unitValOrder, current);


    // Set variables pointed by reassigned to -1.
    updateBafterKnapsacking<indtype> (overloadedAgent, B, nextAgent, overloadedAgentTask, reassign);


    /*
    {
    Rcout << "After KP, stay = \n";
    for(indtype i = 0, iend = overloadedAgent.size(); i < iend; ++i)
    {
    for(indtype j = 0, jend = stay[i].size(); j < jend; ++j)
    {
    Rcout << stay[i][j] << ", ";
    }
    Rcout << "\n";
    }
    Rcout << "After KP, reassign = \n";
    for(indtype i = 0, iend = overloadedAgent.size(); i < iend; ++i)
    {
    for(indtype j = 0, jend = reassign[i].size(); j < jend; ++j)
    {
    Rcout << reassign[i][j] << ", ";
    }
    Rcout << "\n";
    }
    // Rcout << "After KP, B = \n";
    // for(indtype i = 0; i < Nagent; ++i)
    // {
    //   for(indtype j = 0; j < Ntask; ++j)
    //   {
    //     if(B[j][i] >= 0) Rcout << " " << int(B[j][i]) << ", ";
    //     else Rcout << int(B[j][i]) << ", ";
    //   }
    //   Rcout << "\n";
    // }
    }
    */


    valtype totalReve = 0;
    bool thereis = thereIsOverlodedAgent(
      info, B, Nagent, Ntask, residualBudget, &agentCosts[0], totalReve);


    /*
    Rcout << "Knapsack solution, totalReve = " << totalReve << "\n";
    if(thereis) Rcout << "There is overloaded agent\n";
    else Rcout << "There is no overloaded agent\n";
    */


    if(!thereis)
    {
      if(totalReve > currentSolutionRevenue)
      {
        currentSolutionRevenue = totalReve;
        currentSolution.assign(Bcontainer.begin(), Bcontainer.end());
      }
      bool bt = backtrack(T, B, Nagent, info, residualBudget);
      if(!bt) break;
    }
    else
    {
      pushAllBranchingVariableIntoStack<valtype, indtype, greedyBranch> (
          T, B, Nagent, overloadedAgent, stay, overloadedAgentTask,
          overloadedAgentWeight, overloadedAgentPenalty, residualBudget, info);
      indtype firstZeroN = 0;
      while(B[T[firstZeroN].task][T[firstZeroN].agent] == 0) ++firstZeroN;
      if(T.size() - firstZeroN > INT(stackDepthLim))
      {
        stackDepthLim = -1;
        return currentSolutionRevenue;
      }
    }


    /*
     {
     Rcout << "After KP locking, stack = \n";
     for(indtype i = 0, iend = T.size(); i < iend; ++i)
     {
     Rcout << T[i].agent << ", " << T[i].task << ", " << int(B[T[i].task][T[i].agent]) << "    ";
     }
     Rcout << "\n\n\n";
     // Rcout << "After KP locking, B = \n";
     // for(indtype i = 0; i < Nagent; ++i)
     // {
     //   for(indtype j = 0; j < Ntask; ++j)
     //   {
     //     if(B[j][i] >= 0) Rcout << " " << int(B[j][i]) << ", ";
     //     else Rcout << int(B[j][i]) << ", ";
     //   }
     //   Rcout << "\n";
     // }
     }
     */
  }


  return currentSolutionRevenue;
  }




template<typename valtype, typename indtype>
void unstack(vec<signed char> &B, indtype Nagent,
             vec<vec<signed char> > &Bcontainers,
             vec<stackEle<valtype, indtype> > &T,
             indtype stackSizeWanted)
{
  for(indtype i = stackSizeWanted, iend = T.size(); i < iend; ++i)
  {
    indtype t = T[i].task, a = T[i].agent;
    if(B[t * (Nagent + 1) + a] == 2) B[t * (Nagent + 1) + Nagent] = 0;
    B[t * (Nagent + 1) + a] = -1;
  }
  T.resize(stackSizeWanted);


  Bcontainers.resize(INT(1) << T.size());
  Bcontainers[0].assign(B.begin(), B.end());
  indtype i = 1;
  for(indtype iend = Bcontainers.size(); i < iend; ++i)
  {
    // Find the T's last element that is 2, set it 0 and set all subsequent elements to 2.
    indtype t = T.size() - 1;
    for(; t >= 0; --t)
    {
      indtype tk = T[t].task, a = T[t].agent;
      signed char &bta = B[tk * (Nagent + 1) + a];
      signed char &bta2 = B[tk * (Nagent + 1) + Nagent];
      if(bta == 2)
      {
        bta = 0; bta2 = 0;
        break;
      }
      else
      {
        bta = 2; bta2 = 1;
      }
    }
    if(t < 0) break;
    Bcontainers[i].assign(B.begin(), B.end());
  }
  Bcontainers.resize(i);
}




// Plan: there are only maxCore gap objects in memory.
template<typename valtype, typename indtype, typename fmove, bool greedyBranch>
struct gapOBJ
{
  // vec<signed char> *Bctn;
  vec<valtype> residualBudget;
  vec<signed char*> Bv;
  vec<indtype> overloadedAgent;
  vec<vec<indtype> > overloadedAgentTask;
  vec<vec<valtype> > overloadedAgentWeight;
  vec<vec<valtype> > overloadedAgentPenalty; // will be used as values in knapsacking
  vec<vec<indtype> > nextAgent;
  vec<vec<indtype> > reassign;
  vec<vec<indtype> > stay;
  vec<valtype> budgetExceedance;
  vec<stackEle<valtype, indtype> > T;
  vec<kpEle<valtype, indtype> > Xcontain; // aux
  vec<valtype> valuePerWeight;
  vec<indtype> unitValOrder; // aux
  vec<indtype> current; // aux
  vec<valtype> agentCosts; // aux


  void initialize(indtype Nagent, indtype Ntask)
  {
    Bv.resize(Ntask);
    residualBudget.resize(Nagent); // has not been calculated.
    overloadedAgent.reserve(Nagent);
    overloadedAgentTask.resize(Nagent, vec<indtype>(Ntask));
    overloadedAgentWeight.resize(Nagent, vec<valtype>(Ntask));
    overloadedAgentPenalty.resize(Nagent, vec<valtype>(Ntask));
    nextAgent.resize(Nagent, vec<indtype>(Ntask));
    reassign.resize(Nagent, vec<indtype>(Ntask));
    stay.resize(Nagent, vec<indtype>(Ntask));
    budgetExceedance.resize(Nagent);
    T.reserve(INT(Nagent) * Ntask);
    Xcontain.reserve(Ntask + 3); // aux
    valuePerWeight.reserve(Ntask);
    unitValOrder.reserve(Ntask); // aux
    current.reserve(Ntask); // aux
    agentCosts.resize(Nagent); // aux
  }


  // Return false if time is up.
  void run(WV<valtype> **info, vec<signed char> &Bctn,
           std::time_t timer, double tlimit, tbb::spin_mutex *mtx,
           vec<signed char> *currentSolution,
           valtype *currentSolutionRevenue, valtype *budget,
           indtype &Nnode, indtype &Nkp)
  {
    indtype Nagent = residualBudget.size(), Ntask = Bv.size();
    Nnode = 0;
    Nkp = 0;
    vec<signed char> &Bcontainer = Bctn;
    for(indtype j = 0; j < Ntask; ++j)
    {
      INT tmp = j * (INT(Nagent) + 1);
      Bv[j] = &Bcontainer[0] + tmp;
    }
    signed char **B = &Bv[0];
    std::copy(budget, budget + Nagent, residualBudget.begin());
    for(indtype j = 0; j < Ntask; ++j)
    {
      if(!B[j][Nagent]) continue;
      for(indtype i = 0; i < Nagent; ++i)
      {
        if(B[j][i] != 2) continue;
        residualBudget[i] -= info[j][i].weight;
      }
    }
    /*
    // vec<indtype> overloadedAgent(Nagent);
    // vec<vec<indtype> > overloadedAgentTask(Nagent, vec<indtype>(Ntask));
    // vec<vec<valtype> > overloadedAgentWeight(Nagent, vec<valtype>(Ntask));
    // vec<vec<valtype> > overloadedAgentPenalty(Nagent, vec<valtype>(Ntask)); // will be used as values in knapsacking
    // vec<vec<indtype> > nextAgent(Nagent, vec<indtype>(Ntask));
    // vec<vec<indtype> > reassign(Nagent, vec<indtype>(Ntask));
    // vec<vec<indtype> > stay(Nagent, vec<indtype>(Ntask));
    //
    //
    // vec<valtype> budgetExceedance(Nagent);
    // T.reserve(INT(Nagent) * Ntask);
    // currentSolution.resize(INT(Nagent) * Ntask);
    // valtype currentSolutionRevenue = -std::numeric_limits<valtype>::max();


    // Auxiliary containers.
    // maxCore = std::min<int> (maxCore, Nagent);
    // vec<kpEle<valtype, indtype> > Xcontain(Ntask + 2); // aux
    // vec<valtype> valuePerWeight(Ntask);
    // vec<indtype> unitValOrder(Ntask); // aux
    // vec<indtype> current(Ntask); // aux
    // vec<valtype> agentCosts(Nagent);
    */



    while(true)
    {
      std::time_t now; std::time(&now);
      if(std::difftime(now, timer) > tlimit) break;


      while(true) // Repeat backtracking until knapsacking becomes necessary.
      {
        valtype revenueUB = 0;
        bool needNoBacktrack = findOverloadedAgentsPenaltyWeightNextAgent(
          revenueUB, overloadedAgent, info, B,
          Nagent, Ntask, &residualBudget[0], &budgetExceedance[0], &agentCosts[0],
          overloadedAgentTask, overloadedAgentWeight, // of size Nagent
          overloadedAgentPenalty, // will be used as values in knapsacking
          nextAgent, T);
        // Rcout << "needNoBacktrack = " << needNoBacktrack << "\n";
        if(needNoBacktrack and revenueUB > *currentSolutionRevenue)
        {
          if(overloadedAgent.size() > 0) break;
          mtx->lock();
          currentSolution->assign(Bcontainer.begin(), Bcontainer.end());
          *currentSolutionRevenue = revenueUB;
          mtx->unlock();
        }
        /*
        Rcout << "revenueUB = " << revenueUB << "\n";
        Rcout << "before backtrack, T.size() = " << T.size() << "\n";
        for(indtype i = 1, iend = T.size(); i < iend; ++i)
        {
        Rcout << T[i].agent << ", " << T[i].task << ", " <<
          int(B[T[i].task][T[i].agent]) << ",    ";
        }
        Rcout << "\n";
        for(indtype i = 0; i < Nagent; ++i)
        {
        for(indtype j = 0; j < Ntask; ++j)
        {
        if(B[j][i] >= 0) Rcout << " " << int(B[j][i]) << ", ";
        else Rcout << int(B[j][i]) << ", ";
        }
        Rcout << "\n";
        }
        */


        bool bt = backtrack(T, B, Nagent, info, &residualBudget[0]);
        /*
        Rcout << "after backtrack, T.size() = " << T.size() << "\n";
        for(indtype i = 1, iend = T.size(); i < iend; ++i)
        {
        Rcout << T[i].agent << ", " << T[i].task << ", " <<
          int(B[T[i].task][T[i].agent]) << ",    ";
        }
        Rcout << "\n";
        for(indtype i = 0; i < Nagent; ++i)
        {
        for(indtype j = 0; j < Ntask; ++j)
        {
        if(B[j][i] >= 0) Rcout << " " << int(B[j][i]) << ", ";
        else Rcout << int(B[j][i]) << ", ";
        }
        Rcout << "\n";
        }*/
        // if(!bt) return currentSolutionRevenue;
        if(!bt) return;
      }


      /*
      {
      Rcout << "Current max profit = " << currentSolutionRevenue << "\n";
      Rcout << "After initialization, B = \n";
      for(indtype i = 0; i < Nagent; ++i)
      {
      for(indtype j = 0; j < Ntask; ++j)
      {
      if(B[j][i] >= 0) Rcout << " " << int(B[j][i]) << ", ";
      else Rcout << int(B[j][i]) << ", ";
      }
      Rcout << "\n";
      }
      Rcout << "Overloaded agent = \n";
      for(indtype i = 0, iend = overloadedAgent.size(); i < iend; ++i)
      {
      Rcout << overloadedAgent[i] << ", ";
      }
      Rcout << "\n";
      Rcout << "Budget exceedance = \n";
      for(indtype i = 0, iend = overloadedAgent.size(); i < iend; ++i)
      {
      Rcout << budgetExceedance[overloadedAgent[i]] << ", ";
      }
      Rcout << "\n";
      Rcout << "overloadedAgentTask = \n";
      for(indtype i = 0, iend = overloadedAgent.size(); i < iend; ++i)
      {
      indtype a = overloadedAgent[i];
      for(indtype j = 0, jend = overloadedAgentTask[a].size(); j < jend; ++j)
      {
      Rcout << overloadedAgentTask[a][j] << ", ";
      }
      Rcout << "\n";
      }
      Rcout << "overloadedAgentWeight = \n";
      for(indtype i = 0, iend = overloadedAgent.size(); i < iend; ++i)
      {
      indtype a = overloadedAgent[i];
      for(indtype j = 0, jend = overloadedAgentWeight[a].size(); j < jend; ++j)
      {
      Rcout << overloadedAgentWeight[a][j] << ", ";
      }
      Rcout << "\n";
      }
      Rcout << "overloadedAgentPenalty = \n";
      for(indtype i = 0, iend = overloadedAgent.size(); i < iend; ++i)
      {
      indtype a = overloadedAgent[i];
      for(indtype j = 0, jend = overloadedAgentPenalty[a].size(); j < jend; ++j)
      {
      Rcout << overloadedAgentPenalty[a][j] << ", ";
      }
      Rcout << "\n";
      }
      Rcout << "next agent = \n";
      for(indtype i = 0, iend = overloadedAgent.size(); i < iend; ++i)
      {
      indtype a = overloadedAgent[i];
      for(indtype j = 0, jend = nextAgent[a].size(); j < jend; ++j)
      {
      Rcout << nextAgent[a][j] << ", ";
      }
      Rcout << "\n";
      }
      }
      */


      ++Nnode;
      Nkp += overloadedAgent.size();
      // valtype totalPenalty =
      specialBiKpBaBmulti<valtype, indtype, fmove> (
          overloadedAgent, overloadedAgentWeight, overloadedAgentPenalty,
          stay, reassign, budgetExceedance, valuePerWeight,
          Xcontain, unitValOrder, current);


      // Set variables pointed by reassigned to -1.
      updateBafterKnapsacking<indtype> (
          overloadedAgent, B, nextAgent, overloadedAgentTask, reassign);


      /*
      {
      Rcout << "After KP, stay = \n";
      for(indtype i = 0, iend = overloadedAgent.size(); i < iend; ++i)
      {
      for(indtype j = 0, jend = stay[i].size(); j < jend; ++j)
      {
      Rcout << stay[i][j] << ", ";
      }
      Rcout << "\n";
      }
      Rcout << "After KP, reassign = \n";
      for(indtype i = 0, iend = overloadedAgent.size(); i < iend; ++i)
      {
      for(indtype j = 0, jend = reassign[i].size(); j < jend; ++j)
      {
      Rcout << reassign[i][j] << ", ";
      }
      Rcout << "\n";
      }
      // Rcout << "After KP, B = \n";
      // for(indtype i = 0; i < Nagent; ++i)
      // {
      //   for(indtype j = 0; j < Ntask; ++j)
      //   {
      //     if(B[j][i] >= 0) Rcout << " " << int(B[j][i]) << ", ";
      //     else Rcout << int(B[j][i]) << ", ";
      //   }
      //   Rcout << "\n";
      // }
      }
      */


      valtype totalReve = 0;
      bool thereis = thereIsOverlodedAgent(
        info, B, Nagent, Ntask, &residualBudget[0], &agentCosts[0], totalReve);


      /*
      Rcout << "Knapsack solution, totalReve = " << totalReve << "\n";
      if(thereis) Rcout << "There is overloaded agent\n";
      else Rcout << "There is no overloaded agent\n";
      */


      if(!thereis)
      {
        if(totalReve > *currentSolutionRevenue)
        {
          mtx->lock();
          *currentSolutionRevenue = totalReve;
          currentSolution->assign(Bcontainer.begin(), Bcontainer.end());
          mtx->unlock();
        }
        bool bt = backtrack(T, B, Nagent, info, &residualBudget[0]);
        if(!bt) break;
      }
      else
      {
        pushAllBranchingVariableIntoStack<valtype, indtype, greedyBranch> (
            T, B, Nagent, overloadedAgent, stay, overloadedAgentTask,
            overloadedAgentWeight, overloadedAgentPenalty, &residualBudget[0], info);
      }
      /*
      {
      Rcout << "After KP locking, stack = \n";
      for(indtype i = 0, iend = T.size(); i < iend; ++i)
      {
      Rcout << T[i].agent << ", " << T[i].task << ", " << int(B[T[i].task][T[i].agent]) << "    ";
      }
      Rcout << "\n\n\n";
      // Rcout << "After KP locking, B = \n";
      // for(indtype i = 0; i < Nagent; ++i)
      // {
      //   for(indtype j = 0; j < Ntask; ++j)
      //   {
      //     if(B[j][i] >= 0) Rcout << " " << int(B[j][i]) << ", ";
      //     else Rcout << int(B[j][i]) << ", ";
      //   }
      //   Rcout << "\n";
      // }
      }
      */
    }


    // return currentSolutionRevenue;
    // return;
  }
};




template<typename valtype, typename indtype, typename fmove, bool greedyBranch>
struct runGapOBJ: public RcppParallel::Worker
{
  WV<valtype> **info;
  vec<signed char> *Bcontainers;
  std::time_t timer;
  double tlimit;
  tbb::spin_mutex *mtx;
  vec<signed char> *currentSolution;
  valtype *currentSolutionRevenue, *budget;
  indtype *nodes, *kps;
  gapOBJ<valtype, indtype, fmove, greedyBranch> *G;
  dynamicTasking *dT;
  void operator() (std::size_t st, std::size_t end)
  {
    for(;;)
    {
      INT objI = 0;
      if(!dT->nextTaskID(objI)) break;
      indtype Nnode = 0, Nkp = 0;
      G[st].run(info, Bcontainers[objI], timer, tlimit, mtx, currentSolution,
                currentSolutionRevenue, budget, Nnode, Nkp);
      nodes[st] += Nnode;
      kps[st] += Nkp;
    }
  }
  runGapOBJ(indtype Nagent, indtype Ntask,
            WV<valtype> **info, vec<vec<signed char> > &BcontainersV,
            std::time_t timer, double tlimit, vec<signed char> *currentSolution,
            valtype *currentSolutionRevenue,
            valtype *budget, indtype &totalNnodes, indtype &totalNkps,
            std::size_t maxCore): info(info), timer(timer), tlimit(tlimit),
            currentSolution(currentSolution), currentSolutionRevenue(currentSolutionRevenue),
            budget(budget) //, G(&Gvec)
  {
    dynamicTasking dt(maxCore, BcontainersV.size()); dT = &dt;
    tbb::spin_mutex m; mtx = &m;
    Bcontainers = &BcontainersV[0];
    vec<gapOBJ<valtype, indtype, fmove, greedyBranch> > Gvec(maxCore);
    for(int i = 0, iend = maxCore; i < iend; ++i)
    {
      Gvec[i].initialize(Nagent, Ntask);
    }
    G = &Gvec[0];
    vec<indtype> NnodesKps(maxCore * 2, 0);
    nodes = &NnodesKps[0];
    kps = &NnodesKps[maxCore];
      parallelFor(0, dT->NofCore, *this);
    totalNnodes = std::accumulate(nodes, nodes + maxCore, totalNnodes);
    totalNkps = std::accumulate(kps, kps + maxCore, totalNkps);
  }
};




// [[Rcpp::export]]
List auxGAPbbMulthreadNodes(
    NumericMatrix cost, NumericMatrix profitOrLoss, NumericVector budget,
    int maxCore = 7, int threadLoad = 32,
    double tlimit = 60, String ub = "MT",
    bool greedyBranching = true, String optim = "max")
{
  if(maxCore <= 1)
  {
    maxCore = 1;
    threadLoad = 1;
  }
  int Ntask = cost.ncol(), Nagent = cost.nrow();
  vec<WV<double> > costprofitInfo(cost.size());
  double C = -std::numeric_limits<double>::max();
  if(optim == "max")
  {
    for(int i = 0, iend = cost.size(); i < iend; ++i)
    {
      costprofitInfo[i].weight = cost[i];
      costprofitInfo[i].value = profitOrLoss[i];
    }
  }
  else
  {
    C = *std::max_element(profitOrLoss.begin(), profitOrLoss.end()) + 1;
    for(int i = 0, iend = cost.size(); i < iend; ++i)
    {
      costprofitInfo[i].weight = cost[i];
      costprofitInfo[i].value = C - profitOrLoss[i];
    }
  }


  vec<WV<double>*> wvptr(Ntask);
  for(int j = 0; j < Ntask; ++j)
  {
    wvptr[j] = &costprofitInfo[0] + j * INT(Nagent);
  }
  WV<double> **info = &wvptr[0];
  vec<double> residualBudget(budget.begin(), budget.end());
  vec<signed char> solution;


  double solutionRevenue = 0;
  int totalNnodes = 0, totalNkps = 0;
  std::time_t timer; time(&timer);


  // Stage I:
  vec<vec<signed char> > BcontainersV;
  {
    vec<stackEle<double, int> > T(INT(Nagent) * Ntask);
    T.resize(0);
    int stackDepth = std::log2(maxCore * threadLoad + 0.0) + 1;


    vec<signed char> Bcontainer(cost.size() + Ntask, -1);
    int stackDepthResv = stackDepth;
    solutionRevenue = gapBabStageI<double, int, MTfmoveUB<double, int>, false> (
      solution, Bcontainer, Nagent, Ntask, info, &residualBudget[0], timer,
      tlimit, totalNnodes, totalNkps, stackDepth, T);


    if(stackDepth != -1) // It means the function returns not because stack is full, but is actually ended.
    {
      goto wrapUp;
    }


    // Rcout << "stackDepth = " << stackDepth << "\n";
    // Rcout << "stackDepthResv = " << stackDepthResv << "\n";
    // for(int i = 0; i < (int)T.size(); ++i)
    // {
    //   Rcout << T[i].agent << ", " << T[i].task << ", " <<
    //     int(Bcontainer[T[i].task * (Nagent + 1) + T[i].agent]) << "\n";
    // }
    // return 0;


    unstack(Bcontainer, Nagent, BcontainersV, T, stackDepthResv);


    /*
    for(int i = 0, iend = BcontainersV.size(); i < iend; ++i)
    {
      for(int k = 0, kend = T.size(); k < kend; ++k)
      {
        int t = T[k].task, a = T[k].agent;
        Rcout << t << ", " << a << ", " << (int)BcontainersV[i][t * (Nagent + 1) + a] << "    ";
      }
      Rcout << "\n";
      for(int a = 0; a < Nagent + 1; ++a)
      {
        for(int t = 0; t < Ntask; ++t)
        {
          int tmp = BcontainersV[i][t * (Nagent + 1) + a];
          if(tmp >= 0) Rcout << " " << tmp << ", ";
          else Rcout << tmp << ", ";
        }
        Rcout << "\n";
      }
      Rcout << "\n\n";
    }
     */

    // return 0;


  }


  if(greedyBranching)
  {
    if(ub == "MT")
    {
      runGapOBJ<double, int, MTfmoveUB<double, int>, true> (
          Nagent, Ntask, info, BcontainersV, timer, tlimit, &solution,
          &solutionRevenue, &budget[0], totalNnodes, totalNkps, maxCore);
    }
    else
    {
      runGapOBJ<double, int, HSfmoveUB<double, int>, true> (
          Nagent, Ntask, info, BcontainersV, timer, tlimit,
          &solution, &solutionRevenue, &budget[0], totalNnodes, totalNkps, maxCore);
    }
  }
  else
  {
    if(ub == "MT")
    {
      runGapOBJ<double, int, MTfmoveUB<double, int>, false> (
          Nagent, Ntask, info, BcontainersV, timer, tlimit, &solution,
          &solutionRevenue, &budget[0], totalNnodes, totalNkps, maxCore);
    }
    else
    {
      runGapOBJ<double, int, HSfmoveUB<double, int>, false> (
          Nagent, Ntask, info, BcontainersV, timer, tlimit, &solution,
          &solutionRevenue, &budget[0], totalNnodes, totalNkps, maxCore);
    }
  }


  wrapUp:
  if(solutionRevenue == -std::numeric_limits<double>::max()) return List::create();
  if(C != -std::numeric_limits<double>::max())
  {
    solutionRevenue = C * Ntask - solutionRevenue;
  }


  NumericVector agentCost(Nagent);
  IntegerVector assignment(Ntask);
  for(int i = 0; i < Nagent; ++i)
  {
    agentCost[i] = 0;
    for(int j = 0; j < Ntask; ++j)
    {
      if(solution[j * (Nagent + 1) + i] <= 0) continue;
      agentCost[i] += cost[j * Nagent + i];
      assignment[j] = i + 1;
    }
  }
  return List::create(Named("totalProfitOrLoss") = solutionRevenue,
                      Named("agentCost") = agentCost,
                      Named("assignment") = assignment,
                      Named("nodes") = totalNnodes,
                      Named("bkpSolved") = totalNkps);
}
















































