// [[Rcpp::depends(RcppParallel)]]
# include <Rcpp.h>
# include <RcppParallel.h>
// # include "header/BaB01kp.hpp"
# include "header/DP01kp.hpp"
using namespace Rcpp;


template<typename valtype, typename indtype>
struct WV
{
  indtype weight;
  valtype value;
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
bool best2agent(WV<valtype, indtype> *t, signed char *Bt, indtype Nagent,
                indtype &best, indtype &best2, indtype *residualBudget)
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
    WV<valtype, indtype> **info, indtype Nagent,
    indtype *residualBudget, indtype K)
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


/*
// Returning true means no need for backtrack.
template<typename valtype, typename indtype>
bool findOverloadedAgentsPenaltyWeightNextAgent(
    valtype &totalRevenue,
vec<indtype> &overloaded, WV<valtype, indtype> **info, signed char **B,
indtype Nagent, indtype Ntask, indtype *residualBudget,
indtype *budgetExceedance, indtype *agentCost,
vec<vec<indtype> > &agentTask,
vec<vec<indtype> > &agentWeight, // of size Nagent
vec<vec<valtype> > &agentPenalty, // will be used as values in knapsacking
vec<vec<indtype> > &nextAgent)
{
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
indtype i = 0;
while(B[j][i] != 2) ++i;
totalRevenue += info[j][i].value;
continue;
}
indtype best, best2;
bool needNoBacktrack = best2agent( // Also set nonzero variables to -1.
  info[j], B[j], Nagent, best, best2, residualBudget);
if(!needNoBacktrack) return false; // need back track.
B[j][best] = 1;
totalRevenue += info[j][best].value;
valtype w = info[j][best].weight;
agentWeight[best].push_back(w);
agentCost[best] += w;
nextAgent[best].push_back(best2);
agentPenalty[best].push_back(info[j][best].value - info[j][best2].value);
agentTask[best].push_back(j);
}


overloaded.resize(0);
for(indtype i = 0; i < Nagent; ++i)
{
budgetExceedance[i] = agentCost[i] - residualBudget[i];
if(budgetExceedance[i] > 0) overloaded.push_back(i);
}
return true;
}
*/


// Returning true means no need for backtrack.
template<typename valtype, typename indtype>
bool findOverloadedAgentsPenaltyWeightNextAgent(
    valtype &totalRevenue,
    vec<indtype> &overloaded, WV<valtype, indtype> **info,
    signed char **B, indtype Nagent, indtype Ntask,
    indtype *residualBudget, indtype *budgetExceedance, indtype *agentCost,
    vec<vec<indtype> > &agentTask,
    vec<vec<indtype> > &agentWeight, // of size Nagent
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
bool thereIsOverlodedAgent(
    WV<valtype, indtype> **info,signed char **B, indtype Nagent, indtype Ntask,
    indtype *residualBudget, indtype *agentCost, valtype &totalRevenue)
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
    vec<vec<indtype> > &weight, // overloaded agent, task weights
    vec<vec<valtype> > &penalty, // overloaded agent, task values
    indtype *residualBudget, WV<valtype, indtype> **info)
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
    indtype *w = &weight[a][0];
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
               indtype Nagent, WV<valtype, indtype> **info,
               indtype *residualBudget)
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




template<typename valtype, typename indtype>
struct KPinGap: public dp01kp<valtype, indtype, false>
{
  indtype Nitem, cap;
  vec<valtype> content;
  void assign(indtype givenNitem, indtype givenCap,
              indtype minCost_1, indtype *w, valtype *v)
  {
    Nitem = givenNitem;
    cap = std::accumulate(w, w + Nitem, 0) - givenCap;
    INT nrow = cap + 1 - minCost_1, ncol = Nitem + 1;
    // The first ncol slots in content store pointers.
    content.assign(INT(nrow) * ncol + ncol, -std::numeric_limits<valtype>::max());
    dp01kp<valtype, indtype, false>::assign(minCost_1, w, (valtype**)(&content[0]), v, 0);
    valtype *val = &content[ncol];
    for(INT i = 0; i < ncol; ++i)
    {
      dp01kp<valtype, indtype, false>::value[i] = val + nrow * i;
    }
  }


  // There is an run() function in the base class.
  valtype run(vec<indtype> &stay, vec<indtype> &reassign)
  {
    valtype S = dp01kp<valtype, indtype, false>::run(Nitem, cap); // defined in struct dp01kp{}
    indtype i = Nitem, j = cap;
    stay.resize(0);
    reassign.resize(0);
    for(; i > 0; --i)
    {
      if(j >= dp01kp<valtype, indtype, false>::minCost_1 and
           dp01kp<valtype, indtype, false>::value[i][j - dp01kp<valtype, indtype, false>::minCost_1] !=
             dp01kp<valtype, indtype, false>::value[i - 1][j - dp01kp<valtype, indtype, false>::minCost_1] and
           j >= dp01kp<valtype, indtype, false>::w[i - 1])
      {
        stay.push_back(i - 1);
        j -= dp01kp<valtype, indtype, false>::w[i - 1];
      }
      else reassign.push_back(i - 1);
    }
    return S; // as a result, stay and reassign are in reverse order.
  }
};




// return totalPenalty
template<typename valtype, typename indtype>
valtype specialBiKpDPmulti(
    KPinGap<valtype, indtype> &kpg,
    vec<indtype> &overloadedAgent, vec<vec<indtype> > &weight, vec<vec<valtype> > &penalty,
    vec<vec<indtype> > &stay, vec<vec<indtype> > &reassign, vec<indtype> &budgetExceedance)
{
  valtype totalPenalty = 0;
  for(indtype objI = 0, end = overloadedAgent.size(); objI < end; ++objI)
  {
    indtype a = overloadedAgent[objI];
    kpg.assign(
        weight[a].size(), budgetExceedance[a], *std::min_element(weight[a].begin(), weight[a].end()) - 1,
        &weight[a][0], &penalty[a][0]); // re-allocation will not happen every time it runs here.
    totalPenalty += kpg.run(stay[objI], reassign[objI]);
    // Print stay and reassgin
    // {
    //   Rcout << "stay:\n";
    //   for(indtype i = 0, iend = stay[objI].size(); i < iend; ++i)
    //   {
    //     Rcout << stay[objI][i] << ", ";
    //   }
    //   Rcout << "\n";
    //   Rcout << "reassign:\n";
    //   for(indtype i = 0, iend = reassign[objI].size(); i < iend; ++i)
    //   {
    //     Rcout << reassign[objI][i] << ", ";
    //   }
    //   Rcout << "\n\n\n";
    // }
  }
  // Rcout << "-------------------------\n";
  return totalPenalty;
}


// B has been initialized. Elements are either -1 or 1.
template<typename valtype, typename indtype, bool greedyBranch>
valtype gapBabDpStageI(vec<signed char> &currentSolution, vec<signed char> &Bcontainer,
                 indtype Nagent, indtype Ntask, WV<valtype, indtype> **info,
                 indtype *residualBudget, std::time_t timer, double tlimit,
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
  vec<vec<indtype> > overloadedAgentWeight(Nagent, vec<indtype>(Ntask));
  vec<vec<valtype> > overloadedAgentPenalty(Nagent, vec<valtype>(Ntask)); // will be used as values in knapsacking
  vec<vec<indtype> > nextAgent(Nagent, vec<indtype>(Ntask));
  vec<vec<indtype> > reassign(Nagent, vec<indtype>(Ntask));
  vec<vec<indtype> > stay(Nagent, vec<indtype>(Ntask));


  vec<indtype> budgetExceedance(Nagent);
  T.reserve(INT(Nagent) * Ntask);
  currentSolution.resize(INT(Nagent) * Ntask);
  valtype currentSolutionRevenue = -std::numeric_limits<valtype>::max();


  // Auxiliary containers.
  vec<indtype> agentCosts(Nagent);


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
        Nagent, Ntask, residualBudget, &budgetExceedance[0], &agentCosts[0], overloadedAgentTask,
        overloadedAgentWeight, // of size Nagent
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
      }
      */
      if(!bt) return currentSolutionRevenue;
    }


    /*
    {
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
    }
    */


    // valtype totalPenalty = 0;
    ++Nnode; Nkp += overloadedAgent.size();
    KPinGap<valtype, indtype> kpg; // aux;
    specialBiKpDPmulti<valtype, indtype> (kpg,
        overloadedAgent, overloadedAgentWeight,
        overloadedAgentPenalty, stay, reassign, budgetExceedance);


    // Set variables pointed by reassigned to -1.
    updateBafterKnapsacking<indtype> (overloadedAgent, B, nextAgent, overloadedAgentTask, reassign);


    /*
    {
    Rcout << "After KP, reassign = \n";
    for(indtype i = 0, iend = overloadedAgent.size(); i < iend; ++i)
    {
    for(indtype j = 0, jend = reassign[i].size(); j < jend; ++j)
    {
    Rcout << reassign[i][j] << ", ";
    }
    Rcout << "\n";
    }
    Rcout << "After KP, stay = \n";
    for(indtype i = 0, iend = overloadedAgent.size(); i < iend; ++i)
    {
    for(indtype j = 0, jend = stay[i].size(); j < jend; ++j)
    {
    Rcout << stay[i][j] << ", ";
    }
    Rcout << "\n";
    }
    Rcout << "After KP, B = \n";
    for(indtype i = 0; i < Nagent; ++i)
    {
    for(indtype j = 0; j < Ntask; ++j)
    {
    if(B[j][i] >= 0) Rcout << " " << int(B[j][i]) << ", ";
    else Rcout << int(B[j][i]) << ", ";
    }
    Rcout << "\n";
    }
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
     Rcout << "After KP locking, B = \n";
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
template<typename valtype, typename indtype, bool greedyBranch>
struct gapOBJ
{
  // vec<signed char> *Bctn;
  vec<indtype> residualBudget;
  vec<signed char*> Bv;
  vec<indtype> overloadedAgent;
  vec<vec<indtype> > overloadedAgentTask;
  vec<vec<indtype> > overloadedAgentWeight;
  vec<vec<valtype> > overloadedAgentPenalty; // will be used as values in knapsacking
  vec<vec<indtype> > nextAgent;
  vec<vec<indtype> > reassign;
  vec<vec<indtype> > stay;
  vec<indtype> budgetExceedance;
  vec<stackEle<valtype, indtype> > T;
  KPinGap<valtype, indtype> kpg; // aux
  vec<indtype> agentCosts; // aux


  void initialize(indtype Nagent, indtype Ntask)
  {
    Bv.resize(Ntask);
    residualBudget.resize(Nagent); // has not been calculated.
    overloadedAgent.reserve(Nagent);
    overloadedAgentTask.resize(Nagent, vec<indtype>(Ntask));
    overloadedAgentWeight.resize(Nagent, vec<indtype>(Ntask));
    overloadedAgentPenalty.resize(Nagent, vec<valtype>(Ntask));
    nextAgent.resize(Nagent, vec<indtype>(Ntask));
    reassign.resize(Nagent, vec<indtype>(Ntask));
    stay.resize(Nagent, vec<indtype>(Ntask));
    budgetExceedance.resize(Nagent);
    T.reserve(INT(Nagent) * Ntask);
    agentCosts.resize(Nagent); // aux
  }


  // Return false if time is up.
  void run(WV<valtype, indtype> **info, vec<signed char> &Bctn,
           std::time_t timer, double tlimit, tbb::spin_mutex *mtx,
           vec<signed char> *currentSolution,
           valtype *currentSolutionRevenue, indtype *budget,
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


    while(true)
    {
      std::time_t now; std::time(&now);
      if(std::difftime(now, timer) > tlimit) break;


      while(true) // Repeat backtracking until knapsacking becomes necessary.
      {
        valtype revenueUB = 0;
        bool needNoBacktrack = findOverloadedAgentsPenaltyWeightNextAgent(
          revenueUB, overloadedAgent, info, B, Nagent, Ntask,
          &residualBudget[0], &budgetExceedance[0], &agentCosts[0], overloadedAgentTask,
          overloadedAgentWeight, // of size Nagent
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
      Rcout << "Current max profit = " << *currentSolutionRevenue << "\n";
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
      specialBiKpDPmulti<valtype, indtype> (kpg,
          overloadedAgent, overloadedAgentWeight, overloadedAgentPenalty,
          stay, reassign, budgetExceedance);


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
  }
};


template<typename valtype, typename indtype, bool greedyBranch>
struct runGapOBJ: public RcppParallel::Worker
{
  WV<valtype, indtype> **info;
  vec<signed char> *Bcontainers;
  std::time_t timer;
  double tlimit;
  tbb::spin_mutex *mtx;
  vec<signed char> *currentSolution;
  valtype *currentSolutionRevenue;
  indtype *budget;
  indtype *nodes, *kps;
  gapOBJ<valtype, indtype, greedyBranch> *G;
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
            WV<valtype, indtype> **info, vec<vec<signed char> > &BcontainersV,
            std::time_t timer, double tlimit, vec<signed char> *currentSolution,
            valtype *currentSolutionRevenue,
            indtype *budget, indtype &totalNnodes, indtype &totalNkps,
            std::size_t maxCore): info(info), timer(timer), tlimit(tlimit),
            currentSolution(currentSolution), currentSolutionRevenue(currentSolutionRevenue),
            budget(budget) //, G(&Gvec)
  {
    dynamicTasking dt(maxCore, BcontainersV.size()); dT = &dt;
    tbb::spin_mutex m; mtx = &m;
    Bcontainers = &BcontainersV[0];
    vec<gapOBJ<valtype, indtype, greedyBranch> > Gvec(maxCore);
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
List auxGAPbbDpMulthreadNodes(IntegerMatrix cost, NumericMatrix profitOrLoss, IntegerVector budget,
                            int maxCore = 7, int threadLoad = 32, double tlimit = 60,
                            bool greedyBranching = true, String optim = "max")
{
  if(maxCore <= 1)
  {
    maxCore = 1;
    threadLoad = 1;
  }
  int Ntask = cost.ncol(), Nagent = cost.nrow();
  vec<signed char> Bcontainer(cost.size() + Ntask, -1);
  vec<WV<double, int> > costprofitInfo(cost.size());
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
  vec<WV<double, int>*> wvptr(Ntask);
  for(int j = 0; j < Ntask; ++j)
  {
    wvptr[j] = &costprofitInfo[0] + j * INT(Nagent);
  }
  WV<double, int> **info = &wvptr[0];
  vec<int> residualBudget(budget.begin(), budget.end());
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
    solutionRevenue = gapBabDpStageI<double, int, false> (
      solution, Bcontainer, Nagent, Ntask, info, &residualBudget[0], timer,
      tlimit, totalNnodes, totalNkps, stackDepth, T);


    if(stackDepth != -1) // It means the function returns not because stack is full, but is actually ended.
    {
      goto wrapUp;
    }
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
  }



  if(greedyBranching)
  {
    runGapOBJ<double, int, true> (
        Nagent, Ntask, info, BcontainersV, timer, tlimit, &solution,
        &solutionRevenue, &budget[0], totalNnodes, totalNkps, maxCore);
  }
  else
  {
    runGapOBJ<double, int, false> (
        Nagent, Ntask, info, BcontainersV, timer, tlimit, &solution,
        &solutionRevenue, &budget[0], totalNnodes, totalNkps, maxCore);
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





























