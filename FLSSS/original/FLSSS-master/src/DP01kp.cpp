// # pragma once
// [[Rcpp::depends(RcppParallel)]]
# include <Rcpp.h>
# include <RcppParallel.h>
# include <thread>
# include "header/macros.hpp"
# include "header/dnyTasking.hpp"
# include "header/DP01kp.hpp"
using namespace Rcpp;
// # define sINT std::int64_t;


template<typename valtype, typename indtype, bool timeConstraint>
struct paraDp01kpForCaps: public RcppParallel::Worker
{
  indtype minCost_1;
  indtype Nitem, *caps;
  valtype **value;
  indtype *w;
  valtype *v;
  valtype *bestVal;
  vec<unsigned char> *selections;
  double endTime;
  dynamicTasking *dT;
  void operator() (std::size_t st, std::size_t end)
  {
    for(;;)
    {
      INT objI = 0;
      if(!dT->nextTaskID(objI)) break;
      dp01kp<valtype, indtype, timeConstraint> x;
      x.assign(minCost_1, w, value, v, endTime);
      bestVal[objI] = x.run(Nitem, caps[objI]);
      // Rcout << "recursion is finished\n";
      getSelection(value, Nitem, caps[objI], w, v, &selections[objI][0], minCost_1);
      // Rcout << "Selection is finished\n";
    }
  }
  paraDp01kpForCaps(indtype minCost_1, indtype Nitem, indtype *caps,
                    indtype capSize, valtype **value, indtype *w, valtype *v,
                    valtype *bestVal, vec<vec<unsigned char> > &selects,
                    double endTime, std::size_t maxCore):
    minCost_1(minCost_1), Nitem(Nitem), caps(caps), value(value),
    w(w), v(v), bestVal(bestVal), endTime(endTime)
  {
    selects.resize(capSize, vec<unsigned char> (Nitem, false));
    selections = &selects[0];
    dynamicTasking dt(maxCore, capSize); dT = &dt;
      parallelFor(0, dT->NofCore, *this);
  }
};




// [[Rcpp::export]]
List auxKnapsack01dp(IntegerVector weight, NumericVector value,
                       IntegerVector caps, int maxCore = 7,
                       double tlimit = 60, bool simplify = true)
{
  int minCost_1 = *std::min_element(weight.begin(), weight.end()) - 1;
  if(minCost_1 < 0) minCost_1 = 0;
  int cap = *std::max_element(caps.begin(), caps.end());
  NumericMatrix V(cap + 1 - minCost_1 , weight.size() + 1);
  std::fill(V.begin(), V.end(), -std::numeric_limits<double>::max());
  vec<double*> val(weight.size() + 1);
  for(int i = 0, iend = val.size(); i < iend; ++i)
  {
    val[i] = &V[0] + INT(i) * (cap + 1 - minCost_1);
  }


  vec<vec<unsigned char> > selection;
  NumericVector bestVal(caps.size());
  maxCore = std::min<int> (maxCore, caps.size());
  double endtime = (double)std::clock() + tlimit * CLOCKS_PER_SEC;
  paraDp01kpForCaps<double, int, true> (minCost_1,
    weight.size(), &caps[0], caps.size(), &val[0],
    &weight[0], &value[0], &bestVal[0], selection, endtime, maxCore);


  List solutions(caps.size());
  for(int i = 0, iend = caps.size(); i < iend; ++i)
  {
    IntegerVector rst(
        std::accumulate(selection[i].begin(), selection[i].end(), int(0)));
    for(int j = 0, k = 0, jend = selection[i].size(); j < jend; ++j)
    {
      if(!selection[i][j]) continue;
      rst[k] = j + 1;
      ++k;
    }
    solutions[i] = rst;
  }


  if(caps.size() == 1 and simplify)
  {
    IntegerVector theSelection = solutions[0];
    return List::create(Named("maxValue") = bestVal,
                        Named("selection") = theSelection,
                        Named("lookupTable") = V);
  }
  return List::create(Named("maxValue") = bestVal,
                      Named("selection") = solutions,
                      Named("lookupTable") = V);
}


















































