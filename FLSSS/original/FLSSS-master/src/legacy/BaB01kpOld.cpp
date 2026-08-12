// # pragma once
// [[Rcpp::depends(RcppParallel)]]
# include <Rcpp.h>
# include <RcppParallel.h>
# include "header/BaB01kp.hpp"
// # include <fstream>
// # include "header/macros.hpp"
// # include "header/dnyTasking.hpp"
using namespace Rcpp;
// # define sINT std::int64_t;


/*
template<typename valtype, typename indtype>
struct kpEle
{
  // bool B;
  // valtype weight, value;
  valtype minWeightAfter;
  valtype accWeight, accValue;
  valtype valuePerWeight;
};


// HS algorithm: do a forward move and find the upper bound.
// Upon completion, j points to the one on the critical item's right.
// existingSum is the existing total value.
// Last argument ub is assigned by the upper bound found
template<typename valtype, typename indtype>
struct HSfmoveUB
{
  HSfmoveUB(valtype &existingSum, valtype &residualCapacity, valtype &ub,
            indtype &j, kpEle<valtype, indtype> *X, indtype residualNitem)
  {
    indtype jstart = j;
    valtype tmpCap = residualCapacity + X[jstart - 1].accWeight;
    indtype NitemCap = residualNitem + j;
    while(X[j].accWeight <= tmpCap and j < NitemCap) ++j;
    // The constraint j < Xsize is unncessary
    // because X[Xsize].accWeight being greater
    // than cap is guaranteed.
    existingSum += X[j - 1].accValue - X[jstart - 1].accValue;
    residualCapacity -= X[j - 1].accWeight - X[jstart - 1].accWeight;
    ub = existingSum + X[j].valuePerWeight * residualCapacity;
  }
};


// MT1R algorithm: do a forward move and return the upper bound.
// Upon completion, j points to the critical item.
template<typename valtype, typename indtype>
struct MTfmoveUB
{
  MTfmoveUB(valtype &existingSum, valtype &residualCapacity, valtype &ub,
            indtype &j, kpEle<valtype, indtype> *X, indtype residualNitem)
  {
    indtype jstart = j;
    valtype tmpCap = residualCapacity + X[jstart - 1].accWeight;
    indtype NitemCap = residualNitem + j;
    while(X[j].accWeight <= tmpCap and j < NitemCap) ++j;
    // The constraint j < Xsize is unncessary
    // because X[Xsize].accWeight being greater
    // than cap is guaranteed.
    existingSum += X[j - 1].accValue - X[jstart - 1].accValue;
    residualCapacity -= X[j - 1].accWeight - X[jstart - 1].accWeight;
    ub = existingSum + std::max<valtype> (
      X[j + 1].valuePerWeight * residualCapacity,
      X[j].accValue - X[j - 1].accValue -
        (X[j].accWeight - X[j - 1].accWeight - residualCapacity) *
        X[j - 1].valuePerWeight);
  }
  // existingSum does not include value in j.
  // residualCapacity does not include weight in j.
};


template<typename indtype>
void pushSequence(vec<indtype> &x, indtype jstart, indtype j) // j points to the critical item
{
  indtype &t = jstart;
  for(; t < j; ++t) x.push_back(t);
}


template<typename valtype, typename indtype>
bool backTrackOne(
    kpEle<valtype, indtype> *X, indtype &j,
    vec<indtype> &current, valtype &currentVal,
    valtype &residualCapacity)
{
  if(current.size() == 0) return false;
  j = current.back();
  current.pop_back();
  currentVal -= X[j].accValue - X[j - 1].accValue;
  residualCapacity += X[j].accWeight - X[j - 1].accWeight;
  return true;
}




// X[-1] and X[Xsize] have been defined
// best and current have sufficient capacities but zero sizes
template<typename valtype, typename indtype, typename fmove>
valtype bkp(kpEle<valtype, indtype> *X, indtype Xsize, valtype cap,
            indtype lenCap, vec<indtype> &best,
            vec<indtype> &current, double endtime)
{
  best.reserve(Xsize + 2); best.resize(0);
  current.reserve(Xsize + 2); current.resize(0);


  indtype j = -1; // This element has a huge weight, and is the first critical item.
  valtype bestVal = -std::numeric_limits<valtype>::max();
  valtype currentVal = 0;
  valtype residualCapacity = cap;
  indtype jstart = j + 1;


  INT iter = 0;
  while(true)
  {
    ++iter;
    if((iter & (512 - 1)) == 0)
    {
      if(double(std::clock()) > endtime) break;
    }


    if(X[jstart - 1].minWeightAfter <= residualCapacity)
    {
      while(X[jstart].accWeight - X[jstart - 1].accWeight
              > residualCapacity) ++jstart;
    }
    else
    {
      // j is assigned in backTrackOne().
      bool btsuccess = backTrackOne(X, j, current, currentVal, residualCapacity);
      if(!btsuccess) break;
      jstart = j + 1; // j is the critical item that should be excluded.
      continue;
    } // backtrack


    valtype tmpCurrentVal = currentVal;
    valtype tmpRedidualCap = residualCapacity;
    valtype ub;
    j = jstart;


    // Forward move.
    fmove(tmpCurrentVal, tmpRedidualCap, ub, j, X, lenCap - current.size());
    // After forward move, j points to the critical item.


    if(ub <= bestVal)
    {
      bool btsuccess = backTrackOne(X, j, current, currentVal, residualCapacity);
      if(!btsuccess) break;
      jstart = j + 1; // j is the critical item that should be excluded.
      continue;
    }


    // currentVal = tmpCurrentVal;
    // residualCapacity = tmpRedidualCap;
    pushSequence(current, jstart, j);
    if(tmpCurrentVal > bestVal)
    {
      best.assign(current.begin(), current.end());
      bestVal = tmpCurrentVal;
    }


    if(j >= Xsize)
    {
      current.resize(current.size() - (j - jstart));
      bool btsuccess = backTrackOne(X, j, current, currentVal, residualCapacity);
      if(!btsuccess) break;
    }
    else
    {
      currentVal = tmpCurrentVal;
      residualCapacity = tmpRedidualCap;
    }
    jstart = j + 1;
  }
  return bestVal;
}


struct cmp
{
  double *val;
  bool operator() (int i, int j) {return val[i] >= val[j];}
};
*/




template<typename valtype, typename indtype, typename fmove, bool timeConstraint>
struct paraBkpForCaps: public RcppParallel::Worker
{
  indtype Xsize;
  valtype *capV;
  indtype *lenCapV;
  double endtime;
  kpEle<valtype, indtype> *X;
  vec<indtype> *bestV, *currentV;
  valtype *bestVal;
  dynamicTasking *dT;
  void operator() (std::size_t st, std::size_t end)
  {
    for(;;)
    {
      INT objI = 0;
      if(!dT->nextTaskID(objI)) break;
      bestVal[objI] = bkp<valtype, indtype, fmove, timeConstraint> (
          X, Xsize, capV[objI], lenCapV[objI], bestV[objI], currentV[st], endtime);
    }
  }
  paraBkpForCaps(indtype Xsize, valtype *capV, indtype *lenCapV,
                 indtype capVsize, double endtime, kpEle<valtype, indtype> *X,
                 vec<vec<indtype> > &bestVec, // of size capVsize
                 valtype *bestVal, int maxCore): Xsize(Xsize),
                 capV(capV), lenCapV(lenCapV), endtime(endtime),
                 X(X), bestVal(bestVal)
  {
    dynamicTasking dt(maxCore, capVsize); dT = &dt;
    vec<vec<indtype> > current(maxCore, vec<indtype>(Xsize + 2));
    currentV = &current[0];
    bestV = &bestVec[0];
      parallelFor(0, dT->NofCore, *this);
  }
};




// List bkpOrdered(NumericVector weight, NumericVector value, double cap)
// [[Rcpp::export]]
List extra01knapsackBaB(
    NumericVector weight, NumericVector value, NumericVector caps,
    IntegerVector itemNcaps = IntegerVector(0), int maxCore = 7,
    double tlimit = 60, String ub = "MT")
{
  int Xsize = value.size();
  vec<kpEle<double, int> > Xcontain(Xsize + 2);
  kpEle<double, int> *X = &Xcontain[0] + 1;
  // Initialize Xcontain.
  double cap = *std::max_element(caps.begin(), caps.end());
  int *itemNcapV = nullptr;
  vec<int> itemNcapVcontainer;
  if(itemNcaps.size() == 0)
  {
    itemNcapVcontainer.assign(caps.size(), Xsize + 2);
    // Xsize + 2 is unncessary, Xsize is sufficient.
    itemNcapV = &itemNcapVcontainer[0];
  }
  else itemNcapV = &itemNcaps[0];


  vec<int> unitValOrder(Xsize);
  {
    X[-1].accValue = 0;
    // X[-1].accWeight = cap + 1;
    X[-1].accWeight = 0;
    vec<double> valuePerWeight(Xsize);
    for(int i = 0; i < Xsize; ++i)
    {
      unitValOrder[i] = i;
      valuePerWeight[i] = value[i] / weight[i];
    }
    cmp tmpCompare; tmpCompare.val = &valuePerWeight[0];
    std::sort(unitValOrder.begin(), unitValOrder.end(), tmpCompare);


    // Initialize X.accWeight, X.accValue, X.valuePerWeight.
    for(int i = 0; i < Xsize; ++i)
    {
      X[i].accWeight = X[i - 1].accWeight + weight[unitValOrder[i]];
      X[i].accValue = X[i - 1].accValue + value[unitValOrder[i]];
      X[i].valuePerWeight = valuePerWeight[unitValOrder[i]];
    }
    X[Xsize].accWeight = X[Xsize - 1].accWeight + (cap + 1);
    X[Xsize].accValue = X[Xsize - 1].accValue + 0;


    X[Xsize - 1].minWeightAfter = cap + 1;
    for(int i = Xsize - 2; i >= -1; --i)
    {
      X[i].minWeightAfter = std::min<double> (
        weight[unitValOrder[i + 1]], X[i + 1].minWeightAfter);
    }
  }


  maxCore = std::min<int> (maxCore, caps.size());
  vec<vec<int> > best(caps.size(), vec<int> (Xsize + 2));
  NumericVector bestVal(best.size());
  double endtime = (double)std::clock() + tlimit * CLOCKS_PER_SEC;


  if(ub == "MT")
    paraBkpForCaps<double, int, MTfmoveUB<double, int>, true> (
        Xsize, &caps[0], itemNcapV, caps.size(), endtime, X, best, &bestVal[0], maxCore);
  else
    paraBkpForCaps<double, int, HSfmoveUB<double, int>, true> (
        Xsize, &caps[0], itemNcapV, caps.size(), endtime, X, best, &bestVal[0], maxCore);


  List rst(caps.size());
  for(int i = 0, iend = rst.size(); i < iend; ++i)
  {
    IntegerVector v(best[i].size());
    for(int j = 0, jend = best[i].size(); j < jend; ++j)
    {
      v[j] = unitValOrder[best[i][j]] + 1;
    }
    rst[i] = v;
  }


  if(caps.size() > 1)
    return List::create(Named("maxVal") = bestVal, Named("selection") = rst);


  IntegerVector theSolution = rst[0];
  return List::create(Named("maxVal") = bestVal, Named("selection") = theSolution);
}





