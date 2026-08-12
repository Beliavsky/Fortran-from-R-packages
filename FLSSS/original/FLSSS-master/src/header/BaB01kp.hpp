# pragma once
// // [[Rcpp::depends(RcppParallel)]]
// # include <Rcpp.h>
// # include <RcppParallel.h>
# include "macros.hpp"
# include "dnyTasking.hpp"
// using namespace Rcpp;


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
inline void pushSequence(vec<indtype> &x, indtype jstart, indtype j) // j points to the critical item
{
  indtype &t = jstart;
  for(; t < j; ++t) x.push_back(t);
}


template<typename valtype, typename indtype>
inline bool backTrackOne(
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
// Elements in X have been ordered by values / weights.
// best and current will be resized in this function, but it is better to
// reserve the size of Xsize + 2 for each.
template<typename valtype, typename indtype, typename fmove, bool timeConstraint>
inline valtype bkp(
    kpEle<valtype, indtype> *X, indtype Xsize, valtype cap,
    indtype lenCap, vec<indtype> &best, vec<indtype> &current, double endtime)
{
  best.resize(0);
  current.resize(0);


  indtype j = -1; // This element has a huge weight, and is the first critical item.
  valtype bestVal = -std::numeric_limits<valtype>::max();
  valtype currentVal = 0;
  valtype residualCapacity = cap;
  indtype jstart = j + 1;


  INT iter = 0;
  while(true)
  {
    if(timeConstraint)
    {
      ++iter;
      if((iter & (512 - 1)) == 0)
      {
        if(double(std::clock()) > endtime) break;
      }
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




template<typename valtype, typename indtype>
struct cmp
{
  valtype *val;
  cmp(){}
  cmp(valtype *val): val(val){}
  bool operator() (indtype i, indtype j) {return val[i] > val[j];}
};



