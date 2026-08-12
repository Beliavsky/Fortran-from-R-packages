// # pragma once
// [[Rcpp::depends(RcppParallel)]]
# include <Rcpp.h>
# include <RcppParallel.h>
# include <mutex>
# include <fstream>
// # include <iostream>
// # include <ctime>
// # include <setjmp.h>
// # include "../header/oneDoperation.hpp"
// # include "../header/mvalOperation.hpp"
# include "header/macros.hpp"
# include "header/dnyTasking.hpp"
using namespace Rcpp;


template<typename valtype, typename indtype>
struct task // there will be tasks many boxes
{
  indtype *ird; // irregular dimensions
  // ird[i] is the answer to "which dimension is irregular in the i_th row ?"
  valtype *val; // values of irregular dimensions
  valtype *profit;
};


template<typename valtype, typename indtype>
void getV(NumericMatrix dividedV, vec<valtype> &container,
          NumericVector profitV, vec<task<valtype, indtype> > &V)
{
  indtype nagent = dividedV.nrow();
  indtype ntask = dividedV.ncol() / nagent;
  double tmp = (sizeof(indtype) * nagent) / (sizeof(valtype) + 0.0);
  int indexVecSize = tmp;
  if(tmp > indexVecSize) ++indexVecSize;
  container.resize((indexVecSize + nagent + nagent) * ntask);
  V.resize(ntask);
  valtype *x = &container[0];
  for(indtype i = 0; i < ntask; ++i)
  {
    valtype *xi = x + unsigned(indexVecSize + nagent + nagent) * i;
    V[i].ird = (indtype*)xi;
    V[i].val = xi + indexVecSize;
    V[i].profit = V[i].val + nagent;
  }


  double *st = &dividedV[0];
  for(int i = 0; i < ntask; ++i)
  {
    double *x = st + i * nagent * nagent;
    for(int j = 0; j < nagent; ++j) // j_th row in box
    {
      double *y = x + j * nagent;
      for(int k = 0; k < nagent; ++k)
      {
        if(y[k] - j > eps)
        {
          V[i].ird[j] = k;
          V[i].val[j] = y[k];
          break;
        }
      }
    }
  }
  int k = 0;
  for(int i = 0; i < ntask; ++i)
  {
    for(int j = 0; j < nagent; ++j)
    {
      V[i].profit[j] = profitV[k];
      ++k;
    }
  }
}


// Do not delete! Code valuable!
/*
template<typename valtype, typename indtype>
indtype findBound002(indtype nagent, indtype ntask,
                  // vec<valtype> &acntr, // a temp container of size (nagent + 1) * 2
                  task<valtype, indtype> *T,
                  indtype *LB,
                  indtype *UB,
                  valtype *MIN_sumUB, indtype &MIN_sumUB_maxDim,
                  valtype *MAX_sumLB, indtype &MAX_sumLB_minDim
                  )
{
  // sumLB, sumUB, MIN and MAX are of size nagent + 1. The last dimension is the key index
  bool boo = false;
  int LBsum = 0, UBsum = 0, LBUBdiff = 0;


  while(true)
  {
    bool boundChanged = false;
    indtype I = 0;
    LBsum = 0;
    LBUBdiff = 0;
    indtype sumDiff = 0;


    for(; I < ntask; ++I)
    {
      valtype *&val = T[I].val;
      indtype *&ird = T[I].ird;


      // Rcout << "UB[I] = " << UB[I] << "\n";
      // {
      //   for(int i = 0; i <= nagent; ++i)
      //   {
      //     Rcout << MIN_sumUB[i] << ", ";
      //   }
      //   Rcout << "\n";
      // }
      valtype max = MIN_sumUB[MIN_sumUB_maxDim] + UB[I]; // Would-be
      indtype &irrDim = ird[UB[I]];
      valtype tmp = MIN_sumUB[irrDim] + val[UB[I]]; // Would-be
      if(max < tmp) // No need to do max < tmp - eps
      {
        max = tmp;
        // MIN_sumUB_maxDim = irrDim;
      }
      tmp = indtype(max);
      if(std::abs(max - tmp) < eps) // cautious about numeric issues
      {
        max = tmp;
      }
      // Rcout << "max = " << max << ", LB[I] = " << LB[I] << "\n";
      if(max < LB[I] - eps)
      {
        LBUBdiff += UB[I] - LB[I];
        continue;
      }
      indtype LBI = LB[I];
      LB[I] = max;
      // Rcout << "LB[I] = " << LB[I] << "\n";
      if(LB[I] >= nagent) return 0;
      // tmp = val[ird[LB[I]]];
      tmp = val[LB[I]]; // val[LB[I]] is that irregular value in row LB[I]
      // Rcout << "ird[LB[I]] = " << ird[LB[I]] << "\n";
      if(tmp < max - eps) // if max is integer, 'LB[I] = max' has already covered that case.
      {
        // Rcout << "tmp < max - eps and tmp = " << tmp << "\n";
        ++LB[I];
      }
      if(LB[I] > UB[I]) return 0;
      LBUBdiff += UB[I] - LB[I];
      // update MAX_sumLB, with the extra portion
      indtype diff = LB[I] - LBI;
      if(diff != 0)
      {
        MAX_sumLB[ird[LBI]] += val[LBI] - LBI; // add back the decimals of the special dimension, old LBI
        MAX_sumLB[ird[LB[I]]] -= val[LB[I]] - LB[I]; // subtract the decimals of the special dimension, new LB[I]
        boundChanged = true;
      }
      sumDiff += diff;
      LBsum += LB[I];
    }


    // {
    //   Rcout << "LB = ";
    //   for(int i = 0; i < ntask; ++i)
    //   {
    //     Rcout << LB[i] << ", ";
    //   }
    //   Rcout << "\n";
    // }


    // Rcout << nagent << "\n";
    // Update MAX_sumLB and MAX_sumLB_minDim
    {
      MAX_sumLB_minDim = 0;
      valtype currentMin = infinity;
      for(indtype i = 0; i <= nagent; ++i)
      {
        MAX_sumLB[i] -= sumDiff;
        if(MAX_sumLB[i] < currentMin)
        {
          currentMin = MAX_sumLB[i];
          MAX_sumLB_minDim = i;
        }
      }
    }


    // Can I stop now
    {
      if(!boo) boo = true;
      else
      {
        if(!boundChanged)
        {
          if(LBUBdiff == 0)
          {
            // Rcout << "LBUBdiff = " << LBUBdiff << "\n";
            return 2;
          }
          break;
        }
      }
    }


    boundChanged = false;
    sumDiff = 0;
    UBsum = 0;
    LBUBdiff = 0;
    I = 0;
    for(; I < ntask; ++I)
    {
      valtype *&val = T[I].val;
      indtype *&ird = T[I].ird;


      valtype min = MAX_sumLB[MAX_sumLB_minDim] + LB[I];
      indtype &irrDim = ird[LB[I]];
      valtype tmp = MAX_sumLB[irrDim] + val[LB[I]];
      if(tmp < min) // update MAX_sumLB_minDim if possible
      {
        min = tmp;
        // MAX_sumLB_minDim = irrDim;
      }
      tmp = indtype(min);  // cautious about numeric issues
      if(std::abs(tmp - min) < eps) min = tmp;
      if(min > UB[I] + eps)
      {
        LBUBdiff += UB[I] - LB[I];
        continue;
      }
      indtype UBI = UB[I];
      if(min < 0 - eps) return 0;
      UB[I] = min;
      tmp = val[UB[I]];
      if(tmp > min + eps) --UB[I];
      if(LB[I] > UB[I]) return 0;
      LBUBdiff += UB[I] - LB[I];
      // update MAX_sumLB, with the extra portion
      indtype diff = UB[I] - UBI;
      if(diff != 0)
      {
        MIN_sumUB[ird[UBI]] += val[UBI] - UBI; // add back the decimals of the special dimension, old LBI
        MIN_sumUB[ird[UB[I]]] -= val[UB[I]] - UB[I]; // subtract the decimals of the special dimension, new LB[I]
        boundChanged = true;
      }
      sumDiff += diff;
      UBsum += UB[I];
    }


    // Update MAX_sumLB and MAX_sumLB_minDim
    {
      MIN_sumUB_maxDim = 0;
      valtype currentMax = -infinity;
      for(indtype i = 0; i <= nagent; ++i)
      {
        MIN_sumUB[i] -= sumDiff;
        if(MIN_sumUB[i] > currentMax)
        {
          currentMax = MIN_sumUB[i];
          MIN_sumUB_maxDim = i;
        }
      }
    }


    // can I jump out now
    {
      if(!boundChanged)
      {
        if(LBUBdiff == 0)
        {
          return 2;
        }
        break;
      }
    }
  }


  return 1;
}


// [[Rcpp::export]]
List testFindBound002GAP(NumericMatrix dividedV, NumericVector target, NumericVector ME)
{
  // target and ME's dimensionality is nagent + 1
  int nagent = dividedV.nrow();
  int ntask = dividedV.ncol() / nagent;
  vec<double> Vcontainer;
  vec<task<double, int> > V;
  getV<double, int> (dividedV, Vcontainer, V);
  // vec<double> acntr((nagent + 1) * 2);
  task<double, int> *t = &V[0];
  vec<int> LB(ntask, 0), UB(ntask, nagent - 1);
  vec<double> sumLB(nagent + 1, 0), sumUB(nagent + 1, (nagent - 1) * ntask);
  // vec<double> MIN_sumUB(nagent), MAX_sumLB(nagent);
  for(int i = 0; i < ntask; ++i)
  {
    int *&ird = t[i].ird;
    double *&val = t[i].val;
    // Rcout << "=========\n";
    // for(int k = 0; k < nagent; ++k)
    // {
    //   Rcout << ird[k] << ", ";
    //   Rcout << val[k] << "; ";
    // }
    // Rcout << "=========\n";
    int &irdim = ird[LB[i]];
    // Rcout << "irdim = " << irdim << ", " << val[irdim] << ", ";
    sumLB[irdim] += val[LB[i]] - LB[i];
    // Rcout << sumLB[irdim] << "\n";
  }
  // Rcout << "----\n";
  for(int i = 0; i < ntask; ++i)
  {
    int *&ird = t[i].ird;
    // Rcout << "=========\n";
    // for(int k = 0; k < nagent; ++k)
    // {
    //   Rcout << ird[k] << ", ";
    // }
    // Rcout << "=========\n";
    double *&val = t[i].val;
    int &irdim = ird[UB[i]];
    // std::cout << "UB[i] = " << UB[i] << "\n"; // this one checked, correct
    // Rcout << "irdim = " << irdim << ", " << val[irdim] << ", ";
    sumUB[irdim] += val[UB[i]] - UB[i];
    // Rcout << sumUB[irdim] << "\n";
  }


  // Rcout << "sumUB = ";
  // for(int i = 0; i <= nagent; ++i)
  // {
  //   Rcout << sumUB[i] << ", ";
  // }
  // Rcout << "end of sumUB\n";


  vec<double> MIN_sumUB(target.size()), MAX_sumLB(target.size());
  for(int i = 0, iend = MIN_sumUB.size(); i < iend; ++i)
  {
    MIN_sumUB[i] = target[i] - ME[i] - sumUB[i];
    MAX_sumLB[i] = target[i] + ME[i] - sumLB[i];
  }


  // {
  //   for(int i = 0; i <= nagent; ++i)
  //   {
  //     Rcout << MIN_sumUB[i] << " ";
  //   }
  //   Rcout << "\n";
  //   for(int i = 0; i <= nagent; ++i)
  //   {
  //     Rcout << MAX_sumLB[i] << " ";
  //   }
  //   Rcout << "\n";
  // }


  int MIN_sumUB_maxDim = std::max_element(MIN_sumUB.begin(), MIN_sumUB.end()) - MIN_sumUB.begin();
  int MAX_sumLB_minDim = std::min_element(MAX_sumLB.begin(), MAX_sumLB.end()) - MAX_sumLB.begin();


  // Rcout << MIN_sumUB_maxDim << ", " << MAX_sumLB_minDim << "\n";


  int boo = findBound002<double, int> (
    nagent, ntask, &V[0], &LB[0], &UB[0],
    &MIN_sumUB[0], MIN_sumUB_maxDim, &MAX_sumLB[0], MAX_sumLB_minDim);


  return List::create(Named("LB") = LB, Named("UB") = UB, Named("boo") = boo);
}
*/


// The dimensionality of LB or UB equals nagent; the dimensionality of MAX_sumLB is of nagent + 1 ---
// the extra dimension is for the index column.
template<typename valtype, typename indtype>
indtype findBound003(indtype nagent, indtype ntask, // nagent is the dimension and never changes
                     task<valtype, indtype> *T, indtype *taskInd,
                     indtype *LB, indtype *UB,
                     indtype &MIN_sumUBindVal,
                     // Lower bound is not restricted. All we need to know is the information of the index column.
                     // valtype *MIN_sumUB, indtype &MIN_sumUB_maxDim,
                     valtype *MAX_sumLB, indtype &MAX_sumLB_minDim,
                     indtype &MAX_sumLB_2ndMinDim
)
{
  // sumLB, sumUB, MIN and MAX are of size nagent + 1. The last dimension is the key index
  bool boo = false;
  while(true)
  {
    bool boundChanged = false;
    indtype I = 0;
    int LBUBdiff = 0;
    int sumDiff = 0;


    // Rcout << "MIN_sumUBindVal = " << (int)MIN_sumUBindVal << "\n";


    for(; I < ntask; ++I)
    {
      indtype LBI = MIN_sumUBindVal + UB[I];
      if(LBI < LB[I])
      {
        LBUBdiff += UB[I] - LB[I];
        continue;
      }
      if(LBI > UB[I]) return 0;
      std::swap(LBI, LB[I]); // LBI now equals the old LB[I]
      boundChanged = boundChanged or LBI != LB[I];


      LBUBdiff += UB[I] - LB[I];
      valtype *&val = T[taskInd[I]].val;
      indtype *&ird = T[taskInd[I]].ird;
      MAX_sumLB[ird[LBI]] += val[LBI] - LBI; // add back the decimals of the special dimension, old LBI
      MAX_sumLB[ird[LB[I]]] -= val[LB[I]] - LB[I]; // subtract the decimals of the special dimension, new LB[I]
      sumDiff += LB[I] - LBI;
    }


    // Update MAX_sumLB, MAX_sumLB_minDim, MAX_sumLB_2ndMinDim
    {
      for(indtype i = 0; i <= nagent; ++i)
      {
        MAX_sumLB[i] -= sumDiff;
      }


      MAX_sumLB_minDim = 0;
      MAX_sumLB_2ndMinDim = 1;
      if(MAX_sumLB[MAX_sumLB_minDim] > MAX_sumLB[MAX_sumLB_2ndMinDim])
        std::swap(MAX_sumLB_minDim, MAX_sumLB_2ndMinDim);
      for(indtype i = 2; i <= nagent; ++i)
      {
        if(MAX_sumLB[i] < MAX_sumLB[MAX_sumLB_minDim])
        {
          indtype tmp = MAX_sumLB_minDim;
          MAX_sumLB_minDim = i;
          MAX_sumLB_2ndMinDim = tmp;
        }
        else if(MAX_sumLB[i] < MAX_sumLB[MAX_sumLB_2ndMinDim])
        {
          MAX_sumLB_2ndMinDim = i;
        }
      }
    }


    // Can I stop now
    if(!boo) boo = true;
    else
    {
      if(!boundChanged)
      {
        if(LBUBdiff == 0) return 2;
        break;
      }
    }


    boundChanged = false;
    LBUBdiff = 0;
    I = 0;


    for(; I < ntask; ++I)
    {
      valtype *&val = T[taskInd[I]].val;
      indtype *&ird = T[taskInd[I]].ird;


      indtype LBirrDim = ird[LB[I]];


      // There are 2 dimensions and 3 values that matter:
      // valtype MAX_sumLB_minDim_val = 0;
      // valtype MAX_sumLB_irrDim_val = 0;
      valtype min = 0;


      if(LBirrDim == MAX_sumLB_minDim)
      {
        // Rcout << "irrDim == MAX_sumLB_minDim, and irrDim = " << irrDim << "\n";
        valtype MAX_sumLB_minDim_val = MAX_sumLB[LBirrDim] + val[LB[I]];
        min = MAX_sumLB[MAX_sumLB_2ndMinDim] + LB[I];
        min = std::min<valtype> (min, MAX_sumLB_minDim_val);
      }
      else
      {
        // Rcout << "irrDim != MAX_sumLB_minDim, and irrDim = " << irrDim << "\n";
        // valtype MAX_sumLB_minDim_val = MAX_sumLB[MAX_sumLB_minDim] + LB[I];
        // min = MAX_sumLB_minDim_val;
        min = MAX_sumLB[MAX_sumLB_minDim] + LB[I];
      }


      // irrDim = ird[UB[I]];
      valtype tmp = indtype(min); // Cautious about numeric issues
      if(std::abs(tmp - min) < eps) min = tmp;
      if(min < 0 - eps)
      {
        // Rcout << "min < 0 - eps)\n";
        return 0;
      }


      indtype UBI = min;
      if(UBI > UB[I]) // When the new UBI is not as good as the old UB[I]
      {
        LBUBdiff += UB[I] - LB[I];
        continue;
      }


      // if(UBI < 0) return 0;
      // Rcout << "MAX_sumLB[ird[UBI]] = " << MAX_sumLB[ird[UBI]] << ", ";
      // Rcout << "val[UBI] - eps = " << val[UBI] - eps << "\n";
      if(UBI == LB[I]) // equivalently, LBirrDim == ird[UBI]
      {
        // Do not delete this line :// if(MAX_sumLB[ird[UBI]] + val[LB[I]] < val[UBI] - eps) --UBI; return 0;
        if(MAX_sumLB[ird[UBI]] < 0 - eps) --UBI;
      }
      else
      {
        if(MAX_sumLB[ird[UBI]] + LB[I] < val[UBI] - eps) --UBI;
      }


      if(UBI < LB[I])
      {
        // Rcout << "I == " << I << ", ";
        // Rcout << "UBI = " << UBI << ", UB[I] = " << UB[I] << "\n";
        // Rcout << "UBI < LB[I]\n";
        return 0;
      }
      std::swap(UBI, UB[I]);
      boundChanged = boundChanged or UBI != UB[I];


      LBUBdiff += UB[I] - LB[I]; // Difference between the lower and upper bounds.
      // Rcout << "UB[I] = " << UB[I] << ", UBI = " << UBI << "\n";
      MIN_sumUBindVal += UBI - UB[I];
    }


    // Can I jump out now
    if(!boundChanged)
    {
      if(LBUBdiff == 0) return 2;
      break;
    }
  }


  return 1;
}




// [[Rcpp::export]]
List testFindBound003GAP(NumericMatrix dividedV, NumericVector target, NumericVector profit, NumericVector ME)
{
  // target and ME's dimensionality is nagent + 1
  int nagent = dividedV.nrow();
  int ntask = dividedV.ncol() / nagent;
  vec<double> Vcontainer;
  vec<task<double, int> > V;
  getV<double, int> (dividedV, Vcontainer, profit, V);
  task<double, int> *t = &V[0];
  vec<int> LB(ntask, 0), UB(ntask, nagent - 1);
  vec<double> sumLB(nagent + 1, 0), sumUB(nagent + 1, (nagent - 1) * ntask);
  for(int i = 0; i < ntask; ++i)
  {
    int *&ird = t[i].ird;
    double *&val = t[i].val;
    int &irdim = ird[LB[i]];
    sumLB[irdim] += val[LB[i]] - LB[i];
  }
  for(int i = 0; i < ntask; ++i)
  {
    int *&ird = t[i].ird;
    double *&val = t[i].val;
    int &irdim = ird[UB[i]];
    sumUB[irdim] += val[UB[i]] - UB[i];
  }


  vec<double> MIN_sumUB(target.size()), MAX_sumLB(target.size());
  for(int i = 0, iend = MIN_sumUB.size(); i < iend; ++i)
  {
    MIN_sumUB[i] = target[i] - ME[i] - sumUB[i];
    MAX_sumLB[i] = target[i] + ME[i] - sumLB[i];
  }


  // int MIN_sumUB_maxDim = std::max_element(MIN_sumUB.begin(), MIN_sumUB.end()) - MIN_sumUB.begin();
  int MAX_sumLB_minDim = 0;
  int MAX_sumLB_2ndMinDim = 1;
  if(MAX_sumLB[MAX_sumLB_minDim] > MAX_sumLB[MAX_sumLB_2ndMinDim])
    std::swap(MAX_sumLB_minDim, MAX_sumLB_2ndMinDim);
  for(int i = 2, iend = MAX_sumLB.size(); i < iend; ++i)
  {
    if(MAX_sumLB[i] < MAX_sumLB[MAX_sumLB_minDim])
    {
      int tmp = MAX_sumLB_minDim;
      MAX_sumLB_minDim = i;
      MAX_sumLB_2ndMinDim = tmp;
    }
    else if(MAX_sumLB[i] < MAX_sumLB[MAX_sumLB_2ndMinDim])
    {
      MAX_sumLB_2ndMinDim = i;
    }
  }


  for(int i = 0, iend = MAX_sumLB.size(); i < iend; ++i)
  {
    Rcout << MAX_sumLB[i] << ", ";
  }
  Rcout << "\n";
  Rcout << "MAX_sumLB_minDim = " << MAX_sumLB_minDim << "\n";


  int MIN_sumUBindVal = MIN_sumUB.back();
  vec<int> taskInd(ntask);
  for(int i = 0; i < ntask; ++i) taskInd[i] = i;
  int boo = findBound003<double, int> (
    nagent, ntask, &V[0], &taskInd[0], &LB[0], &UB[0], MIN_sumUBindVal, &MAX_sumLB[0],
    MAX_sumLB_minDim, MAX_sumLB_2ndMinDim);


  return List::create(Named("LB") = LB, Named("UB") = UB, Named("boo") = boo);
}


// [[Rcpp::export]]
List testFindBound003GAP2(NumericMatrix dividedV, NumericVector targetMAX)
{
  // target and ME's dimensionality is nagent + 1
  int nagent = dividedV.nrow();
  int ntask = dividedV.ncol() / nagent;
  vec<double> Vcontainer;
  vec<task<double, int> > V;
  NumericVector profit(ntask);
  getV<double, int> (dividedV, Vcontainer, profit, V);
  task<double, int> *t = &V[0];
  vec<int> LB(ntask, 0), UB(ntask, nagent - 1);
  vec<double> sumLB(nagent + 1, 0), sumUB(nagent + 1, (nagent - 1) * ntask);
  for(int i = 0; i < ntask; ++i)
  {
    int *&ird = t[i].ird;
    double *&val = t[i].val;
    int &irdim = ird[LB[i]];
    sumLB[irdim] += val[LB[i]] - LB[i];
  }
  for(int i = 0; i < ntask; ++i)
  {
    int *&ird = t[i].ird;
    double *&val = t[i].val;
    int &irdim = ird[UB[i]];
    sumUB[irdim] += val[UB[i]] - UB[i];
  }


  vec<double> // MIN_sumUB(targetMAX.size()),
    MAX_sumLB(targetMAX.size());
  for(int i = 0, iend = MAX_sumLB.size(); i < iend; ++i)
  {
    // MIN_sumUB[i] = target[i] - ME[i] - sumUB[i];
    MAX_sumLB[i] = targetMAX[i] - sumLB[i];
  }


  // int MIN_sumUB_maxDim = std::max_element(MIN_sumUB.begin(), MIN_sumUB.end()) - MIN_sumUB.begin();
  int MAX_sumLB_minDim = 0;
  int MAX_sumLB_2ndMinDim = 1;
  if(MAX_sumLB[MAX_sumLB_minDim] > MAX_sumLB[MAX_sumLB_2ndMinDim])
    std::swap(MAX_sumLB_minDim, MAX_sumLB_2ndMinDim);
  for(int i = 2, iend = MAX_sumLB.size(); i < iend; ++i)
  {
    if(MAX_sumLB[i] < MAX_sumLB[MAX_sumLB_minDim])
    {
      int tmp = MAX_sumLB_minDim;
      MAX_sumLB_minDim = i;
      MAX_sumLB_2ndMinDim = tmp;
    }
    else if(MAX_sumLB[i] < MAX_sumLB[MAX_sumLB_2ndMinDim])
    {
      MAX_sumLB_2ndMinDim = i;
    }
  }


  for(int i = 0, iend = MAX_sumLB.size(); i < iend; ++i)
  {
    Rcout << MAX_sumLB[i] << ", ";
  }
  Rcout << "\n";
  Rcout << "MAX_sumLB_minDim = " << MAX_sumLB_minDim << "\n";


  int MIN_sumUBindVal = targetMAX[nagent] - sumUB[nagent];
  vec<int> taskInd(ntask);
  for(int i = 0; i < ntask; ++i) taskInd[i] = i;
  int boo = findBound003<double, int> (
    nagent, ntask, &V[0], &taskInd[0], &LB[0], &UB[0], MIN_sumUBindVal, &MAX_sumLB[0],
    MAX_sumLB_minDim, MAX_sumLB_2ndMinDim);


  return List::create(Named("LB") = LB, Named("UB") = UB, Named("boo") = boo);
}




template<typename T>
inline T *adrs(void *current)
{
  return (T*)((INT(current) + (sizeof(INT) - 1)) & ~(sizeof(INT) - 1));
}


template<typename indtype, typename T>
inline void minDim01(indtype &minDim0, indtype &minDim1, T *v, indtype vsize)
{
  minDim0 = 0;
  minDim1 = 1;
  if(v[0] > v[1]) std::swap(minDim0, minDim1);
  for(indtype i = 2; i < vsize; ++i)
  {
    if(v[i] < v[minDim0])
    {
      indtype tmp = minDim0;
      minDim0 = i;
      minDim1 = tmp;
    }
    else if(v[i] < v[minDim1]) minDim1 = i;
  }
}




template<typename valtype, typename indtype>
struct gapPAT
{
  indtype position;
  indtype s, send;
  indtype len; // aka ntask
  indtype MIN_sumUBindVal;
  indtype MAX_sumLB_minDim;
  indtype MAX_sumLB_2ndMinDim;
  indtype positionTask;
  indtype *LB, *UB;
  indtype *taskInd;
  valtype *MAX_sumLB;
  double accProfit;


  inline void copyAnother(gapPAT &X, void *Xcntr, void *cntr, indtype nagent) // nagent is the dimensionality, never changes
  // Xcntr is the head address of the value container X resides
  // cntr is the head address of the container
  {
    position = X.position;
    s = X.s;
    send = X.send;
    len = X.len;
    MIN_sumUBindVal = X.MIN_sumUBindVal;
    MAX_sumLB_minDim = X.MAX_sumLB_minDim;
    MAX_sumLB_2ndMinDim = X.MAX_sumLB_2ndMinDim;
    positionTask = X.positionTask;
    accProfit = X.accProfit;
    LB = X.LB - (indtype*)Xcntr + (indtype*)cntr;
    UB = X.UB - (indtype*)Xcntr + (indtype*)cntr;
    taskInd = X.taskInd - (indtype*)Xcntr + (indtype*)cntr;
    MAX_sumLB = X.MAX_sumLB - (valtype*)Xcntr + (valtype*)cntr;
    std::copy(X.LB, X.LB + len, LB);
    std::copy(X.UB, X.UB + len, UB);
    std::copy(X.MAX_sumLB, X.MAX_sumLB + nagent + 1, MAX_sumLB);
    std::copy(X.taskInd, X.taskInd + len, taskInd);
  }


  inline void print(indtype nagent)
  {
    Rcpp::Rcout << "position = " << (int)position << ", s = " <<
      (int)s << ", send = " << (int)send << ", len = " << (int)len << ", ";
    Rcpp::Rcout << "MIN_sumUBindVal = " << int(MIN_sumUBindVal) << ", ";
    Rcpp::Rcout << "MAX_sumLB_minDim = " << int(MAX_sumLB_minDim) << ", ";
    Rcpp::Rcout << "MAX_sumLB_2ndMinDim = " << int(MAX_sumLB_2ndMinDim) << ", ";
    Rcpp::Rcout << "positionTask = " << int(positionTask) << ", ";
    Rcout << "accProfit = " << accProfit << "\n";
    Rcpp::Rcout << "LB = ";
    for(int i = 0; i < len; ++i)
    {
      Rcpp::Rcout << int(LB[i]) << ", ";
    }
    Rcpp::Rcout << "\nUB = ";
    for(int i = 0; i < len; ++i)
    {
      Rcpp::Rcout << int(UB[i]) << ", ";
    }
    Rcpp::Rcout << "\ntaskInd = ";
    for(int i = 0; i < len; ++i)
    {
      Rcpp::Rcout <<int(taskInd[i]) << ", ";
    }
    Rcpp::Rcout << "\nMAX_sumLB = ";
    for(int i = 0; i <= nagent; ++i)
    {
      Rcpp::Rcout << MAX_sumLB[i] << ", ";
    }
    Rcpp::Rcout << "\n";
  }


  // len is the parent's subset size
  inline void copyParentGene(gapPAT &x, indtype nagent) // x is the parent
  {
    len = x.len;
    MIN_sumUBindVal = x.MIN_sumUBindVal;
    MAX_sumLB_minDim = x.MAX_sumLB_minDim;
    MAX_sumLB_2ndMinDim = x.MAX_sumLB_2ndMinDim;
    LB = adrs<indtype> (x.MAX_sumLB + nagent + 1);
    UB = LB + len;
    taskInd = UB + len;
    MAX_sumLB = adrs<valtype> (taskInd + len);
    std::copy(x.LB, x.LB + len, LB);
    std::copy(x.UB, x.UB + len, UB);
    std::copy(x.taskInd, x.taskInd + len, taskInd);
    std::copy(x.MAX_sumLB, x.MAX_sumLB + nagent + 1, MAX_sumLB);
    accProfit = x.accProfit;
  }


  // equavalent to giveBirth(), and len here is still gene from the parent
  // inline indtype grow(valtype *ME, valtype ***M, indtype d, bool useBiSearch, std::ofstream *outfile = nullptr)
  inline indtype grow(task<valtype, indtype> *T, indtype nagent, double optProfit)
  {
    indtype boo = findBound003<valtype, indtype> (
      nagent, len, T, taskInd, LB, UB, MIN_sumUBindVal, MAX_sumLB,
      MAX_sumLB_minDim, MAX_sumLB_2ndMinDim);


    // Rcout << "Bounds found ___________________________________boo = "
    //       << (int)boo << "\n";
    // print(nagent);


    if(boo == 0) return 0;
    if(len == 1) return 3;
    if(boo == 2) return 2;


    if(optProfit > 0)
    {
      double S = accProfit;
      for(indtype i = 0; i < len; ++i)
      {
        for(indtype j = 0; j < nagent; ++j)
        {
          S += T[taskInd[i]].profit[j];
        }
      }
      if(S <= optProfit) return 0;
    }


    // find the slot that has the least gap
    position = 0;
    indtype Min = *UB - *LB;
    for(indtype i = 1; i < len; ++i)
    {
      indtype tmp = UB[i] - LB[i];
      if(Min > tmp)
      {
        Min = tmp;
        position = i;
      }
    }


    s = UB[position];
    send = LB[position];
    positionTask = taskInd[position];
    // Rcout << "s = " << (int)s << ", send = " << (int)send << ", position = "
    //       << (int)position << ", " << "positionTask = " << (int)positionTask
    //       << "\n";


    // MIN_sumUBindVal needs no change? yes
    // Logic: MIN = MIN - v[s], sumUB = sumUB - v[s]
    // Update MAX_sumLB, MAX_sumLB_minDim, MAX_sumLB_2ndMinDim.
    // Logic: MAX = MAX - v[s], sumLB = sumLB - v[send]
    // Logic: MAX_sumLB = MAX_sumLB - v[s] + v[send]
    {
      indtype *&ird = T[positionTask].ird;
      valtype *&val = T[positionTask].val;


      // for(indtype i = 0; i <= nagent; ++i)
      // {
      //   Rcout << "i = " << (int)i;
      //   Rcout << ", ird[s] = " << (int)ird[s];
      //   Rcout << ", ird[send] = " << (int)ird[send];
      //   Rcout << ", MAX_sumLB[i] = " << MAX_sumLB[i];
      //   Rcout << ", val[i] = " << val[i];
      //
      //
      //   if(i == ird[s]) MAX_sumLB[i] -= val[s];
      //   else MAX_sumLB[i] -= s;
      //   if(i == ird[send]) MAX_sumLB[i] += val[send];
      //   else MAX_sumLB[i] += send;
      //
      //
      //   Rcout << ", after, MAX_sumLB[i] = " << MAX_sumLB[i] << "\n";
      // }


      for(indtype i = 0; i <= nagent; ++i)
      {
        MAX_sumLB[i] += send - s;
      }
      indtype &irdsend = ird[send], &irds = ird[s];
      MAX_sumLB[irdsend] -= send - s;
      MAX_sumLB[irds] -= send - s;
      if(irdsend != irds)
      {
        MAX_sumLB[irds] += send - val[s];
        MAX_sumLB[irdsend] += val[send] - s;
      }


      minDim01<indtype, valtype> (
          MAX_sumLB_minDim, MAX_sumLB_2ndMinDim, MAX_sumLB, nagent + 1);
    }


    // Rcout << "___________________________________\n\n";


    accProfit += T[positionTask].profit[s];


    if(position >= len / 2)
    {
      std::copy(LB + position + 1, LB + len, LB + position);
      std::copy(UB + position + 1, UB + len, UB + position);
      std::copy(taskInd + position + 1, taskInd + len, taskInd + position);
    }
    else
    {
      std::copy_backward(LB, LB + position, LB + position + 1);
      std::copy_backward(UB, UB + position, UB + position + 1);
      std::copy_backward(taskInd, taskInd + position, taskInd + position + 1);
      ++LB;
      ++UB;
      ++taskInd;
    }


    --len;
    return 1;
  }


  inline indtype update(task<valtype, indtype> *T, indtype nagent)
  {
    if(s <= send) return 0;
    --s;
    MIN_sumUBindVal += 1;
    // Update MAX_sumLB, MAX_sumLB_minDim, MAX_sumLB_2ndMinDim.
    {
      indtype *&ird = T[positionTask].ird;
      valtype *&val = T[positionTask].val;
      for(indtype i = 0; i <= nagent; ++i)
      {
        MAX_sumLB[i] += 1;
      }
      indtype &irds1 = ird[s + 1], &irds = ird[s];
      MAX_sumLB[irds1] = MAX_sumLB[irds1] - 1 + val[s + 1] - s;
      MAX_sumLB[irds] = MAX_sumLB[irds] - 1 + s + 1 - val[s];


      minDim01<indtype, valtype> (
          MAX_sumLB_minDim, MAX_sumLB_2ndMinDim, MAX_sumLB, nagent + 1);
    }
    accProfit += T[positionTask].profit[s + 1] - T[positionTask].profit[s];
    return 1;
  }
};




template<typename valtype, typename indtype>
inline double asolutionProfit(
    task<valtype, indtype> *tv, indtype *asolution, indtype len)
{
  double S = 0;
  for(indtype i = 0; i < len; ++i)
  {
    S += tv[i].profit[asolution[i]];
  }
  return S;
}


template<typename valtype, typename indtype>
inline double asolutionProfit(
    task<valtype, indtype> *tv, indtype *asolution, indtype len, indtype nagent)
{
  double S = 0;
  for(indtype i = 0; i < len; ++i)
  {
    indtype q = asolution[i] / nagent, r = asolution[i] % nagent;
    S += tv[q].profit[r];
  }
  return S;
}




// return -1 means time is up
template<typename valtype, typename indtype>
signed char TTTstack(
    indtype LEN, indtype nagent,
    task<valtype, indtype> *originalTV, // The very original tv
    indtype *optimalRst, // length is LEN
    double *optimalSolProfit,
    gapPAT<valtype, indtype> *SK,
    gapPAT<valtype, indtype> *&SKback,
    double endTime, bool verbose, std::mutex *mx,
    vec<indtype> &acntr)
{
  if(SKback <= SK) return 0;


  // std::ofstream outfile("proboutput.csv", std::ios_base::app);
  // outfile << "\n\n\nNew target--------------------------------------------------------\n";
  while(true)
  {
    // std::cout << "rstCurrentSize = " << rstCurrentSize << ", ";
    // (SKback - 1)->print(d, dl, du, outfile);
    // outfile << "parent printed ___________________________________\n\n";


    SKback->copyParentGene(*(SKback - 1), nagent);


    // SKback->print(d, dl, du, outfile);
    // outfile << "parent copied ___________________________________\n\n";


    indtype boo = SKback->grow(originalTV, nagent, *optimalSolProfit);


    // outfile << "returned boo =," << (int)boo << "\n-----------";


    // SKback->print(d, dl, du, outfile);
    // outfile << "child grown ___________________________________\n\n";


    // continue to give birth.
    if(boo == 1)
    {
      ++SKback;
      continue;
    }


    // if len in the child becomes 1
    // if(boo == 3 or boo == 2)
    if(boo != 0)
    {
      indtype *common = &*acntr.begin();
      for(int i = 1, iend = SKback - &SK[0]; i < iend; ++i)
      {
        // common[i - 1] = SK[i].s;
        common[i - 1] = SK[i].positionTask * nagent + SK[i].s;
      }
      indtype *st = common + (SKback - &SK[0] - 1);
      for(int i = 0; i < SKback->len; ++i)
      {
        st[i] = nagent * SKback->taskInd[i] + SKback->UB[i];
      }
      // std::copy(SKback->UB, SKback->UB + SKback->len, common + (SKback - &SK[0] - 1));
      // double tmpProfit = asolutionProfit<valtype, indtype> (originalTV, common, LEN);
      double tmpProfit = asolutionProfit<valtype, indtype> (originalTV, common, LEN, nagent);
      mx->lock();
      {
        if(tmpProfit > *optimalSolProfit)
        {
          std::copy(common, common + LEN, optimalRst);
          *optimalSolProfit = tmpProfit;
          // if(verbose) Rcout << "Updated profit = " << tmpProfit << "\n";
        }
      }
      mx->unlock();
    }


    while(true)
    {
      // bool updateBool = (SKback - 1)->update(V, d, dlst, dl, dust, du, profitVec);
      bool updateBool = (SKback - 1)->update(originalTV, nagent);


      // (SKback - 1)->print(d, dl, du, outfile);
      // outfile << "parent updated ___________________________________\n\n";


      if(updateBool != 0) break;
      --SKback;
      if(SKback - SK <= 1)
      {
        return 0; // all the combinations have been tried
      }
    }


    if(double(std::clock()) > endTime) return -1; // global time up
  }
  return 1;
}




template<typename valtype, typename indtype>
struct parMgap: public RcppParallel::Worker
{
  bool verbose;
  indtype len;
  indtype nagent;
  double endTime;
  vec<vec<gapPAT<valtype, indtype> > > &SKgroup;
  vec<gapPAT<valtype, indtype>*> &SKbackGroup;
  task<valtype, indtype> *originalTV;
  indtype *optimalSolution;
  double *optimalSolProfit;
  vec<indtype> *cntr;
  std::mutex *mx;
  dynamicTasking *dT;


  void operator() (std::size_t st, std::size_t end)
  {
    for(;;)
    {
      INT objI = 0;
      if(!dT->nextTaskID(objI)) break;
      // __________________________________________________________________________________________
      // thread function:
      {
        vec<gapPAT<valtype, indtype> > &SK = SKgroup[objI];
        gapPAT<valtype, indtype> *SKbegin = &SK.front(), *&SKback = SKbackGroup[objI];
        signed char tmp = TTTstack(
          len, nagent, originalTV,
          optimalSolution, // length is LEN
          optimalSolProfit,
          SKbegin, SKback,
          endTime, verbose, mx, cntr[st]);
        if(tmp < 0) break;
      }
      // __________________________________________________________________________________________
    }
  }


  parMgap(bool verbose, indtype len, indtype nagent,
          double endTime, vec<vec<gapPAT<valtype, indtype> > > &SKgroup,
          vec<gapPAT<valtype, indtype>*> &SKbackGroup,
          task<valtype, indtype> *originalTV, indtype *optimalSolution,
          double *optimalSolProfit, std::size_t maxCore, int tasks):
    verbose(verbose), len(len), nagent(nagent), endTime(endTime),
    SKgroup(SKgroup), SKbackGroup(SKbackGroup),originalTV(originalTV),
    optimalSolution(optimalSolution), optimalSolProfit(optimalSolProfit)
  {
    std::mutex mxP; mx = &mxP;
    dynamicTasking dtask(maxCore, tasks); dT = &dtask;
    vec<vec<indtype> > cntrs(maxCore, vec<indtype>(len)); cntr = &cntrs[0];
      parallelFor(0, dT->NofCore, *this);
  }
};




template<typename valtype, typename indtype>
void copySKcouple(vec<gapPAT<valtype, indtype> > &SK,
                  unsigned depth, vec<INT> &content,
                  vec<gapPAT<valtype, indtype> > &SKcopy,
                  vec<INT> &contentCopy, indtype nagent)
{
  SKcopy.resize(SK.size());
  contentCopy.resize(content.size());
  for(unsigned i = 0; i < depth; ++i)
  {
    SKcopy[i].copyAnother(SK[i], &content[0], &contentCopy[0], nagent);
  }
}




// spawn from the first level, BFS
template<typename valtype, typename indtype>
int spawn(
    vec<gapPAT<valtype, indtype> > &SK, vec<INT> &SKcontent,
    vec<vec<gapPAT<valtype, indtype> > > &SKfamily,
    vec<vec<INT> > &SKcontentFamily,
    task<valtype, indtype> *originalTV,
    indtype *optimalSolution,
    double &optimalProfit,
    indtype LEN, indtype nagent,
    int maxCore, int threadLoad)
{
  int back = 1;
  SKfamily.resize(1); SKfamily[0].swap(SK);
  SKcontentFamily.resize(1); SKcontentFamily[0].swap(SKcontent);


  while(SKfamily.size() > 0 and
        SKfamily.size() < unsigned(maxCore * threadLoad)
        and back < LEN)
  {
    vec<vec<vec<gapPAT<valtype, indtype> > > > SKfamilyGroups(SKfamily.size());
    vec<vec<vec<INT> > > SKcontentFamilyGroups(SKfamily.size());


    for(int i = 0, iend = SKfamily.size(); i < iend; ++i)
    {
      vec<gapPAT<valtype, indtype> > &SK = SKfamily[i];
      vec<INT> &contentVec = SKcontentFamily[i];
      SK[back].copyParentGene(SK[back - 1], nagent);
      int boo = SK[back].grow(originalTV, nagent, optimalProfit);


      if(boo == 0)
      {
        continue;
      }
      else if(boo == 3 or boo == 2)
      {
        vec<indtype> commonV(LEN);
        indtype *common = &*commonV.begin();
        for(int k = 1; k < back; ++k)
        {
          // common[k - 1] = SK[k].s;
          common[k - 1] = SK[k].positionTask * nagent + SK[k].s;
        }
        indtype *st = common + back - 1;
        // std::copy(SK[back].UB, SK[back].UB + SK[back].len, common + back - 1);
        for(int i = 0; i < SK[back].len; ++i)
        {
          st[i] = SK[back].UB[i] + SK[back].taskInd[i] * nagent;
        }
        // valtype tmpProfit = asolutionProfit(originalTV, common, LEN);
        valtype tmpProfit = asolutionProfit(originalTV, common, LEN, nagent);
        if(tmpProfit > optimalProfit)
        {
          std::copy(common, common + LEN, optimalSolution);
          optimalProfit = tmpProfit;
        }
        continue;
      }


      int siblings = SK[back].s - SK[back].send + 1;


      SKfamilyGroups[i].resize(siblings);
      SKcontentFamilyGroups[i].resize(siblings);


      SKfamilyGroups[i][0].swap(SK);
      SKcontentFamilyGroups[i][0].swap(contentVec);
      for(int k = 1; k < siblings; ++k)
      {
        copySKcouple(
          SKfamilyGroups[i][k - 1], back + 1,
          SKcontentFamilyGroups[i][k - 1],
          SKfamilyGroups[i][k],
          SKcontentFamilyGroups[i][k], nagent);
        SKfamilyGroups[i][k][back].update(originalTV, nagent);
      }


      for(int k = 0; k < siblings; ++k)
      {
        SKfamilyGroups[i][k][back].send = SKfamilyGroups[i][k][back].s;
      }
    }


    int spawnTotal = 0;
    for(int i = 0, iend = SKfamilyGroups.size(); i < iend; ++i)
    {
      spawnTotal += SKfamilyGroups[i].size();
    }
    SKfamily.resize(spawnTotal);
    SKcontentFamily.resize(spawnTotal);


    int k = 0;
    for(int i = 0, iend = SKfamilyGroups.size(); i < iend; ++i)
    {
      for(int j = 0, jend = SKfamilyGroups[i].size(); j < jend; ++j)
      {
        SKfamilyGroups[i][j].swap(SKfamily[k]);
        SKcontentFamilyGroups[i][j].swap(SKcontentFamily[k]);
        ++k;
      }
    }


    ++back;
  }


  return back;
}




template<typename valtype, typename indtype>
// vr has N rows and _d columns
IntegerVector GAPcpp(
    NumericMatrix dividedV,
    NumericVector profitV,
    NumericMatrix MAXmat,
    IntegerVector zeroBasedLB, IntegerVector zeroBasedUB,
    double endTime, int maxCore, int threadLoad,
    bool verbose, bool heuristic)
{
  if(maxCore <= 1) threadLoad = 0;
  int nagent = dividedV.nrow();
  int len = dividedV.ncol() / nagent;
  vec<task<valtype, indtype> > V;
  vec<valtype> Vcontainer;
  getV(dividedV, Vcontainer, profitV, V);
  task<valtype, indtype> *originalTV = &V[0];


  if(false)
  {
    Rcout << "===========================================\n";
    Rcout << "Print task vector:\n";
    Rcout << "Irregular dimension =\n";
    for(int i = 0; i < len; ++i)
    {
      for(int j = 0; j < nagent; ++j)
      {
        Rcout << (int)originalTV[i].ird[j] << ", ";
      }
      Rcout << "\n";
    }
    Rcout << "\n\n";
    Rcout << "Irregular dimension value =\n";
    for(int i = 0; i < len; ++i)
    {
      for(int j = 0; j < nagent; ++j)
      {
        Rcout << originalTV[i].val[j] << ", ";
      }
      Rcout << "\n";
    }
    Rcout << "\n\n";
    Rcout << "Profits =\n";
    for(int i = 0; i < len; ++i)
    {
      for(int j = 0; j < nagent; ++j)
      {
        Rcout << originalTV[i].profit[j] << ", ";
      }
      Rcout << "\n";
    }
    Rcout << "===========================================\n";
  }


  vec<indtype> optimalSolution(len);
  double optimalProfit = 0;


  // Compute the stack size needed
  INT stackLen = 0;
  {
    for(int i = len + 6; i >= 0; --i)
    {
      unsigned tmp = sizeof(indtype) * i * 3;
      unsigned q = tmp / sizeof(INT), r = tmp % sizeof(INT);
      stackLen += q;
      if(r != 0) ++stackLen;


      tmp = sizeof(valtype) * (nagent + 1);
      q = tmp / sizeof(INT); r = tmp % sizeof(INT);
      stackLen += q;
      if(r != 0) ++stackLen;
    }
  }


  int NofTargets = MAXmat.ncol();
  if(verbose) Rcout << "Mining starts. There are " << NofTargets << " tasks:\n";
  for(int i = 0; i < NofTargets; ++i)
  {
    if(verbose) Rcout << i << ", ";
    // Rcout << "====================================\n";
    if(double(std::clock()) > endTime) break;


    vec<INT> content(stackLen);
    vec<gapPAT<valtype, indtype> > SK(len + 6);


    // Fill the 1st stack slot
    {
      SK[0].LB = (indtype*)&content[0];
      SK[0].UB = SK[0].LB + len;
      SK[0].taskInd = SK[0].UB + len;
      SK[0].MAX_sumLB = adrs<valtype> (SK[0].taskInd + len);


      SK[0].len = len;
      SK[0].accProfit = 0;
      for(int i = 0; i < len; ++i)
      {
        SK[0].LB[i] = zeroBasedLB[i];
        SK[0].UB[i] = zeroBasedUB[i];
      }
      for(int i = 0; i < len; ++i) SK[0].taskInd[i] = i;


      valtype *maxPtr = &MAXmat[0] + i * (nagent + 1);
      SK[0].MIN_sumUBindVal = maxPtr[nagent] -
        std::accumulate(SK[0].UB, SK[0].UB + len, unsigned(0));


      indtype *&lb = SK[0].LB;
      std::copy(maxPtr, maxPtr + nagent + 1, SK[0].MAX_sumLB);
      for(int i = 0; i < len; ++i)
      {
        indtype *&ird = originalTV[i].ird;
        valtype *&val = originalTV[i].val;
        for(int j = 0; j <= nagent; ++j) SK[0].MAX_sumLB[j] -= lb[i];
        SK[0].MAX_sumLB[ird[lb[i]]] += lb[i] - val[lb[i]];
      }
      minDim01<indtype, valtype> (
          SK[0].MAX_sumLB_minDim, SK[0].MAX_sumLB_2ndMinDim,
          SK[0].MAX_sumLB, nagent + 1);
    }


    vec<vec<gapPAT<valtype, indtype> > > SKfamily;
    vec<vec<INT> > contentFamily;


    // int back = 0;
    double previousProfit = optimalProfit;
    int back = spawn<valtype, indtype> (
      SK, content, SKfamily, contentFamily, originalTV,
      &optimalSolution[0], optimalProfit, len, nagent,
      maxCore, threadLoad);


    if(optimalProfit > previousProfit)
    {
      if(verbose) Rcout << "Updated profit = " << optimalProfit << "\n";
      if(heuristic) break;
    }


    if(back >= len) continue;


    vec<gapPAT<valtype, indtype>*> SKfamilyBack(SKfamily.size());
    for(int i = 0, iend = SKfamily.size(); i < iend; ++i)
    {
      SKfamilyBack[i] = &SKfamily[i][0] + back;
    }


    previousProfit = optimalProfit;


    // {
    //   for(int i = 0, iend = SKfamily.size(); i < iend; ++i)
    //   {
    //     SKfamily[i][0].print(nagent);
    //     Rcpp::Rcout << "\n\n";
    //   }
    // }


    // Rcpp::Rcout << "SKfamily.size() = " << SKfamily.size() << "\n";


    parMgap<valtype, indtype> (
        verbose, len, nagent, endTime,
        SKfamily, SKfamilyBack, originalTV, &optimalSolution[0],
        &optimalProfit, maxCore, SKfamily.size());


    if(optimalProfit > previousProfit)
    {
      if(verbose) Rcout << "Updated profit = " << optimalProfit << "\n";
      if(heuristic) break;
    }
    // else
    // {
    //   if(verbose) Rcout << "No profit update\n";
    // }
  }


  IntegerVector rst(optimalSolution.begin(), optimalSolution.end());
  return rst;
}




// [[Rcpp::export]]
IntegerVector z_GAP(int maxCore,
                    NumericMatrix dividedV,
                    NumericVector profitV,
                    NumericMatrix MAXmat,
                    IntegerVector zeroBasedLB,
                    IntegerVector zeroBasedUB,
                    double duration, int threadLoad = 8,
                    bool verbose = true,
                    bool heuristic = false)
{
  int N = dividedV.ncol();
  int nagent = dividedV.nrow();
  double endTime = (double)std::clock() + duration * CLOCKS_PER_SEC;
  IntegerVector result;


  if(std::max(N, nagent + 1) < 127)
  {
    result = GAPcpp<double, signed char> (
      dividedV, profitV, MAXmat, zeroBasedLB, zeroBasedUB,
      endTime, maxCore, threadLoad, verbose, heuristic);
  }
  else if(std::max(N, nagent + 1) < 32767)
  {
    result = GAPcpp<double, short> (
      dividedV, profitV, MAXmat, zeroBasedLB, zeroBasedUB,
      endTime, maxCore, threadLoad, verbose, heuristic);
  }
  else
  {
    result = GAPcpp<double, int> (
      dividedV, profitV, MAXmat, zeroBasedLB, zeroBasedUB,
      endTime, maxCore, threadLoad, verbose, heuristic);
  }
  return result;
}


























