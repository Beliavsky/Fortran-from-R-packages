# pragma once
# include "macros.hpp"
# include "dnyTasking.hpp"
# include <Rcpp.h>


template<typename valtype, typename indtype, bool timeConstraint>
struct dp01kp
{
  indtype tmpj, minCost_1, *w;
  valtype **value, *v;
  INT iter;
  double endTime;
  void assign(indtype minCost_1, indtype *w, valtype **value, valtype *v, double endTime)
  {
    this->minCost_1 = minCost_1;
    this->w = w;
    this->value = value;
    this->v = v;
    this->endTime = endTime;
    iter = 0;
  }
  valtype run(indtype i, indtype j)
  {
    if(timeConstraint)
    {
      ++iter;
      if((iter & (64 - 1)) == 0)
      {
        if(double(std::clock()) > endTime) return 0;
      }
    }


    if(i == 0 or j <= minCost_1) return 0;
    // Rcpp::Rcout << i << ", " << j << ", " << minCost_1 << ", " << j - minCost_1 << "\n";
    valtype tmp1 = 0;
    if(value[i - 1][j - minCost_1] != -std::numeric_limits<valtype>::max())
    {
      tmp1 = value[i - 1][j - minCost_1];
    }
    else tmp1 = run(i - 1, j);
    // Rcpp::Rcout << "tmp1 = " << tmp1 << "\n";


    // Rcpp::Rcout << "j < w[i - 1] == " << int(j < w[i - 1]) << "\n";
    if(j < w[i - 1])
    {
      // if(j - minCost_1 >= 0)
      value[i][j - minCost_1] = tmp1;
      return tmp1;
    }


    tmpj = j - minCost_1 - w[i - 1];
    // Rcpp::Rcout << "tmpj == " << tmpj << "\n";
    if(tmpj >= 0 and value[i - 1][tmpj] != -std::numeric_limits<valtype>::max())
    {
      value[i][j - minCost_1] = std::max<valtype> (
        tmp1, value[i - 1][tmpj] + v[i - 1]);
    }
    else
    {
      value[i][j - minCost_1] = std::max<valtype> (
        tmp1, run(i - 1, j - w[i - 1]) + v[i - 1]);
    }
    // Rcpp::Rcout << "value[i][j - minCost_1] == " << value[i][j - minCost_1] << "\n";
    return value[i][j - minCost_1];
  }
};




template<typename valtype, typename indtype>
inline void getSelection(valtype **value, indtype i, indtype j, indtype *w,
                  valtype *v, unsigned char *selection, indtype minCost_1)
{
  // j is the total weight limit
  for(; i > 0; --i)
  {
    // if comming from tmp1: value[i - 1][j]
    // if coming from tmp2: value[i - 1][j - w[i - 1]] + v[i - 1];
    if(j >= minCost_1 and
       value[i][j - minCost_1] != value[i - 1][j - minCost_1] and
       j >= w[i - 1])
    {
      selection[i - 1] = true;
      j -= w[i - 1];
    }
  }
}



