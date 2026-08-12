# include <Rcpp.h>
# include "header/triMat.hpp"
# include "header/mvalFindBound.hpp"
# include "header/macros.hpp"
using namespace Rcpp;




// Do not delete this!
// [[Rcpp::export]]
List z_findBound(int len, NumericMatrix V, NumericVector target,
                 NumericVector me, IntegerVector initialLB = -1,
                 IntegerVector initialUB = -1, int findBoundTimes = 1,
                 bool useBinarySearch = 0, bool UBfirst = 0)
{
  int d = V.ncol();
  int vlen = V.nrow();


  vec<unsigned char> computeSpace((std::size_t)len * (2 * vlen - len + 1) / 2 *
    d * sizeof(double) + len * sizeof(double**) + (std::size_t)len *
    (2 * vlen - len + 1) / 2 * sizeof(double*));


  triM <double, int> M;

  //longjmp(env,1);

  M.make(&computeSpace.front(), len, V);
  double **v = M[0];


  double *x = &target[0], *ME = &me[0];


  IntegerVector LB, UB;
  if(initialLB[0] < 0 or initialUB[0] < 0)
  {
    IntegerVector tmp(len), tmp2(len);
    for(int i = 0, iend = tmp.size(); i != iend; ++i)
    {
      tmp[i] = i;
      tmp2[i] = vlen - len + i;
    }
    LB = tmp;
    UB = tmp2;
  }
  else
  {
    LB = initialLB - 1;
    UB = initialUB - 1;
  }


  vec<double> sumLB(d, 0), sumUB(d, 0);


  for(int i = 0; i < len; ++i)
  {
    mvalPlus(&sumLB[0], &sumLB[0], v[LB[i]], d);
    mvalPlus(&sumUB[0], &sumUB[0], v[UB[i]], d);
  }


  int findOrNot = 0;
  // unsigned totaltime = 0;


  vec<int> LBres(LB.begin(), LB.end()), UBres(UB.begin(), UB.end());
  vec<double> sumLBres(sumLB.begin(), sumLB.end()), sumUBres(sumUB.begin(), sumUB.end());


  // double Min[d], Max[d];
  vec<double> acntr((int)2 * d);
  double *Min = &*acntr.begin(), *Max = Min + d;
  for(int i = 0, iend = d; i < iend; ++i)
  {
    Min[i] = x[i] - ME[i];
    Max[i] = x[i] + ME[i];
  }
  vec<double> SRVcntr(d);
  for(int I = 0; ;)
  {
    // std::clock_t t = std::clock();
    // totaltime += std::clock() - t;


    if(UBfirst) findOrNot = findBoundUpFirstCpp<double, int, 0, 0> (
      len, d, 0, d, 0, d, Min, Max, &LB[0], &sumLB[0], &UB[0], &sumUB[0], M.mat,
      nullptr, SRVcntr);
    else findOrNot = findBoundCpp<double, int, 0, 0> (
      len, d, 0, d, 0, d, Min, Max, &LB[0], &sumLB[0], &UB[0], &sumUB[0], M.mat,
      nullptr, SRVcntr);


    ++I;
    if(I == findBoundTimes)break;


    std::copy(LBres.begin(), LBres.end(), &LB[0]);
    std::copy(UBres.begin(), UBres.end(), &UB[0]);
    std::copy(sumLBres.begin(), sumLBres.end(), &sumLB[0]);
    std::copy(sumUBres.begin(), sumUBres.end(), &sumUB[0]);
  }


  // Rcout << "timeCost = " << totaltime << "\n";


  return List::create(findOrNot, LB + 1, UB + 1);
}








//-------------------------------------------------------------------------------------------------
// Do not delete this!
// [[Rcpp::export]]
List z_findBoundIntegerized(
    int len, NumericMatrix V, NumericVector mask,
    NumericVector target, NumericVector me,
    IntegerVector initialLB = -1, IntegerVector initialUB = -1,
    int findBoundTimes = 1, bool useBinarySearch = 0, bool UBfirst = 0)
{
  int d = V.ncol();
  int vlen = V.nrow();


  vec<unsigned char> computeSpace((std::size_t)len * (2 * vlen - len + 1) /
    2 * d * sizeof(INT) + len * sizeof(INT**) + (std::size_t)len * (
        2 * vlen - len + 1) / 2 * sizeof(INT*));


  triM<INT, int> M;

  //longjmp(env,1);

  M.make(&computeSpace.front(), len, V);
  INT **v = M[0];


  INT *x = (INT*)&target[0], *ME = (INT*)&me[0];


  IntegerVector LB, UB;
  if(initialLB[0] < 0 or initialUB[0] < 0)
  {
    IntegerVector tmp(len), tmp2(len);
    for(int i = 0, iend = tmp.size(); i != iend; ++i)
    {
      tmp[i] = i;
      tmp2[i] = vlen - len + i;
    }
    LB = tmp;
    UB = tmp2;
  }
  else
  {
    LB = initialLB - 1;
    UB = initialUB - 1;
  }


  vec<INT> sumLB(d, 0), sumUB(d, 0);


  for(int i = 0; i < len; ++i)
  {
    mvalPlus(&sumLB[0], &sumLB[0], v[LB[i]], d);
    mvalPlus(&sumUB[0], &sumUB[0], v[UB[i]], d);
  }


  int findOrNot = 0;
  // unsigned totaltime = 0;


  vec<int> LBres(LB.begin(), LB.end()), UBres(UB.begin(), UB.end());
  vec<INT> sumLBres(sumLB.begin(), sumLB.end()), sumUBres(sumUB.begin(), sumUB.end());


  // INT Min[d], Max[d];
  vec<INT> minmax(d * (int)2);
  INT *Min = &*minmax.begin(), *Max = Min + d;
  for(int i = 0, iend = d; i < iend; ++i)
  {
    Min[i] = x[i] - ME[i];
    Max[i] = x[i] + ME[i];
  }
  vec<INT> SRVcntr(d);
  for(int I = 0; ;)
  {
    // std::clock_t t = std::clock();
    // totaltime += std::clock() - t;


    if(UBfirst) findOrNot = findBoundUpFirstCpp<INT, int, 1, 0> (
      len, d, 0, d, 0, d, Min, Max, &LB[0], &sumLB[0], &UB[0], &sumUB[0], M.mat, (INT*)&mask[0], SRVcntr);
    else findOrNot = findBoundCpp<INT, int, 1, 0> (
      len, d, 0, d, 0, d, Min, Max, &LB[0], &sumLB[0], &UB[0], &sumUB[0], M.mat, (INT*)&mask[0], SRVcntr);


    ++I;
    if(I == findBoundTimes) break;


    std::copy(LBres.begin(), LBres.end(), &LB[0]);
    std::copy(UBres.begin(), UBres.end(), &UB[0]);
    std::copy(sumLBres.begin(), sumLBres.end(), &sumLB[0]);
    std::copy(sumUBres.begin(), sumUBres.end(), &sumUB[0]);
  }


  // Rcout << "timeCost = " << totaltime << "\n";


  return List::create(findOrNot, LB + 1, UB + 1);
}































