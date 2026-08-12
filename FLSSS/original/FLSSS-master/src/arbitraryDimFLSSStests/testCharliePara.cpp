// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
using namespace Rcpp;
// #include "../arbitraryDimFLSSS/charlieParaFor.hpp" // Contain global pointer
// // collection variable `gp`;
#include "../arbitraryDimFLSSS/charlieThreadPool.hpp" // Contain global pointer
// collection variable `gp`;


#ifndef vec
#define vec std::vector
#endif


// // [[Rcpp::export]]
template<typename num>
num testCharlieParaCpp(num *x, num *y, int size,
                       double discount, int Nscaler, int maxCore)
{
  auto emptyFun = [](std::size_t t)->bool { return false; };
  CharlieThreadPool tp(maxCore);


  tp.parFor(0, size, [&](std::size_t i, std::size_t t)->bool
  {
    auto kend = int((1 + x[i] + y[i]) * Nscaler);
    for (int k = 0; k < kend; ++k)
    {
      x[i] = (x[i] + y[i]) * y[i] * discount;
      y[i] = (y[i] + x[i]) * x[i] * discount;
    }
    return false; // No eartly return;
  }, 1000, emptyFun, emptyFun);


  vec<num> thesum(maxCore, 0);
  num *z = &thesum[0];
  tp.parFor(0, size, [&](std::size_t i, std::size_t t)->bool
  {
    auto kend = int((1 + x[i] + y[i]) * Nscaler);
    for (int k = 0; k < kend; ++k)
    {
      z[t] += (x[i] + y[i]) * y[i] * discount;
      z[t] += (y[i] + x[i]) * x[i] * discount;
    }
    return false; // No eartly return;
  }, 1000, emptyFun, emptyFun);


  return std::accumulate(thesum.begin(), thesum.end(), num(0));
}


// [[Rcpp::export]]
double testCharliePara(NumericVector x, NumericVector y,
                       double discount = 10, int Nscaler = 30000,
                       int maxCore = 8)
{
  return testCharlieParaCpp<double> (&x[0], &y[0], x.size(),
                                     discount, Nscaler, maxCore);
}



/*
void **globalPtr = nullptr;
template<typename num>
num testCharlieParaCpp(num *x, num *y, int size, int maxCore)
{
  void *xy[3] = {(void*)x, (void*)y, nullptr};
  globalPtr = &xy[0];
  auto emptyFun = [](std::size_t t)->bool { return false; }; // No early return;
  ParaFor(0, size, [](std::size_t i, std::size_t t)->bool
  {
    num *x = (num*)(globalPtr[0]), *y = (num*)(globalPtr[1]);
    // int kend = (1 + gp.x[i] + gp.y[i]) * 3e4;
    int kend = (1 + x[i] + y[i]) * 3e4;
    for (int k = 0; k < kend; ++k)
    {
      // gp.x[i] = (gp.x[i] + 1) * gp.y[i];
      // gp.y[i] = (gp.y[i] + 1) * gp.x[i];
      x[i] = (x[i] + 1) * y[i];
      y[i] = (y[i] + 1) * x[i];
    }
    return false; // No eartly return;
  }, 1000, emptyFun, emptyFun);


  vec<double> thesum(maxCore, 0);
  // gp.z = &thesum[0];
  xy[2] = (void*)(&thesum[0]);
  ParaFor(0, size, [](std::size_t i, std::size_t t)->bool
  {
    // int kend = (1 + gp.x[i] + gp.y[i]) * 3e4;
    num *x = (num*)(globalPtr[0]),
      *y = (num*)(globalPtr[1]), *z = (num*)(globalPtr[2]);
    int kend = (1 + x[i] + y[i]) * 3e4;
    for (int k = 0; k < kend; ++k)
    {
      // gp.z[t] += (gp.x[i] + 1) * gp.y[i];
      // gp.z[t] += (gp.y[i] + 1) * gp.x[i];
      z[t] += (x[i] + 1) * y[i];
      z[t] += (y[i] + 1) * x[i];
    }
    return false; // No eartly return;
  }, 1000, emptyFun, emptyFun);


  return std::accumulate(thesum.begin(), thesum.end(), num(0));
}


// [[Rcpp::export]]
double testCharliePara(NumericVector x, NumericVector y,
                       int maxCore = 8)
{
  CharlieThreadPool tpool(maxCore); // Is a global variable baby.
  return testCharlieParaCpp(&x[0], &y[0], x.size(), maxCore);
}
*/




/*
// [[Rcpp::export]]
double testCharliePara(NumericVector x, NumericVector y,
                       int maxCore = 8)
{
  // Uncomment double *x, *y, *z in the global pointer collection.
  CharlieThreadPool tpool(maxCore); // Is a global variable baby.
  // while (true);
  gp.x = &x[0]; // Global shared pointer.
  gp.y = &y[0];
  auto emptyFun = [](std::size_t t)->bool { return false; }; // No early return;
  ParaFor(0, x.size(), [](std::size_t i, std::size_t t)->bool
  {
    int kend = (1 + gp.x[i] + gp.y[i]) * 3e4;
    for (int k = 0; k < kend; ++k)
    {
      gp.x[i] = (gp.x[i] + 1) * gp.y[i];
      gp.y[i] = (gp.y[i] + 1) * gp.x[i];
    }
    return false; // No eartly return;
  }, 1000, emptyFun, emptyFun);


  vec<double> thesum(maxCore, 0);
  gp.z = &thesum[0];
  ParaFor(0, x.size(), [](std::size_t i, std::size_t t)->bool
  {
    int kend = (1 + gp.x[i] + gp.y[i]) * 3e4;
    for (int k = 0; k < kend; ++k)
    {
      gp.z[t] += (gp.x[i] + 1) * gp.y[i];
      gp.z[t] += (gp.y[i] + 1) * gp.x[i];
    }
    return false; // No eartly return;
  }, 1000, emptyFun, emptyFun);


  return std::accumulate(thesum.begin(), thesum.end(), 0.0);


}
*/









