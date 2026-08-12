# include <Rcpp.h>
# include <bitset>
# include "header/integerize.hpp"
// # include "gapH/GAPfindBound.hpp"
using namespace Rcpp;


// [[Rcpp::export]]
List z_integerize(int len, NumericMatrix V, NumericVector target,
                  NumericVector ME, IntegerVector precisionLevel)
{
  int d = target.size();
  int N = V.nrow();
  vec<INT> rst(d * N);
  IntegerVector rstTarget(d), rstME(d);
  double *v = &V[0];
  integerize<double, int> (
      &rst[0], &rstTarget[0], &rstME[0], v, len, N, d,
      &target[0], &ME[0], &precisionLevel[0]);
  IntegerMatrix rstint(N, d);
  std::copy(rst.begin(), rst.end(), rstint.begin());
  return List::create(Named("integerized") = rstint,
                      Named("target") = rstTarget,
                      Named("ME") = rstME);
}


inline int z_bitSizeWithAdditionalSign(int x)
{
  x = std::abs(x);
  int size = 0;
  while(x != 0)
  {
    x >>= 1;
    ++size;
  }
  return size + 1;
}


// [[Rcpp::export]]
List z_which64intAndSize(IntegerVector largestSubsetSum)
{
  int *v = &largestSubsetSum[0];
  IntegerVector which64int(largestSubsetSum.size());
  IntegerVector bitSize(which64int.size());
  int current64intWhich = 0;
  int one64intBits = sizeof(INT) * 8;
  int accBits = 0;
  int dim = largestSubsetSum.size();
  for(int i = 0; i < dim; ++i)
  {
    bitSize[i] = z_bitSizeWithAdditionalSign(v[i]);
    accBits += bitSize[i];
    if(accBits > one64intBits)
    {
      ++current64intWhich;
      accBits = bitSize[i];
    }
    which64int[i] = current64intWhich;
  }
  return List::create(Named("which64int") = which64int, Named("bitSize") = bitSize);
}


// [[Rcpp::export]]
NumericMatrix z_collapseTo64int(IntegerMatrix x, IntegerVector which64int, IntegerVector bitSize)
{
  int dim = which64int.size();
  int N = x.ncol();
  int dim64int = *(which64int.end() - 1) + 1;
  NumericMatrix rst(dim64int, N);
  int int64bits = sizeof(INT) * 8;
  for(int i = 0; i < N; ++i)
  {
    INT *v64int = (INT*)&rst[0] + i * dim64int;
    int *vint = &x[0] + i * dim;
    int accBits = 0;
    for(int k = 0; k < dim; ++k)
    {
      accBits += bitSize[k];
      if(accBits > int64bits)
      {
        accBits = bitSize[k];
      }


      std::size_t tmpINT = 0;
      if(vint[k] >= 0) tmpINT = vint[k];
      else tmpINT = INT(0) - INT(0 - vint[k]);
      v64int[which64int[k]] += tmpINT << (int64bits - accBits); // examined, should be this
    }
  }
  return rst;
}


// [[Rcpp::export]]
NumericVector z_mask(IntegerVector which64int, IntegerVector bitSize)
{
  int dim = which64int.size();
  int dim64int = *(which64int.end() - 1) + 1;
  int int64bits = sizeof(INT) * 8;
  NumericVector masks(dim64int, 0.0);
  INT *m = (INT*)&masks[0];
  int accBits = 0;
  for(int k = 0; k < dim; ++k)
  {
    accBits += bitSize[k];
    if(accBits > int64bits)
    {
      accBits = bitSize[k];
    }
    m[which64int[k]] += INT(1) << (bitSize[k] - 1) << (int64bits - accBits);
  }
  return masks;
}


// // [[Rcpp::export]]
// int archr()
// {
//   return sizeof(INT);
// }


// // [[Rcpp::export]]
// void z_printBinary(NumericMatrix X)
// {
//   INT *x = (INT*)&X[0];
//   int dim = X.nrow();
//   int N = X.ncol();
//   for(int i = 0; i < N; ++i)
//   {
//     INT *v = x + i * dim;
//     for(int j = 0; j < dim; ++j)
//     {
//       Rcout << std::bitset<64>(v[j]) << " ";
//     }
//     Rcout << "\n";
//   }
// }




// // [[Rcpp::export]]
// IntegerVector z_fromStringtoInt(std::string x, IntegerVector sizes)
// {
//   std::size_t v = 0;
//   std::size_t mp = 1;
//   for(int i = x.size() - 1; i >= 0; --i)
//   {
//     v += mp * (x[i] - 48);
//     mp *= 10;
//   }
//   vec<int> tmp(sizes.size());
//   std::partial_sum(sizes.begin(), sizes.end(), tmp.begin());
//   IntegerVector rst(tmp.size());
//   rst[0] = v >> (64 - sizes[0]);
//   if((rst[0] & (1 << sizes[0])) != 0) rst[0] = (~rst[0]) + 1;
//   for(int i = 1, iend = rst.size(); i < iend; ++i)
//   {
//     rst[i] = (v << tmp[i - 1]) >> (64 - sizes[i]);
//     if((rst[i] & (1 << sizes[i])) != 0) rst[i] = (~rst[i]) + 1;
//   }
//   return rst;
// }




// // [[Rcpp::export]]
// IntegerVector z_fromNumtoInt(double x, IntegerVector sizes)
// {
//   std::size_t *y = (std::size_t*)&x, v = *y;
//   vec<int> tmp(sizes.size());
//   std::partial_sum(sizes.begin(), sizes.end(), tmp.begin());
//   IntegerVector rst(tmp.size());
//   rst[0] = v >> (64 - sizes[0]);
//   for(int i = 1, iend = rst.size(); i < iend; ++i)
//   {
//     rst[i] = (v << tmp[i - 1]) >> (64 - sizes[i]);
//   }
//   return rst;
// }































