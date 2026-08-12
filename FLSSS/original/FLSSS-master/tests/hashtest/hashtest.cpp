// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
#define XXH_STATIC_LINKING_ONLY   /* access advanced declarations */
#define XXH_IMPLEMENTATION   /* access definitions */
#include "../../src/header/xxhash.h"
#include "../../src/header/tiktok.hpp"
#include "../../src/header/flattenUnorderedSet.hpp"
using namespace Rcpp;
#ifndef vec
#define vec std::vector
#endif
#include <bitset>


// [[Rcpp::export]]
void printdoublebits(NumericVector x)
{
  uint64_t y;
  for(int i = 0, iend = x.size(); i < iend; ++i)
  {
    std::memcpy(&y, &x, 8);
    Rcout << std::bitset<64>(y) << ", ";
  }
}



template<typename Hashtable>
void printCollisionsSummary(Hashtable H)
{
  vec<int> rst;
  for(int i = 0, iend = H.bucket_count(); i < iend; ++i)
  {
    int k = H.bucket_size(i);
    if(k > 1) rst.push_back(k);
  }
  std::sort(rst.begin(), rst.end());
  Rcout << "Number of buckets = " << H.bucket_count() << "\n";
  Rcout << "Number of nonsingleton buckets = " << rst.size() << "\n";
  Rcout << "Proportion of nonsingleton buckets = " <<
    rst.size() / (H.bucket_count() + 0.0)  << "\n";
  if(rst.size() > 0)
  {
    Rcout << "Nonsingleton bucket size Min. = " << rst[0] << "\n";
    Rcout << "Nonsingleton bucket size 1stQu. = " <<
      rst[std::min<int> (rst.size() - 1, (int)(std::round(rst.size() * 0.25)))] << "\n";
    Rcout << "Nonsingleton bucket size Median. = " <<
      rst[std::min<int> (rst.size() - 1, (int)(std::round(rst.size() * 0.5)))] << "\n";
    Rcout << "Nonsingleton bucket size Mean. = " <<
      std::accumulate(rst.begin(), rst.end(), 0.0) / rst.size() << "\n";
    Rcout << "Nonsingleton bucket size 3rdQu. = " <<
      rst[std::min<int> (rst.size() - 1, (int)(std::round(rst.size() * 0.75)))] << "\n";
    Rcout << "Nonsingleton bucket size Max. = " <<
      rst.back() << "\n";
  }
}


// Hash an array of Ele.
template<typename Ele>
struct XXH64hashArray
{
  int size;
  uint64_t seed;
  XXH64hashArray(int size, uint64_t seed): size(size), seed(seed){}
  uint64_t operator()(Ele *x) const
  {
    return XXH64(x, sizeof(Ele) * size, seed);
  }
};


// Array equal.
template<typename Ele>
struct ArrayEqual
{
  int arraySize;
  ArrayEqual(int size): arraySize(size){}
  int size() { return arraySize; }
  bool operator()(Ele *x, Ele *y) const
  {
    for(int i = 0; i < arraySize; ++i)
    {
      if(x[i] != y[i]) return false;
    }
    return true;
  }
};



// [[Rcpp::export]]
List test(IntegerMatrix x, IntegerMatrix y, bool verbose = true)
{
  int N = x.ncol(), segSize = x.nrow();
  std::unordered_set<int*, XXH64hashArray<int>, ArrayEqual<int> >
    XXH64Hashset(N * 1.3, XXH64hashArray<int> (segSize, 42),
                 ArrayEqual<int> (segSize));
  tiktok<std::chrono::microseconds> timer;
  timer.tik();
  for(int i = 0, iend = x.size(); i < iend; i += segSize)
    XXH64Hashset.insert(&x[i]);


  if(verbose) Rcout << "XXH64 hash insertion time = " << timer.tok() << "\n";
  if(verbose) Rcout << "Number of elements = " << N << "\n";
  if(verbose) printCollisionsSummary(XXH64Hashset);
  if(verbose) Rcout << "\n";


  if(verbose) Rcout << "See addresses of elements:\n";
  for(int i = 0, iend = XXH64Hashset.bucket_count(); i < iend; ++i)
  {
    if(verbose) Rcout << "Bucket " << i << "'s distance to the 1st bucket = " <<
      (&*XXH64Hashset.begin(i) - &*XXH64Hashset.begin(0)) / sizeof(int) << "\n";
    int k = 0;
    for(auto it = XXH64Hashset.begin(i); it != XXH64Hashset.end(i); ++it, ++k)
    {
      if(verbose) Rcout << "  Element " << k << "'s distance to the 1st element = " <<
        (&*it - &*XXH64Hashset.begin(i)) / sizeof(int) << "\n";
      std::bitset<64> tmpbits(**it);
      if(verbose) Rcout << **it << ", " << tmpbits.to_string() << "\n";
    }
    if(verbose) Rcout << "\n";
  }


  IntegerVector isYin(y.ncol());
  timer.tok();
  for(int i = 0, iend = y.ncol(); i < iend; ++i)
  {
    int *tmp = &y[0] + i * segSize;
    // std::bitset<64> tmpbits(*tmp);
    // Rcout << tmpbits.to_string() << "\n";
    isYin[i] =  XXH64Hashset.find(tmp) != XXH64Hashset.end();
    // if(*tmp == 0 and verbose)
    // {
    //   Rcout << "i = " << i << ", ";
    //   Rcout << "Found addr = " << &*XXH64Hashset.find(tmp) << ", ";
    //   Rcout << "End addr = " << &*XXH64Hashset.end() << ", ";
    //   Rcout << (XXH64Hashset.find(tmp) != XXH64Hashset.end()) << "\n";
    // }
  }
  Rcout << "flattened query time = " << timer.tok() << "\n";


  // template<typename ing, typename num, typename hashFun, typename equalFun>
  FlatUnorderedSetOfArrays<int, int, XXH64hashArray<int>, ArrayEqual<int> >
    flattenH(XXH64hashArray<int>(segSize, 42), ArrayEqual<int> (segSize));
  List flattenedH = flattenH.exportFlattened(XXH64Hashset);


  flattenH.importFlattened(flattenedH);
  IntegerVector isYinViaFlattened(y.ncol());
  timer.tik();
  for(int i = 0, iend = y.ncol(); i < iend; ++i)
  {
    int *tmp = &y[0] + i * y.nrow();
    isYinViaFlattened[i] = flattenH.isin(tmp);
  }
  Rcout << "flattened query time = " << timer.tok() << "\n";


  return List::create(Named("isYin") = isYin,
                      Named("isYinViaFlattened") = isYinViaFlattened,
                      Named("flattenedH") = flattenedH);
}


#undef XXH_STATIC_LINKING_ONLY   /* access advanced declarations */
#undef XXH_IMPLEMENTATION   /* access definitions */










