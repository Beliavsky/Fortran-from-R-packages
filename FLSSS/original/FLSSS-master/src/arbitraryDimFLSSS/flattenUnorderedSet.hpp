#pragma once
#include <Rcpp.h>
#include <unordered_set>
using namespace Rcpp;


// Flatten unordered_set of fixed-size arrays.
// H is an unordered_set with num* as the key type.
// hashFun and equalFun are template parameters of the unordered_set.
// equalFun MUST have member function .size() returning the array size.
// ing MUST chosen according to the total number of nums in table, namely
// arraySize * the number of arrays.
template<typename ing, typename num, typename hashFun, typename equalFun>
struct FlatUnorderedSetOfArrays
{
  std::size_t indexSize, contentSize; // # buckets, # arrays.
  std::size_t arraySize; // # elements in an array.
  std::pair<ing, ing> *index;
  num *content;
  hashFun hf;
  equalFun ef;
  FlatUnorderedSetOfArrays(hashFun &&hf, equalFun &&ef): hf(hf), ef(ef) {}


  List exportFlattened(std::unordered_set<num*, hashFun, equalFun> &h)
  {
    indexSize = h.bucket_count();
    arraySize = ef.size();
    contentSize = h.size();
    std::size_t indexVbytes = sizeof(std::pair<ing, ing>) * indexSize;
    std::size_t contentBytes = sizeof(num) * arraySize * contentSize;
    RawVector Hstatic(indexVbytes + contentBytes);
    index = (std::pair<ing, ing>*)(&Hstatic[0]);
    content = (num*)(index + indexSize);
    for(ing i = 0, iend = h.bucket_count(), k = 0; i < iend; ++i)
    {
      index[i].first = k * arraySize;
      for(auto it = h.begin(i); it != h.end(i); ++it, ++k)
      {
        num *a = *it;
        memcpy(content + arraySize * k, a, sizeof(num) * arraySize);
      }
      index[i].second = k * arraySize;


      // Rcout << "index[" << i << "].first = " << index[i].first << ", ";
      // Rcout << "index[" << i << "].second = " << index[i].second << "\n";
    }


    // Rcout << "\ncontent = ";
    // for(int u = 0, uend = contentSize * arraySize; u < uend; ++u)
    //   Rcout << content[u] << ", ";
    // Rcout << "\n";


    index = nullptr;
    content = nullptr;
    std::size_t tmp[3] = {indexSize, contentSize, arraySize};
    NumericVector sizeInfo(3);
    memcpy(&sizeInfo[0], tmp, sizeof(std::size_t) * 3);
    sizeInfo.attr("class") = "integer64";
    return List::create(Named("sizeInfo") = sizeInfo,
                        Named("Hstatic") = Hstatic);
  }


  void importFlattened(List X)
  {
    NumericVector sizeInfo = X[0];
    std::size_t *tmp = (std::size_t*)(&sizeInfo[0]);
    indexSize = tmp[0]; contentSize = tmp[1]; arraySize = tmp[2];
    RawVector H = X[1];
    index = (std::pair<ing, ing>*)(&H[0]);
    content = (num*)(index + indexSize);
  }


  bool isin(num *x)
  {
    std::size_t I = hf(x) % indexSize;
    if(index[I].first == index[I].second) return false;
    num *begin = content + index[I].first, *end = content + index[I].second;
    for(; begin < end; begin += arraySize)
    {
      if(ef(begin, x)) return true;
    }
    return false;
  }
};











































