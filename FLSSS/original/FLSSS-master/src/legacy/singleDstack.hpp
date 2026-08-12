// # include <Rcpp.h>
// # include <fstream>
# include "PATclass.hpp"
# include "oneDoperation.hpp"
// using namespace Rcpp;


namespace legacy {


template<typename valtype, typename indtype>
inline int TTTstack(
    indtype LEN, indtype N, valtype **M, valtype ME,
    vec<vec<indtype> > &result, int sizeNeed, std::size_t durationLimit,
    PAT<valtype, indtype> *SK, PAT<valtype, indtype> *SKback,
    bool useBisearchInFindBounds)
{
  if(SKback <= SK) return SKback - SK;


  valtype *V = M[0];
  if(LEN == 1)
  {
    for(indtype i = 0; i < N; ++i)
    {
      if(V[i] <= SK->target + ME and V[i] >= SK->target - ME)
      {
        result.push_back(vec<indtype>(1, i));
      }
    }
    return SKback - SK;
  }


  std::clock_t timeEnd = std::clock() + durationLimit;


  // std::ofstream outfile("proboutput.csv");
  while(true)
  {
    // (SKback - 1)->print(outfile);
    // outfile << "parent printed ___________________________________\n\n";
    SKback->copyParentGene(*(SKback - 1));
    // SKback->print(outfile);
    // outfile << "parent copied ___________________________________\n\n";


    indtype boo = SKback->grow(ME, M, useBisearchInFindBounds);
    // SKback->print(outfile);
    // outfile << "child grown ___________________________________\n\n";


    // continue to give birth.
    if(boo == 1)
    {
      ++SKback;
      continue;
    }


    // if len in the child becomes 1
    if(boo == 3)
    {
      vec<indtype> common(LEN);
      indtype *commoni = &common.front();


      for(PAT<valtype, indtype> *SKi = SK + 1; SKi < SKback; ++SKi, ++commoni)
      {
        *commoni = SKi->s;
      }


      indtype i = SKback->LB[0], iend = SKback -> UB[0] + 1;


      for(; i < iend; ++i)
      {
        common.back() = i;
        result.push_back(common);
      }
    }
    else if(boo == 2)
    {
      result.resize(result.size() + 1);
      result.back().resize(LEN);
      indtype *x = &result.back().front();


      for(PAT<valtype, indtype> *SKi = SK + 1; SKi < SKback; ++SKi, ++x)
      {
        *x = SKi->s;
      }
      for(indtype i = 0, iend = SKback->len; i < iend; ++i)
      {
        x[i] = SKback->UB[i];
      }
    }


    while(true)
    {
      bool updateBool = (SKback - 1)->update(M);
      // (SKback - 1)->print(d, outfile);
      // outfile << "parent updated ___________________________________\n\n";
      if(updateBool != 0) break;
      --SKback;
      if(SKback - SK <= 1)
      {
        return 0; // All the combinations have been tried
      }
    }


    if(result.size() >= (unsigned)sizeNeed or std::clock() > timeEnd) break;
  }


  return SKback - SK;
}




}
