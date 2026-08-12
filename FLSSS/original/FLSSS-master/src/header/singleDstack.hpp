// # include <Rcpp.h>
# include <fstream>
# include "PATclass.hpp"
# include "oneDoperation.hpp"




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
      // if(V[i] <= SK->target + ME and V[i] >= SK->target - ME)
      if(std::abs(V[i] - SK->target) <= ME)
      {
        result.push_back(vec<indtype> (1, i));
      }
    }
    return SKback - SK;
  }


  std::clock_t timeEnd = std::clock() + durationLimit;


  // std::ofstream outfile("proboutput.csv");
  vec<indtype> hopeV(LEN);
  indtype *hope = &hopeV[0];
  while(true)
  {
    // (SKback - 1)->print();
    // Rcpp::Rcout << "parent printed ___________________________________\n\n";
    SKback->copyParentGene(*(SKback - 1));
    // SKback->print();
    // Rcpp::Rcout << "parent copied ___________________________________\n\n";


    indtype boo = SKback->grow(M, ME, hope, useBisearchInFindBounds // , &outfile
                               );
    // SKback->print();
    // Rcpp::Rcout << "child grown ___________________________________\n\n";


    // continue to give birth.
    if(boo == 1)
    {
      ++SKback;
      continue;
    }


    // if len in the child becomes 1
    if(boo == 3)
    {
      indtype i = SKback->LB[0], iend = SKback -> UB[0] + 1;
      for(; i < iend; ++i)
      {
        hopeV.back() = i;
        result.push_back(hopeV);
      }
    }
    else if(boo == 2)
    {
      std::copy(SKback->UB, SKback->UB + SKback->len, hope);
      result.push_back(hopeV);
    }


    while(true)
    {
      bool updateBool = (SKback - 1)->update(M);
      // (SKback - 1)->print();
      // Rcpp::Rcout << "parent updated ___________________________________\n\n";


      if(updateBool != 0) break;
      hope -= (SKback - 1)->Nzeroed;
      --SKback;
      if(SKback - SK <= 1)
      {
        return 0; // All the combinations have been tried
      }
    }
    if(result.size() >= unsigned(sizeNeed) or std::clock() > timeEnd) break;
  }


  return SKback - SK;
}



