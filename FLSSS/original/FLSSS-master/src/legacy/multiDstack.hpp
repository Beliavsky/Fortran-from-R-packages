# include <RcppParallel.h>
# include "mPATclass.hpp"




/*
template<typename valtype, typename indtype>
inline std::size_t TTTstack(
    indtype LEN, indtype N,
    indtype d, indtype dlst, indtype dl, indtype dust, indtype du,
    valtype ***M, vec<vec<indtype> > &result, int sizeNeeded,
    mPAT<valtype, indtype> *SK, mPAT<valtype, indtype> *SKback,
    bool useBisearchInFindBounds,
    std::atomic<int> &totalSize, double endTime)
{
  if(SKback <= SK) return SKback - SK;


  int rstCurrentSize = result.size();


  valtype **V = M[0];
  if(LEN == 1)
  {
    for(indtype i = 0; i < N; ++i)
    {
      bool allBetween = 1;
      for(indtype k = 0; k < dl and allBetween; ++k)
      {
        if(V[i][k + dlst] < SK->MIN[k]) allBetween = 0;
      }
      for(indtype k = 0; k < du and allBetween; ++k)
      {
        if(V[i][k + dust] > SK->MAX[k]) allBetween = 0;
      }


      if(allBetween)
      {
        result.push_back(vec<indtype>(1, i));
      }
    }
    // update totalSize
    {
      int addSize = result.size() - rstCurrentSize;
      if(addSize > 0) totalSize.fetch_add(addSize);
    }
    return SKback - SK;
  }


  // std::ofstream outfile("proboutput.csv");
  // outfile << "output\n\n\n";
  // return 0;
  vec<indtype> hopeV(LEN);
  // indtype *hope = &hopeV[0];
  while(true)
  {
    // outfile << "rstCurrentSize = " << rstCurrentSize << ", ";
    // (SKback - 1)->print(d, dl, du, outfile);
    // outfile << "parent printed ___________________________________\n\n";


    SKback->copyParentGene(*(SKback - 1), d, dl, du);


    // SKback->print(d, dl, du, outfile);
    // outfile << "parent copied ___________________________________\n\n";


    indtype boo = SKback->grow(M, d, dlst, dl, dust, du, useBisearchInFindBounds // , &outfile
                                 );


    // SKback->print(d, dl, du, outfile);
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
      indtype i = SKback->LB[0], iend = SKback -> UB[0] + 1;
      for(; i < iend; ++i)
      {
        hopeV.back() = i;
        result.push_back(hopeV);
      }
    }
    else if(boo == 2)
    {
      std::copy(SKback->UB, SKback->UB + SKback->len, &hopeV[0]);
      result.push_back(hopeV);
    }


    while(true)
    {
      bool updateBool = (SKback - 1)->update(M, d, dlst, dl, dust, du);
      if(updateBool != 0) break;
      --SKback;


      if(SKback - SK <= 1)
      {
        int addSize = result.size() - rstCurrentSize;
        if(addSize > 0) totalSize.fetch_add(addSize);
      }
      return 0; // all the combinations have been tried
    }


    // update totalSize
    {
      int addSize = result.size() - rstCurrentSize;
      if(addSize > 0) totalSize.fetch_add(addSize);
      rstCurrentSize += addSize;
    }


    if(totalSize >= sizeNeeded or (double)std::clock() > endTime)
    {
      break;
    }
  }


  return SKback - SK;
}
*/



template<typename valtype, typename indtype>
inline std::size_t TTTstack(
    indtype LEN, indtype N,
    indtype d, indtype dlst, indtype dl, indtype dust, indtype du,
    valtype ***M, vec<vec<indtype> > &result, int sizeNeeded,
    mPAT<valtype, indtype> *SK, mPAT<valtype, indtype> *SKback,
    bool useBisearchInFindBounds,
    std::atomic<int> &totalSize, double endTime)
{
  if(SKback <= SK) return SKback - SK;
  int rstCurrentSize = result.size();


  valtype **V = M[0];
  if(LEN == 1)
  {
    for(indtype i = 0; i < N; ++i)
    {
      bool allBetween = 1;
      for(indtype k = 0; k < dl and allBetween; ++k)
      {
        if(V[i][k + dlst] < SK->MIN[k]) allBetween = 0;
      }
      for(indtype k = 0; k < du and allBetween; ++k)
      {
        if(V[i][k + dust] > SK->MAX[k]) allBetween = 0;
      }


      if(allBetween)
      {
        result.push_back(vec<indtype>(1, i));
      }
    }
    // update totalSize
    {
      int addSize = result.size() - rstCurrentSize;
      if(addSize > 0) totalSize.fetch_add(addSize);
    }
    return SKback - SK;
  }


  // std::ofstream outfile("proboutput.csv");
  // outfile << "output\n\n\n";
  vec<indtype> hopeV(LEN);
  indtype *hope = &hopeV[0];
  while(true)
  {
    // outfile << "rstCurrentSize = " << rstCurrentSize << ", ";
    // (SKback - 1)->print(d, dl, du, outfile);
    // outfile << "parent printed ___________________________________\n\n";


    SKback->copyParentGene(*(SKback - 1), d, dl, du);


    // SKback->print(d, dl, du, outfile);
    // outfile << "parent copied ___________________________________\n\n";


    indtype boo = SKback->grow(M, d, dlst, dl, dust, du, useBisearchInFindBounds // , &outfile
                              );


    // SKback->print(d, dl, du, outfile);
    // outfile << "child grown ___________________________________\n\n";


    // continue to give birth.
    if(boo == 1)
    {
      *hope = SKback->s;
      ++hope;
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
      bool updateBool = (SKback - 1)->update(M, d, dlst, dl, dust, du);



      // (SKback - 1)->print(d, dl, du, outfile);
      // outfile << "parent updated ___________________________________\n\n";


      if(updateBool != 0)
      {
        *(hope - 1) = (SKback - 1)->s;
        break;
      }


      --hope;
      --SKback;


      if(SKback - SK <= 1)
      {
        // update totalSize
        {
          int addSize = result.size() - rstCurrentSize;
          if(addSize > 0) totalSize.fetch_add(addSize);
        }
        return 0; // all the combinations have been tried
      }
    }


    // update totalSize
    {
      int addSize = result.size() - rstCurrentSize;
      if(addSize > 0) totalSize.fetch_add(addSize);
      rstCurrentSize += addSize;
    }


    if(totalSize >= sizeNeeded or (double)std::clock() > endTime)
    {
      break;
    }
  }


  return SKback - SK;
}



