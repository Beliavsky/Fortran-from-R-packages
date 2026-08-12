# include <RcppParallel.h>
# include "mPATclass.hpp"




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
      if(addSize > 0) totalSize.fetch_and_add(addSize);
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


    indtype boo = SKback->grow(M, d, dlst, dl, dust, du, hope, useBisearchInFindBounds // , &outfile
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
      std::copy(SKback->UB, SKback->UB + SKback->len, hope);
      result.push_back(hopeV);
    }


    while(true)
    {
      bool updateBool = (SKback - 1)->update(M, d, dlst, dl, dust, du);


      // (SKback - 1)->print(d, dl, du, outfile);
      // outfile << "parent updated ___________________________________\n\n";


      if(updateBool != 0) break;


      hope -= (SKback - 1)->Nzeroed;
      --SKback;


      if(SKback - SK <= 1)
      {
        // update totalSize
        {
          int addSize = result.size() - rstCurrentSize;
          if(addSize > 0) totalSize.fetch_and_add(addSize);
        }
        return 0; // all the combinations have been tried
      }
    }


    // update totalSize
    {
      int addSize = result.size() - rstCurrentSize;
      if(addSize > 0) totalSize.fetch_and_add(addSize);
      rstCurrentSize += addSize;
    }


    if(totalSize >= sizeNeeded or (double)std::clock() > endTime)
    {
      break;
    }
  }


  return SKback - SK;
}




/*
// ------------------------------------------------------------------------------------------------
inline int findIdle(std::atomic<unsigned> &totalFreeThreads,
                    std::atomic<unsigned char> *idle0prepare1ready2, unsigned maxCore)
{
  if(totalFreeThreads == 0) return -1;
  for(unsigned i = 0; i < maxCore; ++i)
  {
    bool idle = idle0prepare1ready2[i].compare_and_swap(1, 0);
    if(!idle)
    {
      --totalFreeThreads;
      return i;
    }
  }
  return -1;
}




template<typename valtype, typename indtype>
inline std::size_t TTTstackPar(
    indtype LEN, indtype N, // LEN is the original subset size
    indtype d, indtype dlst, indtype dl, indtype dust, indtype du,
    valtype ***M, vec<vec<indtype> > &result, int sizeNeeded,
    mPAT<valtype, indtype> *SK, mPAT<valtype, indtype> *SKback,
    bool useBisearchInFindBounds,
    std::atomic<int> &totalSize, double endTime,
    vec<indtype> &hopeV, // this is one element in hopeGroup


    // std::atomic<unsigned char> &currentThread, // this is an element in idle0prepare1ready2
    std::atomic<unsigned> &totalFreeThreads,
    std::atomic<unsigned char> *idle0prepare1ready2, unsigned maxCore,
    vec<vec<mPAT<valtype, indtype> > > &SKgroup,
    vec<vec<indtype> > &indGroup,
    vec<vec<valtype> > &valGroup,
    vec<vec<indtype> > &hopeGroup)
{
  if(SKback <= SK)
  {
    return SKback - SK;
  }


  int rstCurrentSize = result.size();


  valtype **V = M[0];
  if(SKback->len == 1)
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
        hopeV.push_back(i);
        result.push_back(hopeV);
        hopeV.pop_back();
      }
    }
    // update totalSize
    {
      int addSize = result.size() - rstCurrentSize;
      if(addSize > 0) totalSize.fetch_and_add(addSize);
    }
    return SKback - SK;
  }


  // std::ofstream outfile("proboutput.csv");
  // outfile << "output\n\n\n";
  indtype *hope = &*hopeV.end();
  hopeV.resize(LEN);
  while(true)
  {
    // outfile << "rstCurrentSize = " << rstCurrentSize << ", ";
    // (SKback - 1)->print(d, dl, du, outfile);
    // outfile << "parent printed ___________________________________\n\n";


    SKback->copyParentGene(*(SKback - 1), d, dl, du);


    // SKback->print(d, dl, du, outfile);
    // outfile << "parent copied ___________________________________\n\n";


    int idle = findIdle(totalFreeThreads, idle0prepare1ready2, maxCore);
    indtype boo = 0;
    if(idle < 0) boo = SKback->grow(M, d, dlst, dl, dust, du, hope, useBisearchInFindBounds //, &outfile
                                    );
    else
    {
      mPAT<valtype, indtype> *idleSK = &SKgroup[idle][1];
      valtype *idleVal = &valGroup[idle][0];
      indtype *idleInd = &indGroup[idle][0];
      vec<indtype> &idleHope = hopeGroup[idle];
      idleHope.resize(hope - &*hopeV.begin());
      std::copy(&*hopeV.begin(), hope, &*idleHope.begin()); // project existing hope to new thread
      std::atomic<unsigned char> &ready = idle0prepare1ready2[idle];
      boo = SKback->growHalfProjectHalf(
        M, d, dlst, dl, dust, du, hope, useBisearchInFindBounds,
        idleSK, idleVal, idleInd, ready);
    }


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
      std::copy(SKback->UB, SKback->UB + SKback->len, hope);
      result.push_back(hopeV);
    }


    while(true)
    {
      bool updateBool = (SKback - 1)->update(M, d, dlst, dl, dust, du);


      // (SKback - 1)->print(d, dl, du, outfile);
      // outfile << "parent updated ___________________________________\n\n";


      if(updateBool != 0) break;


      hope -= (SKback - 1)->Nzeroed;
      --SKback;


      if(SKback - SK <= 1)
      {
        // update totalSize
        {
          int addSize = result.size() - rstCurrentSize;
          if(addSize > 0) totalSize.fetch_and_add(addSize);
        }
        return 0; // all the combinations have been tried
      }
    }


    // update totalSize
    {
      int addSize = result.size() - rstCurrentSize;
      if(addSize > 0) totalSize.fetch_and_add(addSize);
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



