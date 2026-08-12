#pragma once
#include <chrono>
#include "mPATclass.hpp"
#include "convertDataCppR.hpp"
#include "makeKsumIndicatorTables.hpp"


#ifndef TimePoint
#define TimePoint std::chrono::time_point<std::chrono::steady_clock>
#endif


template<typename ing>
struct Shared
{
  ing subsetSize;
  ing N, d;
  int sizeNeed;
  std::atomic<int> totalSize;
  TimePoint endTime;
  uint64_t ***M;
  KsumLookUpTable<ing> *ksumtable;


  Shared()
  {
    subsetSize = 0; 
    N = 0; d = 0; 
    sizeNeed = 0;
    totalSize = 0; 
    endTime = std::chrono::time_point<std::chrono::steady_clock>{};
    M = nullptr; ksumtable = nullptr;
  }


  void reset(ing subsetSize, ing N, ing d, int sizeNeed,
             TimePoint endTime, uint64_t ***M,
             KsumLookUpTable<ing> *ksumtable)
  {
    this->subsetSize = subsetSize;
    this->N = N;
    this->d = d;
    this->sizeNeed = sizeNeed;
    this->M = M;
    totalSize = 0;
    this->endTime = endTime;
    this->ksumtable = ksumtable;
  }


  Shared(ing subsetSize, ing N, ing d, int sizeNeed,
         TimePoint endTime, uint64_t ***M,
         KsumLookUpTable<ing> *ksumtable)
  {
    reset(subsetSize, N, d, sizeNeed, endTime, M, ksumtable);
  }


  Rcpp::RawVector save()
  {
    return copy2rRaw(*this);
  }


  void read(Rcpp::RawVector x)
  {
    std::memcpy((void*)(this), &x[0], x.size());
  }
};


// Return true if time is up.
struct Timer
{
  std::size_t i;
  TimePoint endtime;
  Timer(TimePoint endtime) { i = 0; this->endtime = endtime; }
  bool operator() ()
  {
    ++i;
    if (i % 64 != 0) return false;
    return std::chrono::steady_clock::now() > endtime;
  }
};


template<typename ing>
struct mflsssOBJ
{
  ing *hope; // '*hope' points to the first element to write in 'hopeV'
  Shared<ing> *f;
  vec<ing> hopeV;


  // In this version, use a monolithic container to store everything.
  mPAT<ing> *SKback;
  vec<uint64_t> SKvec;
  vec<uint64_t> SRVcntr;
  vec<vec<ing> > result;


  void swap(mflsssOBJ &X)
  {
    std::swap(X.hope, hope);
    std::swap(X.f, f);
    std::swap(X.hopeV, hopeV);
    std::swap(X.SKvec, SKvec);
    std::swap(X.SKback, SKback);
    std::swap(X.SRVcntr, SRVcntr);
    std::swap(X.result, result);
  }


  void setFirstSK(mPAT<ing> *sk, ing len)
  {
    sk->len = len;
    mPATarrangePtrs<ing> (sk, f->d);
    sk->parent = nullptr;
    sk->beenUpdated = true;
  }


  void initialize ( Shared<ing> *fixedInfo, uint64_t *target, ing *LB, ing *UB )
  {
    f = fixedInfo;
    auto stackLen = uint64_t(f->subsetSize + 2);
    auto biscaleFactor = uint64_t(std::log2(f->N + 1.0 - f->subsetSize) + 3);
    // about '+3': once, there was a failed test case given '+0'


    uint64_t indvecSize = stackLen * (stackLen + 1) / 2 * 3 * biscaleFactor;
    uint64_t valvecSize = stackLen * biscaleFactor * f->d * 4;
    uint64_t SKvecSize  = biscaleFactor * f->subsetSize;
    uint64_t totalBytesUpperBound = 64 +
      SKvecSize * sizeof(mPAT<ing>) +
      SKvecSize * sizeof(uint64_t) * 3 + // This is for bedding 3 uin64_t for every mPAT object.
      indvecSize * sizeof(ing) +
      valvecSize * sizeof(uint64_t);
    SKvec.resize((totalBytesUpperBound + sizeof(uint64_t) - 1) / sizeof(uint64_t), 0);
    // Rcpp::Rcout << "SKvec.size() = " << SKvec.size() << "\n\n";


    hopeV.assign(f->subsetSize, 0);
    SRVcntr.assign(f->d + f->subsetSize, 0);
    hope = &hopeV[0];
    mPAT<ing> *SKbegin = (mPAT<ing>*)(&SKvec.front());


    setFirstSK(SKbegin, f->subsetSize);
    SKback = SKbegin->createChildPtr(); // SKbegin's child is set.
    SKback->parent = SKbegin;


    for(ing i = 0; i < SKbegin->len; ++i)
    {
      SKbegin->LB[i] = LB[i];
      SKbegin->UB[i] = UB[i];
    }


    std::memcpy(SKbegin->target, target, sizeof(uint64_t) * f->d);
    iterSum<ing> (SKbegin->sumLB, f->M[0], SKbegin->LB, SKbegin->len, f->d);
    iterSum<ing> (SKbegin->sumUB, f->M[0], SKbegin->UB, SKbegin->len, f->d);
  }


  int64_t TTTstackRun()
  {
    auto SK = SKback->parent;
    int rstCurrentSize = result.size();


    if ( SK->len == 1 ) // Will return in this branch.
    {
      uint64_t **V = f->M[0];
      for (ing i = SK->LB[0], iend = SK->UB[0]; i <= iend; ++i)
      {
        if ( equal(V[i], SK->target, f->d) )
        {
          *hope = i;
          result.push_back(hopeV); // hopeV NEEDS TO BE COPIED!
        }
      }
      // Update totalSize.
      {
        int addSize = result.size() - rstCurrentSize;
        if ( addSize > 0 ) f->totalSize.fetch_add(addSize);
      }
      return SK - (mPAT<ing>*)&SKvec[0];
    }


    auto Skchild = SK->createChildPtr();
    Timer timer(f->endTime);
    while (true)
    {
      SKback->copyParentGene ( f->d );
      ing boo = SKback->grow ( f->M, f->d, hope, SRVcntr, f->ksumtable);


      // Continue to give birth.
      if (boo == 1)
      {
        auto child = SKback->createChildPtr();
        child->parent = SKback;
        SKback = child;
        continue;
      }


      if (boo == 3) // If len (subset size) in the child becomes 1:
      {
        ing i = SKback->LB[0], iend = SKback->UB[0];
        for(; i <= iend; ++i)
        {
          hopeV.back() = i;
          result.push_back(hopeV);
        }
      }
      else if (boo == 2) // If lower bounds and upper bounds overlap
      {
        std::copy(SKback->UB, SKback->UB + SKback->len, hope);


        // std::cout << "target = " << std::flush;
        // for (int u = 0, uend = f->d; u < uend; ++u)
        //   std::cout << SKback->target[u] << ", ";
        // std::cout << std::endl;


        result.push_back(hopeV);
      }


      while(true)
      {
        ing updateBool = SKback->parent->update( f->M, f->d );
        if ( updateBool != 0 ) break;


        hope -= SKback->parent->Nzeroed;
        SKback = SKback->parent;


        if (SKback <= Skchild)
        {
          int addSize = result.size() - rstCurrentSize;
          if (addSize > 0) f->totalSize.fetch_add(addSize);
          return 0; // All the combinations have been tried.
        }
      }


      // Update totalSize
      {
        int addSize = result.size() - rstCurrentSize;
        if(addSize > 0) f->totalSize.fetch_add(addSize);
        rstCurrentSize += addSize;
      }


      if ( f->totalSize >= f->sizeNeed or timer() ) break;
    }


    return SKback - SK;
  }


  Rcpp::List save()
  {
    auto hopeV_ = copyVec2rRaw(hopeV);
    auto tmp2 = ptrInt(hope) - ptrInt(&hopeV[0]);
    auto hopeOffset = copy2rRaw(tmp2);
    auto tmp0 = ptrInt(&SKvec[0]);
    auto SKvec0addr = copy2rRaw(tmp0);
    auto tmp1 = ptrInt(SKback) - ptrInt(&SKvec[0]);
    auto SKbackOffset = copy2rRaw(tmp1);
    auto SKvec_ = copyVec2rRaw(SKvec);
    int SRVcntrCapacity = SRVcntr.capacity();


    return Rcpp::List::create(
      Rcpp::Named("hopeOffset") = hopeOffset,
      Rcpp::Named("hopeV") = hopeV_,
      Rcpp::Named("SKvec0addr") = SKvec0addr,
      Rcpp::Named("SKvec") = SKvec_,
      Rcpp::Named("SKbackOffset") = SKbackOffset,
      Rcpp::Named("SRVcntrCapacity") = SRVcntrCapacity);
  }


  void read(Rcpp::List X, Shared<ing> *f)
  {
    this->f = f;
    typedef Rcpp::RawVector raw;
    raw x1 = X["hopeV"]; copyRraw2vec(hopeV, x1);
    raw x2 = X["SKvec"]; copyRraw2vec(SKvec, x2);
    raw x3 = X["hopeOffset"];
    auto hopeOffset = copyRraw<std::uintptr_t> (x3);
    hope = &hopeV[0] + hopeOffset / sizeof(ing);
    raw x4 = X["SKvec0addr"];
    auto SKvec0addr = copyRraw<std::uintptr_t> (x4);
    raw x5 = X["SKbackOffset"];
    auto SKbackOffset = copyRraw<std::uintptr_t> (x5);
    SKback = intPtr<mPAT<ing>*> (ptrInt(&SKvec[0]) + SKbackOffset);
    int SRVcntrCapacity = X["SRVcntrCapacity"];
    SRVcntr.resize(SRVcntrCapacity);


    for (auto i = SKback; ; i = i->parent)
    {
      i->target = intPtr<uint64_t*> (ptrInt(i->target) - SKvec0addr + ptrInt(&SKvec[0]));
      i->sumLB = intPtr<uint64_t*> (ptrInt(i->sumLB) - SKvec0addr + ptrInt(&SKvec[0]));
      i->sumUB = intPtr<uint64_t*> (ptrInt(i->sumUB) - SKvec0addr + ptrInt(&SKvec[0]));
      i->sumBresv = intPtr<uint64_t*> (ptrInt(i->sumBresv) - SKvec0addr + ptrInt(&SKvec[0]));
      i->LB = intPtr<ing*> (ptrInt(i->LB) - SKvec0addr + ptrInt(&SKvec[0]));
      i->UB = intPtr<ing*> (ptrInt(i->UB) - SKvec0addr + ptrInt(&SKvec[0]));
      i->Bresv = intPtr<ing*> (ptrInt(i->Bresv) - SKvec0addr + ptrInt(&SKvec[0]));
      if (i->parent == 0) break;
      i->parent = intPtr<mPAT<ing>*> (ptrInt(i->parent) - SKvec0addr + ptrInt(&SKvec[0]));
    }
  }




};




template<typename ing>
struct GrowTwin
{
  vec<ing> acntr;
  vec<uint64_t> cntrS;


  ing operator() ( mflsssOBJ<ing> &Xmflsss, mflsssOBJ<ing> &Ymflsss )
  {
    Shared<ing> &f = *Xmflsss.f;
    mPAT<ing> &X = *Xmflsss.SKback;
    X.copyParentGene ( f.d );


    ing boo = findBoundCpp<ing> (
      X.len, f.d, X.target, X.LB, X.sumLB, X.UB, X.sumUB, f.M,
      &Xmflsss.SRVcntr[0], f.ksumtable);
    // Rcpp::Rcout << "f.ksumtable = " << uint64_t(f.ksumtable) << "\n";


    if(boo == 0) return 0;
    if(X.len == 1) return 3;
    if(boo == 2) return 2;


    // Find the slot that has the least gap
    X.position = 0;
    ing nonzeroMin = -1;


    acntr.resize(X.len);
    ing *overlapPosition = &*acntr.begin(), *olpend = overlapPosition;


    for (ing i = 0; i < X.len; ++i)
    {
      ing tmp = X.UB[i] - X.LB[i];
      if(tmp == 0)
      {
        *Xmflsss.hope = X.UB[i];
        ++Xmflsss.hope;
        *olpend = i;
        ++olpend;
      }
      else if(nonzeroMin > tmp or nonzeroMin < 0)
      {
        nonzeroMin = tmp;
        X.position = i;
      }
    }


    // Erase all positions where LB and UB meet.
    ing &Nzeroed = X.Nzeroed;
    Nzeroed = olpend - overlapPosition;
    if (Nzeroed > 0)
    {
      cntrS.assign(f.d, 0);
      uint64_t *S = &*cntrS.begin();
      *olpend = X.len;
      for (ing i = 0; i < Nzeroed; ++i)
      {
        ing &st = overlapPosition[i], &end = overlapPosition[i + 1];
        mvalPlus(S, S, f.M[0][X.UB[st]], f.d);
        std::copy(X.LB + st + 1, X.LB + end, X.LB + st - i);
        std::copy(X.UB + st + 1, X.UB + end, X.UB + st - i);
      }
      X.len -= Nzeroed;


      mvalMinus(X.target, X.target, S, f.d);
      mvalMinus(X.sumLB,  X.sumLB,  S, f.d);
      mvalMinus(X.sumUB,  X.sumUB,  S, f.d);


      // After erasion, position may change. Adjust position.
      {
        ing tmp = 0;
        for(ing *i = overlapPosition; i < olpend; ++i)
        {
          if(X.position > *i) ++tmp;
          else break;
        }
        X.position -= tmp;
      }
    }


    // x, pass wisdom to your twin!
    {
      Ymflsss.f = &f;
      Ymflsss.hopeV.assign(f.subsetSize, 0);
      std::copy(&Xmflsss.hopeV[0], Xmflsss.hope, &Ymflsss.hopeV[0]);
      Ymflsss.hope = &Ymflsss.hopeV[0] + (Xmflsss.hope - &Xmflsss.hopeV[0]);
      Ymflsss.SKvec.resize(Xmflsss.SKvec.size());
      Ymflsss.SRVcntr.assign(f.d + f.subsetSize, 0);
    }


    auto YmflsssSK = (mPAT<ing>*)&Ymflsss.SKvec[0];
    Ymflsss.setFirstSK(YmflsssSK, Xmflsss.SKback->len);
    Ymflsss.SKback = YmflsssSK->createChildPtr();
    Ymflsss.SKback->parent = YmflsssSK;
    mPAT<ing> &Y = *(Ymflsss.SKback->parent);
    X.beenUpdated = true;
    Y.beenUpdated = true;
    std::copy(X.target, X.target + f.d, Y.target);


    std::copy(X.sumUB, X.sumUB + f.d, Y.sumUB);
    std::copy(X.UB, X.UB + X.len, Y.UB);
    ing cap = (X.UB[X.position] + X.LB[X.position]) / 2;
    ing capResv = cap;
    ing i = X.position;
    for(; i >= 0; --i, --cap)
    {
      if(X.UB[i] <= cap) break;
      mvalMinus(X.sumUB, X.sumUB, f.M[0][X.UB[i]], f.d);
      X.UB[i] = cap;
    }
    mvalPlus(X.sumUB, X.sumUB, f.M[X.position - i - 1][X.UB[i + 1]], f.d);


    cap = capResv;
    ++cap;
    i = X.position;
    std::copy(X.LB, X.LB + i, Y.LB);
    std::copy(X.sumLB, X.sumLB + f.d, Y.sumLB);
    for(; i < X.len; ++i, ++cap)
    {
      if(X.LB[i] >= cap)
      {
        std::copy(X.LB + i, X.LB + X.len, Y.LB + i);
        break;
      }
      mvalMinus(Y.sumLB, Y.sumLB, f.M[0][X.LB[i]], f.d);
      Y.LB[i] = cap;
    }
    mvalPlus(Y.sumLB, Y.sumLB, f.M[i - X.position - 1][Y.LB[X.position]], f.d);


    auto child = Xmflsss.SKback->createChildPtr();
    child->parent = Xmflsss.SKback;
    Xmflsss.SKback = child;
    return 1;
  }
};




template<typename ing>
struct Mitosis
{
  vec<unsigned char> acntr;
  GrowTwin<ing> growTwin;


  void operator() (
      vec<mflsssOBJ<ing> > &descendants, Shared<ing> &f,
      vec<vec<ing> > &rstCollection, ing *LB, ing *UB,
      uint64_t *target, int Ndescendants)
  {

    Ndescendants = Ndescendants <= 1 ? 1 :
      1 << int(std::round(std::log2(Ndescendants + 0.0)) + 1);
    descendants.resize(Ndescendants);
    descendants[0].initialize(&f, target, LB, UB);


    acntr.assign(Ndescendants, 0);
    unsigned char *dead = &*acntr.begin();
    int j = 1;


    while (j < Ndescendants)
    {
      int iend = j;
      for (int i = 0; i < iend; ++i, ++j)
      {
        // if (f.totalSize >= f.sizeNeed) return;


        if (dead[i])
        {
          dead[j] = 1;
          continue;
        }


        ing boo = growTwin(descendants[i], descendants[j]);


        mPAT<ing> &tmp = *descendants[i].SKback;
        if(boo == 0)
        {
          dead[i] = 1;
          dead[j] = 1;
        }
        else if(boo == 3)
        {
          ing k = tmp.LB[0], kend = tmp.UB[0];
          for(; k <= kend; ++k)
          {
            descendants[i].hopeV.back() = k;
            rstCollection.push_back(descendants[i].hopeV);
            ++f.totalSize;
          }
          dead[i] = 1;
          dead[j] = 1;
        }
        else if(boo == 2)
        {
          std::copy(tmp.UB, tmp.UB + tmp.len, descendants[i].hope);
          rstCollection.push_back(descendants[i].hopeV);
          ++f.totalSize;
          dead[i] = 1;
          dead[j] = 1;
        }
      }
    }


    // if (f.totalSize >= f.sizeNeed) return;


    int validDescendents = Ndescendants -
      std::accumulate(dead, dead + Ndescendants, (int)0);


    if (validDescendents >= Ndescendants) return;


    vec<mflsssOBJ<ing> > descendantsRemain(validDescendents);
    for (int i = 0, k = 0; i < Ndescendants; ++i)
    {
      if (dead[i]) continue;
      descendants[i].swap(descendantsRemain[k]);
      ++k;
    }
    descendants.swap(descendantsRemain);
  }


};














