#pragma once
#include "triMat.hpp"
#include "mvalFindBound.hpp"
#include "raiseSupressBound.hpp"
#include "makeKsumIndicatorTables.hpp"


template<typename ing>
struct mPAT
{
  bool beenUpdated;
  ing position; // position would also be the length of UBleftReserve
  ing len;
  ing Nzeroed;
  mPAT<ing> *parent;


  // Members below are stored in another container
  uint64_t *target, *sumLB, *sumUB, *sumBresv;
  // sumBresv could be swapped with sumLB or sumUB in update() !
  ing *LB, *UB, *Bresv; // Bresv always needs to be initialized.


  mPAT<ing> *createChildPtr()
  {
    mPAT<ing> *child;
    if (beenUpdated) child = (mPAT<ing>*)N8BAA(Bresv); // Bresv is unnecessary now.
    else
    {
      ing *end = position <= (len - 1) / 2 ?
        Bresv + position + 1 : Bresv + (len - position);
      child = (mPAT<ing>*)N8BAA(end);
    }
    return child;
  }


  // len is the parent's subset size
  void copyParentGene(ing d) // *parent must have already been set!
  {
    mPAT<ing> &x = *parent;
    beenUpdated = false;
    Nzeroed = 0;
    len = x.len;
    mPATarrangePtrs(this, d);
    std::memcpy( target, x.target, sizeof(uint64_t) * d );
    std::memcpy( sumLB,  x.sumLB,  sizeof(uint64_t) * d );
    std::memcpy( sumUB,  x.sumUB,  sizeof(uint64_t) * d );


    // ! Do not think x.LB and x.UB are continous !
    std::memcpy( LB, x.LB, sizeof(ing) * len );
    std::memcpy( UB, x.UB, sizeof(ing) * len );
  }


  // Equivalent to giveBirth(), and len here is still gene from the parent.
  // SRVcntr has sufficient size during mFLSSS object initialization.
  ing grow ( uint64_t ***M, ing d, ing *&hope, vec<uint64_t> &SRVcntr,
             KsumLookUpTable<ing> *ksumtable)
  {
    ing boo = findBoundCpp<ing> (
      len, d, target, LB, sumLB, UB, sumUB, M, &SRVcntr[0], ksumtable);


    if (boo == 0) return 0;
    if (len == 1) return 3;
    if (boo == 2) return 2;


    // Find the slot that has the least gap
    position = 0;
    ing nonzeroMin = -1;
    ing *overlapPosition = (ing*)(&SRVcntr[0]), *olpend = overlapPosition;
    for (ing i = 0; i < len; ++i)
    {
      ing tmp = UB[i] - LB[i];
      if (tmp == 0)
      {
        *hope = UB[i];
        ++hope;
        *olpend = i;
        ++olpend;
      }
      else if ( nonzeroMin > tmp or nonzeroMin < 0 )
      {
        nonzeroMin = tmp;
        position = i;
      }
    }


    // Erase all positions where LB and UB meet.
    Nzeroed = olpend - overlapPosition;
    if (Nzeroed > 0)
    {
      uint64_t *S = &SRVcntr[0] + len; // The first len entries are occupied.
      std::fill(S, S + d, 0);
      *olpend = len;
      for(ing i = 0; i < Nzeroed; ++i)
      {
        ing &st = overlapPosition[i], &end = overlapPosition[i + 1];
        mvalPlus(S, S, M[0][UB[st]], d);
        std::copy(LB + st + 1, LB + end, LB + st - i);
        std::copy(UB + st + 1, UB + end, UB + st - i);
      }
      len -= Nzeroed;


      mvalMinus<ing> (target, target, S, d);
      mvalMinus<ing> (sumLB,  sumLB,  S, d);
      mvalMinus<ing> (sumUB,  sumUB,  S, d);


      // After erasion, position may change. Adjust position
      {
        ing tmp = 0;
        for(ing *i = overlapPosition; i < olpend; ++i)
        {
          if (position > *i) ++tmp;
          else break;
        }
        position -= tmp;
      }
    }


    // Reserve sumUB
    // Halve space at position
    beenUpdated = false;
    if (position <= (len - 1) / 2)
    {
      ing *&UBleftResv = Bresv;
      uint64_t *&sumUBresv = sumBresv;
      ing cap = (UB[position] + LB[position]) / 2;
      std::copy(UB, UB + position + 1, UBleftResv); // Including UB at position, Bresv is of size position + 1!
      std::copy(sumUB, sumUB + d, sumUBresv);
      ing i = position;
      for (; i >= 0; --i, --cap)
      {
        if (UB[i] <= cap) break;
        mvalMinus(sumUB, sumUB, M[0][UB[i]], d);
        UB[i] = cap;
      }
      mvalPlus(sumUB, sumUB, M[position - i - 1][UB[i + 1]], d);
    }
    else
    {
      ing *&LBrightResv = Bresv;
      uint64_t *sumLBresv = sumBresv;
      ing cap = (UB[position] + LB[position]) / 2 + 1;
      std::copy(LB + position, LB + len, LBrightResv);
      std::copy(sumLB, sumLB + d, sumLBresv);
      ing i = position;
      for (; i < len; ++i, ++cap)
      {
        if(LB[i] >= cap) break;
        mvalMinus(sumLB, sumLB, M[0][LB[i]], d);
        LB[i] = cap;
      }
      mvalPlus(sumLB, sumLB, M[i - position - 1][LB[position]], d);
    }


    return 1;
  }


  ing update(uint64_t ***M, ing d)
  {
    if ( beenUpdated ) return 0;


    if (position <= (len - 1) / 2)
    {
      ing cap = UB[position] + 1;
      ing *&UBleftResv = Bresv;
      uint64_t *&sumUBresv = sumBresv;
      std::copy(UBleftResv, UBleftResv + position + 1, UB);
      std::swap(sumUB, sumUBresv);
      ing i = position;
      for(; i < len; ++i, ++cap)
      {
        if(LB[i] >= cap) break;
        mvalMinus(sumLB, sumLB, M[0][LB[i]], d);
        LB[i] = cap;
      }
      mvalPlus(sumLB, sumLB, M[i - position - 1][LB[position]], d);
    }
    else
    {
      ing cap = LB[position] - 1;
      ing *&LBrightResv = Bresv;
      uint64_t *&sumLBresv = sumBresv;
      std::copy(LBrightResv, LBrightResv + (len - position), LB + position);
      std::swap(sumLB, sumLBresv);
      ing i = position;
      for(; i >= 0; --i, --cap)
      {
        if(UB[i] <= cap) break;
        mvalMinus(sumUB, sumUB, M[0][UB[i]], d);
        UB[i] = cap;
      }
      mvalPlus(sumUB, sumUB, M[position - i - 1][UB[i + 1]], d);
    }


    beenUpdated = true;
    return 1;
  }


};




template<typename ing>
inline void mPATarrangePtrs(mPAT<ing> *x, ing d)
{
  x->target = (uint64_t*)N8BAA(x + 1);
  x->sumLB = x->target + d;
  x->sumUB = x->sumLB + d;
  x->sumBresv = x->sumUB + d;
  x->LB = (ing*)N8BAA(x->sumBresv + d);
  x->UB = x->LB + x->len;
  x->Bresv = x->UB + x->len;
}






