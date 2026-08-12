# pragma once
# include <fstream>
# include "triMat.hpp"
# include "macros.hpp"
# include "mvalFindBound.hpp"
# include "raiseSupressBound.hpp"
# include <RcppParallel.h>




template<typename valtype, typename indtype, bool mk, bool useBiSearch>
struct mPAT
{
  bool beenUpdated;
  indtype position; // position would also be the length of UBleftReserve
  indtype len;
  indtype Nzeroed;


  // members below are stored in another container
  indtype *LB, *UB, *Bresv; // Bresv always needs to be initialized
  valtype *MIN, *MAX; // Min is of size dl, Max is of size du
  valtype *sumLB, *sumUB, *sumBresv;


  inline void print(indtype d, indtype dl, indtype du, std::ofstream &outfile)
  {
    outfile << "position =, " << (int)position << // ", s =, " << (int)s <<", send =, " << (int)send <<
      ", len =, " << (int)len << ",beenUpdated =," << beenUpdated << "\n";
    outfile << "MIN and MAX =,";
    for(int i = 0; i < dl; ++i) outfile <<  MIN[i] << ",";
    outfile << ",,";
    for(int i = 0; i < du; ++i) outfile <<  MAX[i] << ", ";
    outfile << "\n";


    outfile << "sumLB =, ";
    for(int i = 0, iend = d; i < iend; ++i)
    {
      outfile << sumLB[i] << ", ";
    }
    outfile << "\n";


    outfile << "sumUB =, ";
    for(int i = 0, iend = d; i < iend; ++i)
    {
      outfile << sumUB[i] << ", ";
    }
    outfile << "\n";


    outfile << "sumBresv =, ";
    for(int i = 0, iend = d; i < iend; ++i)
    {
      outfile << sumBresv[i] << ", ";
    }
    outfile << "\n";


    outfile << "LB =, ";
    for(int i = 0, iend = len; i < iend; ++i)
    {
      outfile << (int)LB[i] << ", ";
    }
    outfile << "\n";


    outfile << "UB =, ";
    for(int i = 0, iend = len; i < iend; ++i)
    {
      outfile << (int)UB[i] << ", ";
    }
    outfile << "\n";


    outfile << "Bresv =, ";
    indtype BresvSize = position + 1;
    // if(position != len) BresvSize = position + 1;
    for(int i = 0, iend = BresvSize; i < iend; ++i)
    {
      outfile << (int)Bresv[i] <<", ";
    }
    outfile << "\n";
  }


  // len is the parent's subset size
  inline void copyParentGene(mPAT &x, indtype d, indtype dl, indtype du) // x is the parent
  {
    beenUpdated = 0;
    Nzeroed = 0;
    len = x.len;


    valtype *&parentValEnd = MIN;
    indtype *&parentIndEnd = LB;
    {
      if(x.beenUpdated)
      {
        parentValEnd = x.sumUB + d;
        if(x.position <= len / 2) parentIndEnd = x.Bresv + x.position + 1;
        else parentIndEnd = x.Bresv + (len - x.position);
      }
      else
      {
        parentValEnd = x.sumBresv + d;
        if(x.position <= len / 2) parentIndEnd = x.Bresv + x.position + 1;
        else parentIndEnd = x.Bresv + (len - x.position);
      }
    }


    MAX = MIN + dl;
    sumLB = MAX + du;
    sumUB = sumLB + d;
    sumBresv = sumUB + d;


    UB = LB + len;
    Bresv = UB + len;


    std::memcpy(MIN, x.MIN, sizeof(valtype) * (dl + du));
    // {
    //   {
    //     std::ofstream of("debug.csv", std::ios::app);
    //     std::bitset<64> tmp(*x.MIN);
    //     of << "\ncopying parent genes, parent MIN = " << tmp << "\n";
    //     of.close();
    //   }
    //   {
    //     std::ofstream of("debug.csv", std::ios::app);
    //     std::bitset<64> tmp(*(x.MIN + 1));
    //     of << "copying parent genes, parent MAX = " << tmp << "\n";
    //     of.close();
    //   }
    // }
    std::memcpy(sumLB, x.sumLB, sizeof(valtype) * d);
    std::memcpy(sumUB, x.sumUB, sizeof(valtype) * d);


    std::memcpy(LB, x.LB, sizeof(indtype) * len); // ! do not think x.LB and x.UB are continous !
    std::memcpy(UB, x.UB, sizeof(indtype) * len);
  }


  // equavalent to giveBirth(), and len here is still gene from the parent
  // inline indtype grow(valtype *ME, valtype ***M, indtype d, bool useBiSearch, std::ofstream *outfile = nullptr)
  inline indtype grow(
      valtype ***M, indtype d, indtype dlst, indtype dl, indtype dust, indtype du,
      indtype *&hope, INT *mask, vec<valtype> &SRVcntr, std::ofstream *outfile = nullptr)
  {
    // std::clock_t t = std::clock();
    indtype boo = findBoundCpp<valtype, indtype, mk, useBiSearch> (
      len, d, dlst, dl, dust, du, MIN, MAX, LB, sumLB, UB, sumUB, M, mask, SRVcntr);


    if(outfile != nullptr)
    {
      print(d, dl, du, *outfile);
      *outfile << "Bounds found ___________________________________, boo = " << (int)boo << "\n\n";
    }


    if(boo == 0) return 0;
    if(len == 1) return 3;
    if(boo == 2) return 2;


    // find the slot that has the least gap
    position = 0;
    indtype nonzeroMin = -1;
    // indtype overlapPosition[len], *olpend = overlapPosition;
    vec<indtype> acntr(len);
    indtype *overlapPosition = &*acntr.begin(), *olpend = overlapPosition;
    for(indtype i = 0; i < len; ++i)
    {
      indtype tmp = UB[i] - LB[i];
      if(tmp == 0)
      {
        *hope = UB[i];
        ++hope;
        *olpend = i;
        ++olpend;
      }
      else if(nonzeroMin > tmp or nonzeroMin < 0)
      {
        nonzeroMin = tmp;
        position = i;
      }
    }


    // erase all positions where LB and UB meet.
    Nzeroed = olpend - overlapPosition;
    if(Nzeroed > 0)
    {
      // valtype S[d];
      // std::fill(S, S + d, 0);
      vec<valtype> Scntr(d); valtype *S = &*Scntr.begin();
      *olpend = len;
      for(indtype i = 0; i < Nzeroed; ++i)
      {
        indtype &st = overlapPosition[i], &end = overlapPosition[i + 1];
        mvalPlus(S, S, M[0][UB[st]], d);
        std::copy(LB + st + 1, LB + end, LB + st - i);
        std::copy(UB + st + 1, UB + end, UB + st - i);
      }
      len -= Nzeroed;


      mvalMinus(MIN, MIN, S + dlst, dl); // target changes, so MIN and MAX change
      mvalMinus(MAX, MAX, S + dust, du);
      mvalMinus(sumLB, sumLB, S, d);
      mvalMinus(sumUB, sumUB, S, d);


      // after erasion, position may change. Adjust position
      {
        indtype tmp = 0;
        for(indtype *i = overlapPosition; i < olpend; ++i)
        {
          if(position > *i) ++tmp;
          else break;
        }
        position -= tmp;
      }
    }


    // reserve sumUB
    // halve space at position
    beenUpdated = 0;
    if(position <= len / 2)
    {
      indtype *&UBleftResv = Bresv;
      valtype *&sumUBresv = sumBresv;
      indtype cap = (UB[position] + LB[position]) / 2;
      UBleftResv = UB + len;
      std::copy(UB, UB + position + 1, UBleftResv); // including UB at position, Bresv is of size position + 1!
      std::copy(sumUB, sumUB + d, sumUBresv);
      indtype i = position;
      for(; i >= 0; --i, --cap)
      {
        if(UB[i] <= cap) break;
        mvalMinus(sumUB, sumUB, M[0][UB[i]], d);
        UB[i] = cap;
      }
      mvalPlus(sumUB, sumUB, M[position - i - 1][UB[i + 1]], d);
    }
    else
    {
      indtype *&LBrightResv = Bresv;
      valtype *sumLBresv = sumBresv;
      indtype cap = (UB[position] + LB[position]) / 2 + 1;
      LBrightResv = UB + len;
      std::copy(LB + position, LB + len, LBrightResv);
      std::copy(sumLB, sumLB + d, sumLBresv);
      indtype i = position;
      for(; i < len; ++i, ++cap)
      {
        if(LB[i] >= cap) break;
        mvalMinus(sumLB, sumLB, M[0][LB[i]], d);
        LB[i] = cap;
      }
      mvalPlus(sumLB, sumLB, M[i - position - 1][LB[position]], d);
    }


    return 1;
  }




  inline indtype update(valtype ***M, indtype d, indtype dlst, indtype dl, indtype dust, indtype du)
  {
    if(beenUpdated) return 0;


    if(position <= len / 2)
    {
      indtype cap = UB[position] + 1;
      indtype *&UBleftResv = Bresv;
      valtype *&sumUBresv = sumBresv;
      std::copy(UBleftResv, UBleftResv + position + 1, UB);
      std::swap(sumUB, sumUBresv);
      indtype i = position;
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
      indtype cap = LB[position] - 1;
      indtype *&LBrightResv = Bresv;
      valtype *&sumLBresv = sumBresv;
      std::copy(LBrightResv, LBrightResv + (len - position), LB + position);
      std::swap(sumLB, sumLBresv);
      indtype i = position;
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




  // Special grow() function for knapsack problem
  inline indtype growForKnapsack(
      valtype ***M, indtype d, indtype dlst, indtype dl, indtype dust, indtype du,
      indtype *&hope, INT *mask,
      double *profitVec, double &existingSum, double optimalProfit, vec<valtype> &SRVcntr,
      std::ofstream *outfile = nullptr)
  {
    // std::clock_t t = std::clock();
    indtype boo = findBoundCpp<valtype, indtype, mk, useBiSearch> (
      len, d, dlst, dl, dust, du, MIN, MAX, LB, sumLB, UB, sumUB, M, mask, SRVcntr);


    // See if the sum of upper bounds is less than optimalProfit
    if(optimalProfit > 0 and boo != 0)
    {
      double S = existingSum;
      for(indtype i = 0; i < len; ++i)
      {
        S += profitVec[UB[i]];
      }


      if(S <= optimalProfit)
      {
        return 0;
      }
    }


    if(outfile != nullptr)
    {
      print(d, dl, du, *outfile);
      *outfile << "Bounds found ___________________________________, boo = " << (int)boo << "\n\n";
    }


    if(boo == 0) return 0;
    if(len == 1) return 3;
    if(boo == 2) return 2;


    // find the slot that has the least gap
    position = 0;
    indtype nonzeroMin = -1;
    // indtype overlapPosition[len], *olpend = overlapPosition;
    vec<indtype> acntr(len);
    indtype *overlapPosition = &*acntr.begin(), *olpend = overlapPosition;
    for(indtype i = 0; i < len; ++i)
    {
      indtype tmp = UB[i] - LB[i];
      if(tmp == 0)
      {
        *hope = UB[i];
        existingSum += profitVec[*hope]; // update existing extra dim sum
        ++hope;
        *olpend = i;
        ++olpend;
      }
      else if(nonzeroMin > tmp or nonzeroMin < 0)
      {
        nonzeroMin = tmp;
        position = i;
      }
    }


    // erase all positions where LB and UB meet.
    Nzeroed = olpend - overlapPosition;
    if(Nzeroed > 0)
    {
      // valtype S[d];
      // std::fill(S, S + d, 0);
      vec<valtype> Scntr(d); valtype *S = &*Scntr.begin();
      *olpend = len;
      for(indtype i = 0; i < Nzeroed; ++i)
      {
        indtype &st = overlapPosition[i], &end = overlapPosition[i + 1];
        mvalPlus(S, S, M[0][UB[st]], d);
        std::copy(LB + st + 1, LB + end, LB + st - i);
        std::copy(UB + st + 1, UB + end, UB + st - i);
      }
      len -= Nzeroed;


      mvalMinus(MIN, MIN, S + dlst, dl); // target changes, so MIN and MAX change
      mvalMinus(MAX, MAX, S + dust, du);
      mvalMinus(sumLB, sumLB, S, d);
      mvalMinus(sumUB, sumUB, S, d);


      // after erasion, position may change. Adjust position
      {
        indtype tmp = 0;
        for(indtype *i = overlapPosition; i < olpend; ++i)
        {
          if(position > *i) ++tmp;
          else break;
        }
        position -= tmp;
      }
    }


    // reserve sumUB
    // halve space at position
    beenUpdated = 0;
    if(position <= len / 2)
    {
      indtype *&UBleftResv = Bresv;
      valtype *&sumUBresv = sumBresv;
      indtype cap = (UB[position] + LB[position]) / 2;
      UBleftResv = UB + len;
      std::copy(UB, UB + position + 1, UBleftResv); // including UB at position, Bresv is of size position + 1!
      std::copy(sumUB, sumUB + d, sumUBresv);
      indtype i = position;
      for(; i >= 0; --i, --cap)
      {
        if(UB[i] <= cap) break;
        mvalMinus(sumUB, sumUB, M[0][UB[i]], d);
        UB[i] = cap;
      }
      mvalPlus(sumUB, sumUB, M[position - i - 1][UB[i + 1]], d);
    }
    else
    {
      indtype *&LBrightResv = Bresv;
      valtype *sumLBresv = sumBresv;
      indtype cap = (UB[position] + LB[position]) / 2 + 1;
      LBrightResv = UB + len;
      std::copy(LB + position, LB + len, LBrightResv);
      std::copy(sumLB, sumLB + d, sumLBresv);
      indtype i = position;
      for(; i < len; ++i, ++cap)
      {
        if(LB[i] >= cap) break;
        mvalMinus(sumLB, sumLB, M[0][LB[i]], d);
        LB[i] = cap;
      }
      mvalPlus(sumLB, sumLB, M[i - position - 1][LB[position]], d);
    }


    return 1;
  }
};



