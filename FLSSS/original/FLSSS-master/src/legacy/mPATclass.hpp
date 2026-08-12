# pragma once
# include <fstream>
# include "triMat.hpp"
# include "macros.hpp"
# include "mvalFindBound.hpp"
# include "raiseSupressBound.hpp"




template<typename valtype, typename indtype>
struct mPAT
{
  indtype position; // position would also be the length of UBleftReserve
  indtype s, send;
  indtype len;


  // members below are stored in another container
  indtype *LB, *UB, *Bresv; // Bresv always needs to be initialized
  valtype *MIN, *MAX; // Min is of size dl, Max is of size du
  valtype *sumLB, *sumUB;


  inline void print(indtype d, indtype dl, indtype du)
  {
    Rcpp::Rcout << "position = " << (int)position << ", s = "<< (int)s <<", send = " << (int)send <<
      ", len = " << (int)len << "\n";
    Rcpp::Rcout << "target LB and UB = ";
    for(int i = 0; i < dl; ++i) Rcpp::Rcout <<  MIN[i] << ",";
    Rcpp::Rcout << ",,";
    for(int i = 0; i < du; ++i) Rcpp::Rcout <<  MAX[i] << ", ";
    Rcpp::Rcout << "\n";


    Rcpp::Rcout << "sumLB = ";
    for(int i = 0, iend = d; i < iend; ++i)
    {
      Rcpp::Rcout << sumLB[i] << ", ";
    }
    Rcpp::Rcout << "\n";


    Rcpp::Rcout << "sumUB = ";
    for(int i = 0, iend = d; i < iend; ++i)
    {
      Rcpp::Rcout << sumUB[i] << ", ";
    }
    Rcpp::Rcout << "\n";


    Rcpp::Rcout << "LB = ";
    for(int i = 0, iend = len; i < iend; ++i)
    {
      Rcpp::Rcout << (int)LB[i] << ", ";
    }
    Rcpp::Rcout<<"\n";


    Rcpp::Rcout<<"UB = ";
    for(int i = 0, iend = len; i < iend; ++i)
    {
      Rcpp::Rcout << (int)UB[i] << ", ";
    }
    Rcpp::Rcout << "\n";


    Rcpp::Rcout << "Bresv = ";
    indtype BresvSize = 0;
    if(position != len) BresvSize = position;

    for(int i = 0, iend = BresvSize; i < iend; ++i)
    {
      Rcpp::Rcout << (int)Bresv[i] <<", ";
    }
    Rcpp::Rcout << "\n\n";
  }


  inline void print(indtype d, indtype dl, indtype du, std::ofstream &outfile)
  {
    outfile << "position =, " << (int)position << ", s =, "<< (int)s <<", send =, " << (int)send <<
      ", len =, " << (int)len << "\n";
    outfile << "target LB and UB = ";
    for(int i = 0; i < dl; ++i) Rcpp::Rcout <<  MIN[i] << ",";
    outfile << ",,";
    for(int i = 0; i < du; ++i) Rcpp::Rcout <<  MAX[i] << ", ";
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
    indtype BresvSize = 0;
    if(position != len) BresvSize = position;
    for(int i = 0, iend = BresvSize; i < iend; ++i)
    {
      outfile << (int)Bresv[i] <<", ";
    }
    outfile << "\n";
  }


  // len is the parent's subset size
  inline void copyParentGene(mPAT &x, indtype d, indtype dl, indtype du) // x is the parent
  {
    len = x.len;


    indtype BresvSize = 0;
    if(x.position != len) BresvSize = x.position;


    MIN = x.sumUB + d;
    MAX = MIN + dl;
    sumLB = MAX + du;
    sumUB = sumLB + d;


    LB = x.Bresv + BresvSize;
    UB = LB + len;


    std::memcpy(MIN, x.MIN, sizeof(valtype) * ((std::size_t)dl + (std::size_t)du + (std::size_t)d * 2));
    std::memcpy(LB, x.LB, sizeof(indtype) * len); // ! do not think x.LB and x.UB are continous !
    std::memcpy(UB, x.UB, sizeof(indtype) * len);
  }


  // equavalent to giveBirth(), and len here is still gene from the parent
  // inline indtype grow(valtype *ME, valtype ***M, indtype d, bool useBiSearch, std::ofstream *outfile = nullptr)
  inline indtype grow(valtype ***M, indtype d, indtype dlst, indtype dl, indtype dust, indtype du,
                      bool useBiSearch, std::ofstream *outfile = nullptr)
  {
    // std::clock_t t = std::clock();
    // indtype boo = findBoundCpp(len, d, target, ME, LB, sumLB, UB, sumUB, M, useBiSearch);
    indtype boo = findBoundCpp(len, d, dlst, dl, dust, du, MIN, MAX, LB, sumLB, UB, sumUB, M, useBiSearch);
    // findBoundTime += std::clock() - t;


    // this->print(d, *outfile);
    // *outfile << "Bounds found ___________________________________\n";


    if(boo == 0) return 0;
    if(len == 1) return 3;
    if(boo == 2) return 2;


    // find the slot that has the least gap
    position = 0;
    indtype Min = *UB - *LB;
    for(indtype i = 1; i < len; ++i)
    {
      indtype tmp = UB[i] - LB[i];
      if(Min > tmp)
      {
        Min = tmp;
        position = i;
      }
    }


    if(position == 0)
    {
      s = LB[position];
      send = UB[position];


      mvalMinus(MIN, MIN, M[0][s] + dlst, dl);
      mvalMinus(MAX, MAX, M[0][s] + dust, du);


      mvalMinus(sumLB, sumLB, M[0][s], d);
      mvalMinus(sumUB, sumUB, M[0][send], d);


      //erase LB's and UB's first elements
      ++LB;
      ++UB;
      --len;
      Bresv = UB + len;
    }
    else if(position == len - 1)
    {
      s = UB[position];
      send = LB[position];


      mvalMinus(MIN, MIN, M[0][s] + dlst, dl);
      mvalMinus(MAX, MAX, M[0][s] + dust, du);


      mvalMinus(sumLB, sumLB, M[0][send], d);
      mvalMinus(sumUB, sumUB, M[0][s], d);


      --len;
      Bresv = UB + len;
    }
    else
    {
      s = LB[position];
      send = UB[position];


      mvalMinus(MIN, MIN, M[0][s] + dlst, dl);
      mvalMinus(MAX, MAX, M[0][s] + dust, du);


      mvalMinus(sumLB, sumLB, M[0][s], d);
      mvalMinus(sumUB, sumUB, M[0][send], d);


      std::copy(LB + position + 1, LB + len, LB + position);
      std::copy(UB + position + 1, UB + len, UB + position);


      --len;
      Bresv = UB + len;
      std::memcpy(Bresv, UB, sizeof(indtype) * position);


      // supress the left part of UB
      indtype cap = s - 1;
      indtype i = position - 1;
      for(; i >= 0; --i, --cap)
      {
        if(UB[i] <= cap) break;
        else
        {
          mvalMinus(sumUB, sumUB, M[0][UB[i]], d);
          UB[i] = cap;
        }
      }
      if(i != position - 1) mvalPlus(sumUB, sumUB, M[position - 1 - i - 1][UB[i + 1]], d);
    }
    return 1;
  }


  inline indtype update(valtype ***M, indtype d, indtype dlst, indtype dl, indtype dust, indtype du)
  {
    if(s == send) return 0;
    mvalPlus(MIN, MIN, M[0][s] + dlst, dl);
    mvalPlus(MAX, MAX, M[0][s] + dust, du);


    if(position == 0)
    {
      ++s;
      mvalMinus(MIN, MIN, M[0][s] + dlst, dl);
      mvalMinus(MAX, MAX, M[0][s] + dust, du);


      indtype bottom = s + 1;
      indtype i = 0;
      for(; i < len; ++i, ++bottom)
      {
        if(LB[i] >= bottom) break;
        else
        {
          LB[i] = bottom;
        }
      }
      if(i != 0) mvalMinusPlus(sumLB, sumLB, M[i - 1][LB[0] - 1], M[i - 1][LB[0]], d);
    }
    else if(position == len)
    {
      --s;
      mvalMinus(MIN, MIN, M[0][s] + dlst, dl);
      mvalMinus(MAX, MAX, M[0][s] + dust, du);


      indtype cap = s - 1;
      indtype i = len - 1;
      for(; i >= 0; --i, --cap)
      {
        if(UB[i] <= cap) break;
        else
        {
          UB[i] = cap;
        }
      }
      if(i != len - 1) mvalMinusPlus(sumUB, sumUB, M[len - 1 - i - 1][UB[i + 1] + 1], M[len - 1 - i - 1][UB[i + 1]], d);
    }
    else
    {
      ++s;
      mvalMinus(MIN, MIN, M[0][s] + dlst, dl);
      mvalMinus(MAX, MAX, M[0][s] + dust, du);


      indtype bottom = s + 1;
      indtype i = position;
      for(; i < len; ++i, ++bottom)
      {
        if(LB[i] >= bottom) break;
        else
        {
          LB[i] = bottom;
        }
      }
      if(i != position) mvalMinusPlus(sumLB, sumLB, M[i - position - 1][LB[position] - 1], M[i - position - 1][LB[position]], d);


      i = position - 1;
      for(; i >= 0; --i)
      {
        if(UB[i] >= Bresv[i]) break;
        else
        {
          ++UB[i];
        }
      }
      if(i != position - 1) mvalMinusPlus(sumUB, sumUB, M[position - 1 - i - 1][UB[i + 1] - 1], M[position - 1 - i - 1][UB[i + 1]], d);
    }
    return 1;
  }
};




/*
template<typename valtype, typename indtype>
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
        if((unsigned)x.position * 2 <= (unsigned)len) parentIndEnd = x.Bresv + x.position + 1;
        else parentIndEnd = x.Bresv + (len - x.position);
      }
      else
      {
        parentValEnd = x.sumBresv + d;
        if((unsigned)x.position * 2 <= (unsigned)len) parentIndEnd = x.Bresv + x.position + 1;
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
    std::memcpy(sumLB, x.sumLB, sizeof(valtype) * d);
    std::memcpy(sumUB, x.sumUB, sizeof(valtype) * d);


    std::memcpy(LB, x.LB, sizeof(indtype) * len); // ! do not think x.LB and x.UB are continous !
    std::memcpy(UB, x.UB, sizeof(indtype) * len);
  }


  // equavalent to giveBirth(), and len here is still gene from the parent
  // inline indtype grow(valtype *ME, valtype ***M, indtype d, bool useBiSearch, std::ofstream *outfile = nullptr)
  inline indtype grow(
      valtype ***M, indtype d, indtype dlst, indtype dl, indtype dust, indtype du, indtype *&hope,
      bool useBiSearch, std::ofstream *outfile = nullptr)
  {
    // std::clock_t t = std::clock();
    indtype boo = findBoundCpp(len, d, dlst, dl, dust, du, MIN, MAX, LB, sumLB, UB, sumUB, M, useBiSearch);
    // findBoundTime += std::clock() - t;
    // return 0;


    // print(d, dl, du, *outfile);
    // *outfile << "Bounds found ___________________________________, boo = " << (int)boo << "\n\n";


    if(boo == 0) return 0;
    if(len == 1) return 3;
    if(boo == 2) return 2;


    // find the slot that has the least gap
    position = 0;
    indtype nonzeroMin = -1;
    indtype overlapPosition[len], *olpend = overlapPosition;
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


    // *outfile << "\n\n";
    // {
    //   *outfile << "overlapPosition =,";
    //   for(int i = 0, iend = olpend - overlapPosition; i < iend; ++i)
    //   {
    //     *outfile << (int)overlapPosition[i] << ",";
    //   }
    //   *outfile << "\n\n";
    // }
    // return 0;


    // erase all positions where LB and UB meet.
    Nzeroed = olpend - overlapPosition;
    if(Nzeroed > 0)
    {
      valtype S[d];
      std::fill(S, S + d, 0);
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
    if((unsigned)position * 2 <= (unsigned)len)
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


    if((unsigned)position * 2 <= (unsigned)len)
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
};
*/



