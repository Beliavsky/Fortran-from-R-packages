# pragma once
// # include <Rcpp.h>
// #include <ctime>
// #include <setjmp.h>
# include "oneDoperation.hpp"
# include "mvalOperation.hpp"
# include "macros.hpp"
# include <fstream>
# include <bitset>
// using namespace Rcpp;

// jmp_buf env;




template<typename valtype, typename indtype, bool mk, bool useBiSearch>
// inline unsigned char LBiFind(
//     indtype &ciLB, valtype ***M, indtype ci_1LB, valtype *SR,
//     indtype d, indtype I, indtype &J, indtype *UB, bool useBinarySearch)
inline unsigned char LBiFind(
    indtype &ciLB, valtype ***M, indtype ci_1LB, valtype *SR,
    indtype cmpst, indtype dcmp,
    indtype I, indtype &J, indtype *UB, INT *mask)
{

  if(ciLB < ci_1LB + 1) ciLB = ci_1LB + 1;
  valtype **v = M[0];
  // mvalPlus(SR, SR, v[UB[I]], d);
  mvalPlus(SR + cmpst, SR + cmpst, v[UB[I]] + cmpst, dcmp);


  // std::cout << (int)cmpst << ", " << (int)dcmp << "\n";


  while(true)
  {
    if(UB[J] >= ciLB - (I - J)) break; // when I == J, certainly UB[I] >= ciLB and it breaks
    // mvalMinus(SR, SR, v[UB[J]], d);
    mvalMinus(SR + cmpst, SR + cmpst, v[UB[J]] + cmpst, dcmp);
    ++J;
  }


  indtype ic = 0; // , icPre = ic;


  // the 2nd force is that v[UB[J]]+v[UB[J]+1]+...+v[UB[J]+I-J] must >= SR
  while(true)
  {
    // what is the upper bound of I? UB[I], and the actual value is [UB[I]]
    // icPre = ic;
    // std::cout << (int)I << ", " << (int)J << "\n";
    if(J >= I)
    {
      // if(notAllGreaterEqual(v[UB[I]], SR, ic, d))
      if(notAllGreaterEqual<valtype, indtype, mk> (
          v[UB[I]] + cmpst, SR + cmpst, ic, dcmp, mask))
      {
        // std::cout << "J >= I, return 0\n";
        // std::cout << (int)I << ", " << (int)J << "\n";
        return 0;
      }
      break;
    }


    // which column in M?: M[I-J]
    // what is the value we are looking for? v[UB[J]]+v[UB[J]+1]+...+v[UB[J]+I-J]
    // which row in that column?: M[I-J][UB[J]]
    // if(notAllGreaterEqual(M[I - J][UB[J]], SR, ic, d))
    if(notAllGreaterEqual<valtype, indtype, mk> (
        M[I - J][UB[J]] + cmpst, SR + cmpst, ic, dcmp, mask))
    {
      // mvalMinus(SR, SR, v[UB[J]], d);
      mvalMinus(SR + cmpst, SR + cmpst, v[UB[J]] + cmpst, dcmp);
      ++J;
    }
    else break;
  }


  // the free point J is now located, conduct binary search
  indtype I_J = I - J;
  if(useBiSearch)
  {
    // ciLB = mvalLowerBoundBiMan(&M[I_J][ciLB - I_J], &M[I_J][UB[J]] + 1, SR, d) - &M[I_J][0] + I_J;
    ciLB = mvalLowerBoundBiMan<valtype, indtype, mk> (
      &M[I_J][ciLB - I_J], &M[I_J][UB[J]] + 1, SR, cmpst, dcmp, mask) - &M[I_J][0] + I_J;
  }
  else
  {
    ic = 0;
    // ciLB = mvalLowerBoundLr(&M[I_J][ciLB - I_J], &M[I_J][UB[J]] + 1, SR, ic, d) - &M[I_J][0] + I_J;
    ciLB = mvalLowerBoundLr<valtype, indtype, mk> (
      &M[I_J][ciLB - I_J], &M[I_J][UB[J]] + 1, SR, ic, cmpst, dcmp, mask) - &M[I_J][0] + I_J;
  }


  return 1;
}




template<typename valtype, typename indtype, bool mk, bool useBiSearch>
// inline unsigned char UBiFind(indtype &ciUB, valtype ***M, indtype ciP1UB, valtype *SR,
//                              indtype d, indtype I, indtype &J, indtype *LB, bool useBinarySearch)
inline unsigned char UBiFind(indtype &ciUB, valtype ***M, indtype ciP1UB, valtype *SR,
                             indtype cmpst, indtype dcmp,
                             indtype I, indtype &J, indtype *LB, INT *mask)
{
  if(ciUB > ciP1UB - 1) ciUB = ciP1UB - 1;
  valtype **v = M[0];
  // mvalPlus(SR, SR, v[LB[I]], d);
  mvalPlus(SR + cmpst, SR + cmpst, v[LB[I]] + cmpst, dcmp);


  while(true)
  {
    if(LB[J] <= ciUB + (J - I)) break;
    // mvalMinus(SR, SR, v[LB[J]], d);
    mvalMinus(SR + cmpst, SR + cmpst, v[LB[J]] + cmpst, dcmp);
    --J;
  }


  indtype ic = 0; // icPre = ic;


  // the 2nd force is that v[UB[J]]+v[UB[J]+1]+...+v[UB[J]+I-J] must >= SR
  while(true)
  {
    // icPre = ic;
    if(I == J)
    {
      // if(notAllLessEqual(v[LB[I]], SR, ic, d))
      if(notAllLessEqual<valtype, indtype, mk> (
          v[LB[I]] + cmpst, SR + cmpst, ic, dcmp, mask))
      {
        return 0;
      }
      break;
    }


    // if(notAllLessEqual(M[J - I][LB[J] - (J - I)], SR, ic, d)) // M[J-I][LB[J]-(J-I)]>SR
    if(notAllLessEqual<valtype, indtype, mk> (
        M[J - I][LB[J] - (J - I)] + cmpst, SR + cmpst, ic, dcmp, mask)) // M[J-I][LB[J]-(J-I)] > SR
    {
      // mvalMinus(SR, SR, v[LB[J]], d);
      mvalMinus(SR + cmpst, SR + cmpst, v[LB[J]] + cmpst, dcmp);
      --J;
    }
    else break;
  }


  // the free point J is now located, conduct binary search
  indtype J_I = J - I;
  if(useBiSearch)
  {
    // ciUB = mvalUpperBoundBiMan(&M[J_I][LB[J] - J_I], &M[J_I][ciUB] + 1, SR, d) - 1 - &M[J_I][0];
    ciUB = mvalUpperBoundBiMan<valtype, indtype, mk> (
      &M[J_I][LB[J] - J_I], &M[J_I][ciUB] + 1, SR, cmpst, dcmp, mask) - 1 - &M[J_I][0];
  }
  else
  {
    ic = 0;
    // ciUB = mvalUpperBoundLr(&M[J_I][LB[J] - J_I], &M[J_I][ciUB] + 1, SR, ic, d) - &M[J_I][0];
    ciUB = mvalUpperBoundLr<valtype, indtype, mk> (
      &M[J_I][LB[J] - J_I], &M[J_I][ciUB] + 1, SR, ic, cmpst, dcmp, mask) - &M[J_I][0];
  }


  return 1;
}




// dlend is the end of dimensions bounded from below
// dustart is the start of dimensions bounded from above
// for a multidimensional subset sum, dlend = d and dustart = 0
// LBiFind focuses on [0, dlend)


//_________________________________________________________________________________________________
template<typename valtype, typename indtype, bool mk, bool useBiSearch>
// inline indtype findBoundCpp(
//     indtype len, indtype d, valtype *x, valtype *ME, indtype *LB, valtype *sumLB, indtype *UB, valtype *sumUB,
//     valtype ***M, bool useBinarySearch)
inline indtype findBoundCpp(
    indtype len, indtype d, indtype dlst, indtype dl, indtype dust, indtype du,
    valtype *Min, valtype *Max, // valtype *ME,
    indtype *LB, valtype *sumLB, // int &sumLBind,
    indtype *UB, valtype *sumUB, // int &sumUBind,
    valtype ***M, INT *mask, vec<valtype> &SRVcntr)
// Min is of size dl, Max is of size du
// Min, Min + dl is aligned with value + dlst ~ value + dlst + dl
// Max, Max + du is aligned with value + dust ~ value + dust + du
{
  bool boo = 0;
  valtype **v = M[0];


  // vec<valtype> SRVcntr(d);
  // SRVcntr.resize(d);
  valtype *SRV = &*SRVcntr.begin();
  unsigned LBsum = 0, UBsum = 0;
  while(true)
  {
    bool boundChanged = 0;


    indtype I = 0, J = 0;
    valtype *SR = SRV;


    mvalPlusMinus(SR + dlst, Min + dlst, v[UB[I]] + dlst, sumUB + dlst, dl);


    // Rcpp::Rcout << "LB = ";
    // for(int u = 0; u < len; ++u) Rcpp::Rcout << LB[u] << ", ";
    // Rcpp::Rcout << "\nUB = ";
    // for(int u = 0; u < len; ++u) Rcpp::Rcout << UB[u] << ", ";
    // Rcpp::Rcout << "\n\n";


    // update first lower bound, 0;
    LBsum = 0;
    {
      indtype tmpLB = LB[I];
      if(useBiSearch)
      {
        LB[I] = mvalLowerBoundBiMan<valtype, indtype, mk> (
          &v[LB[I]], &v[UB[I]] + 1, SR, dlst, dl, mask) - v;
      }
      else
      {
        LB[I] = mvalLowerBoundLr<valtype, indtype, mk> (
          &v[LB[I]], &v[UB[I]] + 1, SR, dlst, dl, mask) - v;
      }
      if(LB[I] > UB[I]) return 0;
      boundChanged = boundChanged or (tmpLB != LB[I]);
      std::memcpy(sumLB, v[LB[I]], sizeof(valtype) * d); // every sumLB's dimension matters!
    }
    LBsum += LB[I];


    ++I;


    // update lower bounds from 1 to the end
    {
      for(; I < len; ++I)
      {
        indtype LBI = LB[I];
        unsigned char b = LBiFind<valtype, indtype, mk, useBiSearch> (
          LB[I], M, LB[I - 1], SR, dlst, dl, I, J, UB, mask);
        if(b == 0) return 0;
        boundChanged = boundChanged or (LBI != LB[I]);
        mvalPlus(sumLB, sumLB, v[LB[I]], d); // every sumLB's dimension matters!
        LBsum += LB[I];
      }
    }


    {
      if(!boo) boo = 1; // not the first time go through
      else
      {
        if(!boundChanged)
        {
          if(LBsum == UBsum) return 2;
          break;
        }
      }
    }


    // Rcpp::Rcout << "LB = ";
    // for(int u = 0; u < len; ++u) Rcpp::Rcout << LB[u] << ", ";
    // Rcpp::Rcout << "\nUB = ";
    // for(int u = 0; u < len; ++u) Rcpp::Rcout << UB[u] << ", ";
    // Rcpp::Rcout << "\n\n\n\n";



    boundChanged = 0;
    J = len - 1; I = J;
    mvalPlusMinus(SR + dust, Max, v[LB[I]] + dust, sumLB + dust, du);


    UBsum = 0;
    {
      indtype tmpUB = UB[I];
      if(useBiSearch)
      {
        UB[I] = mvalUpperBoundBiMan<valtype, indtype, mk> (
          &v[LB[I]], &v[UB[I]] + 1, SR, dust, du, mask) - v - 1;
      }
      else
      {
        UB[I] = mvalUpperBoundLr<valtype, indtype, mk> (
          &v[LB[I]], &v[UB[I]] + 1, SR, dust, du, mask) - v;
      }
      if(LB[I] > UB[I]) return 0;
      boundChanged = boundChanged or (tmpUB != UB[I]);
      std::memcpy(sumUB, v[UB[I]], sizeof(valtype) * d);
      UBsum += UB[I];
    }


    --I;


    // refresh the rest UB bounds
    {
      for(; I >= 0; --I)
      {
        indtype UBI = UB[I];
        unsigned char b = UBiFind<valtype, indtype, mk, useBiSearch> (
          UB[I], M, UB[I + 1], SR, dust, du, I, J, LB, mask);
        if(b == 0) return 0;
        boundChanged = boundChanged or (UBI != UB[I]);
        mvalPlus(sumUB, sumUB, v[UB[I]], d); // every dimension of sumUB matters.
        UBsum += UB[I];
      }
    }


    // can I jump out now
    {
      if(!boundChanged)
      {
        if(LBsum == UBsum) return 2;
        break;
      }
    }
  }


  return 1;
}





//_________________________________________________________________________________________________
template<typename valtype, typename indtype, bool mk, bool useBiSearch>
inline indtype findBoundUpFirstCpp(
    indtype len, indtype d, indtype dlst, indtype dl, indtype dust, indtype du,
    valtype *Min, valtype *Max, // valtype *ME,
    indtype *LB, valtype *sumLB, // int &sumLBind,
    indtype *UB, valtype *sumUB, // int &sumUBind,
    valtype ***M, INT *mask, vec<valtype> &SRVcntr)
{
  bool boo = 0;
  valtype **v = M[0];


  unsigned LBsum = 0, UBsum = 0;
  bool boundChanged;
  indtype I, J;
  // vec<valtype> SRVcntr(d);
  valtype *SRV = &*SRVcntr.begin(), *SR = SRV;
  while(true)
  {
    J = len - 1;
    I = J;
    mvalPlusMinus(SR + dust, Max, v[LB[I]] + dust, sumLB + dust, du);
    boundChanged = 0;


    UBsum = 0;
    {
      indtype tmpUB = UB[I];
      if(useBiSearch)
      {
        UB[I] = mvalUpperBoundBiMan<valtype, indtype, mk>(
          &v[LB[I]], &v[UB[I]] + 1, SR, dust, du, mask) - v - 1;
      }
      else
      {
        UB[I] = mvalUpperBoundLr<valtype, indtype, mk>(
          &v[LB[I]], &v[UB[I]] + 1, SR, dust, du, mask) - v;
      }
      if(LB[I] > UB[I]) return 0;
      boundChanged = boundChanged or (tmpUB != UB[I]);
      std::memcpy(sumUB, v[UB[I]], sizeof(valtype) * d);
      UBsum += UB[I];
    }


    --I;


    // refresh the rest UB bounds
    {
      for(; I >= 0; --I)
      {
        indtype UBI = UB[I];
        unsigned char b = UBiFind<valtype, indtype, mk, useBiSearch> (
          UB[I], M, UB[I + 1], SR, dust, du, I, J, LB, mask);
        if(b == 0) return 0;
        boundChanged = boundChanged or (UBI != UB[I]);
        mvalPlus(sumUB, sumUB, v[UB[I]], d); // every dimension of sumUB matters.
        UBsum += UB[I];
      }
    }


    // middle life
    {
      if(!boo) boo = 1; // not the first time go through
      else
      {
        if(!boundChanged)
        {
          if(LBsum == UBsum) return 2;
          break;
        }
      }
    }


    I = 0;
    J = 0;
    mvalPlusMinus(SR + dlst, Min + dlst, v[UB[I]] + dlst, sumUB + dlst, dl);
    boundChanged = 0;


    // update first lower bound, 0;
    LBsum = 0;
    {
      indtype tmpLB = LB[I];
      if(useBiSearch)
      {
        LB[I] = mvalLowerBoundBiMan<valtype, indtype, mk> (
          &v[LB[I]], &v[UB[I]] + 1, SR, dlst, dl, mask) - v;
      }
      else
      {
        LB[I] = mvalLowerBoundLr<valtype, indtype, mk> (
          &v[LB[I]], &v[UB[I]] + 1, SR, dlst, dl, mask) - v;
      }
      if(LB[I] > UB[I]) return 0;
      boundChanged = boundChanged or (tmpLB != LB[I]);
      std::memcpy(sumLB, v[LB[I]], sizeof(valtype) * d); // every sumLB's dimension matters!
    }
    LBsum += LB[I];


    ++I;


    // update lower bounds from 1 to the end
    {
      for(; I < len; ++I)
      {
        indtype LBI = LB[I];
        unsigned char b = LBiFind<valtype, indtype, mk, useBiSearch> (
          LB[I], M, LB[I - 1], SR, dlst, dl, I, J, UB, mask);
        if(b == 0) return 0;
        boundChanged = boundChanged or (LBI != LB[I]);
        mvalPlus(sumLB, sumLB, v[LB[I]], d); // every sumLB's dimension matters!
        LBsum += LB[I];
      }
    }


    // can I jump out now
    {
      if(!boundChanged)
      {
        if(LBsum == UBsum) return 2;
        break;
      }
    }


  }
  return 1;
}



