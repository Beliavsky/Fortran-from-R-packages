# pragma once
// # include <Rcpp.h>
// #include <ctime>
// #include <setjmp.h>
# include "oneDoperation.hpp"
# include "mvalOperation.hpp"
# include "macros.hpp"
// using namespace Rcpp;
// jmp_buf env;


// int BBB = 1;


template<typename valtype, typename indtype>
// inline unsigned char LBiFind(
//     indtype &ciLB, valtype ***M, indtype ci_1LB, valtype *SR,
//     indtype d, indtype I, indtype &J, indtype *UB, bool useBinarySearch)
inline unsigned char LBiFind(
    indtype &ciLB, valtype ***M, indtype ci_1LB, valtype *SR,
    indtype cmpst, indtype dcmp,
    indtype I, indtype &J, indtype *UB, bool useBinarySearch)
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
    if(J >= I)
    {
      // if(notAllGreaterEqual(v[UB[I]], SR, ic, d))
      if(notAllGreaterEqual(v[UB[I]] + cmpst, SR + cmpst, ic, dcmp))
      {
        return 0;
      }
      break;
    }


    // which column in M?: M[I-J]
    // what is the value we are looking for? v[UB[J]]+v[UB[J]+1]+...+v[UB[J]+I-J]
    // which row in that column?: M[I-J][UB[J]]
    // if(notAllGreaterEqual(M[I - J][UB[J]], SR, ic, d))
    if(notAllGreaterEqual(M[I - J][UB[J]] + cmpst, SR + cmpst, ic, dcmp))
    {
      // mvalMinus(SR, SR, v[UB[J]], d);
      mvalMinus(SR + cmpst, SR + cmpst, v[UB[J]] + cmpst, dcmp);
      ++J;
    }
    else break;
  }


  // the free point J is now located, conduct binary search
  indtype I_J = I - J;
  if(useBinarySearch)
  {
    // ciLB = mvalLowerBoundBiMan(&M[I_J][ciLB - I_J], &M[I_J][UB[J]] + 1, SR, d) - &M[I_J][0] + I_J;
    ciLB = mvalLowerBoundBiMan(&M[I_J][ciLB - I_J], &M[I_J][UB[J]] + 1, SR, cmpst, dcmp) - &M[I_J][0] + I_J;
  }
  else
  {
    ic = 0;
    // ciLB = mvalLowerBoundLr(&M[I_J][ciLB - I_J], &M[I_J][UB[J]] + 1, SR, ic, d) - &M[I_J][0] + I_J;
    ciLB = mvalLowerBoundLr(&M[I_J][ciLB - I_J], &M[I_J][UB[J]] + 1, SR, ic, cmpst, dcmp) - &M[I_J][0] + I_J;
  }


  return 1;
}




template<typename valtype, typename indtype>
// inline unsigned char UBiFind(indtype &ciUB, valtype ***M, indtype ciP1UB, valtype *SR,
//                              indtype d, indtype I, indtype &J, indtype *LB, bool useBinarySearch)
inline unsigned char UBiFind(indtype &ciUB, valtype ***M, indtype ciP1UB, valtype *SR,
                             indtype cmpst, indtype dcmp,
                             indtype I, indtype &J, indtype *LB, bool useBinarySearch)
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
      if(notAllLessEqual(v[LB[I]] + cmpst, SR + cmpst, ic, dcmp))
      {
        return 0;
      }
      break;
    }


    // if(notAllLessEqual(M[J - I][LB[J] - (J - I)], SR, ic, d)) // M[J-I][LB[J]-(J-I)]>SR
    if(notAllLessEqual(M[J - I][LB[J] - (J - I)] + cmpst, SR + cmpst, ic, dcmp)) // M[J-I][LB[J]-(J-I)] > SR
    {
      // mvalMinus(SR, SR, v[LB[J]], d);
      mvalMinus(SR + cmpst, SR + cmpst, v[LB[J]] + cmpst, dcmp);
      --J;
    }
    else break;
  }


  // the free point J is now located, conduct binary search
  indtype J_I = J - I;
  if(useBinarySearch)
  {
    // ciUB = mvalUpperBoundBiMan(&M[J_I][LB[J] - J_I], &M[J_I][ciUB] + 1, SR, d) - 1 - &M[J_I][0];
    ciUB = mvalUpperBoundBiMan(&M[J_I][LB[J] - J_I], &M[J_I][ciUB] + 1, SR, cmpst, dcmp) - 1 - &M[J_I][0];
  }
  else
  {
    ic = 0;
    // ciUB = mvalUpperBoundLr(&M[J_I][LB[J] - J_I], &M[J_I][ciUB] + 1, SR, ic, d) - &M[J_I][0];
    ciUB = mvalUpperBoundLr(&M[J_I][LB[J] - J_I], &M[J_I][ciUB] + 1, SR, ic, cmpst, dcmp) - &M[J_I][0];
  }


  return 1;
}




// dlend is the end of dimensions bounded from below
// dustart is the start of dimensions bounded from above
// for a multidimensional subset sum, dlend = d and dustart = 0
// LBiFind focuses on [0, dlend)


//_________________________________________________________________________________________________
template<typename valtype, typename indtype>
// inline indtype findBoundCpp(
//     indtype len, indtype d, valtype *x, valtype *ME, indtype *LB, valtype *sumLB, indtype *UB, valtype *sumUB,
//     valtype ***M, bool useBinarySearch)
inline indtype findBoundCpp(
    indtype len, indtype d, indtype dlst, indtype dl, indtype dust, indtype du,
    valtype *Min, valtype *Max, // valtype *ME,
    indtype *LB, valtype *sumLB, // int &sumLBind,
    indtype *UB, valtype *sumUB, // int &sumUBind,
    valtype ***M, bool useBinarySearch)
// Min is of size dl, Max is of size du
// Min, Min + dl is aligned with value + dlst ~ value + dlst + dl
// Max, Max + du is aligned with value + dust ~ value + dust + du
{
  // valtype MAX[d], *Max = MAX, MIN[d], *Min = MIN;
  // {
  //   for(indtype k = 0; k < d; ++k)
  //   {
  //     Max[k] = x[k] + ME[k];
  //     Min[k] = x[k] - ME[k];
  //   }
  // }


  // // test if it's worth finding the bounds
  // {
  //   indtype k = 0;
  //   for(; k < d; ++k)
  //   {
  //     if(absDiff(sumUB[k], sumLB[k]) > relaEps) break;
  //     if(sumUB[k] < Min[k] or sumLB[k] > Max[k]) return 0;
  //   }
  //   if(k == d) return 2;
  //   for(; k < d; ++k)
  //   {
  //     if(sumUB[k] < Min[k] or sumLB[k] > Max[k]) return 0;
  //   }
  // }


  bool boo = 0;
  valtype **v = M[0];


  vec<valtype> SRcontainer(d);
  valtype *SRV = &SRcontainer[0], *SR = SRV;
  unsigned LBsum = 0, UBsum = 0;
  while(true)
  {
    bool boundChanged = 0;


    indtype I = 0, J = 0;
    // valtype SRV[d], *SR = SRV;


    // mvalPlusMinus(SR, Min, v[UB[I]], sumUB, d);
    mvalPlusMinus(SR + dlst, Min + dlst, v[UB[I]] + dlst, sumUB + dlst, dl);


    // update first lower bound, 0;
    LBsum = 0;
    {
      indtype tmpLB = LB[I];
      if(useBinarySearch)
      {
        // LB[I] = mvalLowerBoundBiMan(&v[LB[I]], &v[UB[I]] + 1, SR, d) - v;
        LB[I] = mvalLowerBoundBiMan(&v[LB[I]], &v[UB[I]] + 1, SR, dlst, dl) - v;
      }
      else
      {
        // LB[I] = mvalLowerBoundLr(&v[LB[I]], &v[UB[I]] + 1, SR, d) - v;
        LB[I] = mvalLowerBoundLr(&v[LB[I]], &v[UB[I]] + 1, SR, dlst, dl) - v;
      }
      if(LB[I] > UB[I]) return 0;
      if(tmpLB != LB[I]) boundChanged = 1;
      // std::copy(v[LB[I]], v[LB[I]] + d, sumLB);
      std::memcpy(sumLB, v[LB[I]], sizeof(valtype) * d); // every sumLB's dimension matters!
    }
    LBsum += LB[I];


    ++I;


    // update lower bounds from 1 to the end
    {
      for(; I < len; ++I)
      {
        indtype LBI = LB[I];
        // unsigned char b = LBiFind(LB[I], M, LB[I - 1], SR, d, I, J, UB, useBinarySearch);
        unsigned char b = LBiFind(LB[I], M, LB[I - 1], SR, dlst, dl, I, J, UB, useBinarySearch);
        if(b == 0) return 0;
        if(!boundChanged)
        {
          boundChanged = (LBI != LB[I]);
        }
        mvalPlus(sumLB, sumLB, v[LB[I]], d); // every sumLB's dimension matters!
        LBsum += LB[I];
      }
    }


    // {
    //   if(!boo) boo = 1;
    //   else
    //   {
    //     if(!boundChanged)
    //     {
    //       if(LBsum == UBsum) // double insurance and faster
    //       {
    //         bool LBequalUB = 1;
    //         for(indtype k = 0; k < d; ++k)
    //         {
    //           if(absDiff(sumUB[k], sumLB[k]) > relaEps)
    //           {
    //             LBequalUB = 0;
    //             break;
    //           }
    //         }
    //         if(LBequalUB) return 2;
    //       }
    //       break;
    //     }
    //   }
    // }


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


    J = len - 1; I = J;
    // mvalPlusMinus(SR, Max, v[LB[I]], sumLB, d);
    mvalPlusMinus(SR + dust, Max, v[LB[I]] + dust, sumLB + dust, du);
    boundChanged = 0;


    UBsum = 0;
    {
      indtype tmpUB = UB[I];
      if(useBinarySearch)
      {
        // UB[I] = mvalUpperBoundBiMan(&v[LB[I]], &v[UB[I]] + 1, SR, d) - v - 1;
        UB[I] = mvalUpperBoundBiMan(&v[LB[I]], &v[UB[I]] + 1, SR, dust, du) - v - 1;
      }
      else
      {
        // UB[I] = mvalUpperBoundLr(&v[LB[I]], &v[UB[I]] + 1, SR, d) - v;
        UB[I] = mvalUpperBoundLr(&v[LB[I]], &v[UB[I]] + 1, SR, dust, du) - v;
      }
      if(LB[I] > UB[I]) return 0;
      if(tmpUB != UB[I]) boundChanged = 1;
      std::memcpy(sumUB, v[UB[I]], sizeof(valtype) * d);
      UBsum += UB[I];
    }


    --I;


    // refresh the rest UB bounds
    {
      for(; I >= 0; --I)
      {
        indtype UBI = UB[I];


        // unsigned char b = UBiFind(UB[I], M, UB[I + 1], SR, d, I, J, LB, useBinarySearch);
        unsigned char b = UBiFind(UB[I], M, UB[I + 1], SR, dust, du, I, J, LB, useBinarySearch);


        if(b == 0) return 0;
        if(!boundChanged)
        {
          boundChanged = (UBI != UB[I]);
        }
        mvalPlus(sumUB, sumUB, v[UB[I]], d); // every dimension of sumUB matters.
        UBsum += UB[I];
      }
    }


    // can I jump out now
    {
      if(!boundChanged)
      {
        if(LBsum == UBsum)
        // {
        //   bool LBequalUB = 1;
        //   for(indtype k = 0; k < d; ++k)
        //   {
        //     if(absDiff(sumUB[k], sumLB[k]) > relaEps)
        //     {
        //       LBequalUB = 0;
        //       break;
        //     }
        //   }
        //   if(LBequalUB) return 2;
        // }
          return 2;
        break;
      }
    }


    // {
    //   if(!boundChanged)
    //   {
    //     if(LBsum == UBsum) return 2;
    //     break;
    //   }
    // }
  }


  return 1;
}




//_________________________________________________________________________________________________
template<typename valtype, typename indtype>
// inline indtype findBoundCpp(
//     indtype len, indtype d, valtype *x, valtype *ME, indtype *LB, valtype *sumLB, indtype *UB, valtype *sumUB,
//     valtype ***M, bool useBinarySearch)
inline indtype findBoundUpFirstCpp(
    indtype len, indtype d, indtype dlst, indtype dl, indtype dust, indtype du,
    valtype *Min, valtype *Max, // valtype *ME,
    indtype *LB, valtype *sumLB, // int &sumLBind,
    indtype *UB, valtype *sumUB, // int &sumUBind,
    valtype ***M, bool useBinarySearch)
{
  bool boo = 0;
  valtype **v = M[0];


  unsigned LBsum = 0, UBsum = 0;
  bool boundChanged;
  indtype I, J;
  valtype SRcontainer(d);
  valtype *SRV = &SRcontainer[0], *SR = SRV;
  // valtype SRV[d], *SR = SRV;
  while(true)
  {
    J = len - 1;
    I = J;
    mvalPlusMinus(SR + dust, Max, v[LB[I]] + dust, sumLB + dust, du);
    boundChanged = 0;


    UBsum = 0;
    {
      indtype tmpUB = UB[I];
      if(useBinarySearch)
      {
        UB[I] = mvalUpperBoundBiMan(&v[LB[I]], &v[UB[I]] + 1, SR, dust, du) - v - 1;
      }
      else
      {
        UB[I] = mvalUpperBoundLr(&v[LB[I]], &v[UB[I]] + 1, SR, dust, du) - v;
      }
      if(LB[I] > UB[I]) return 0;
      if(tmpUB != UB[I]) boundChanged = 1;
      std::memcpy(sumUB, v[UB[I]], sizeof(valtype) * d);
      UBsum += UB[I];
    }


    --I;


    // refresh the rest UB bounds
    {
      for(; I >= 0; --I)
      {
        indtype UBI = UB[I];
        unsigned char b = UBiFind(UB[I], M, UB[I + 1], SR, dust, du, I, J, LB, useBinarySearch);
        if(b == 0) return 0;
        if(!boundChanged)
        {
          boundChanged = (UBI != UB[I]);
        }
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
      if(useBinarySearch)
      {
        LB[I] = mvalLowerBoundBiMan(&v[LB[I]], &v[UB[I]] + 1, SR, dlst, dl) - v;
      }
      else
      {
        LB[I] = mvalLowerBoundLr(&v[LB[I]], &v[UB[I]] + 1, SR, dlst, dl) - v;
      }
      if(LB[I] > UB[I]) return 0;
      if(tmpLB != LB[I]) boundChanged = 1;
      std::memcpy(sumLB, v[LB[I]], sizeof(valtype) * d); // every sumLB's dimension matters!
    }
    LBsum += LB[I];


    ++I;


    // update lower bounds from 1 to the end
    {
      for(; I < len; ++I)
      {
        indtype LBI = LB[I];
        unsigned char b = LBiFind(LB[I], M, LB[I - 1], SR, dlst, dl, I, J, UB, useBinarySearch);
        if(b == 0) return 0;
        if(!boundChanged)
        {
          boundChanged = (LBI != LB[I]);
        }
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



