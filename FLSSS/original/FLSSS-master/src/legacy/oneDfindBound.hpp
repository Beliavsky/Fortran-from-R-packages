# pragma once
// #include <ctime>
// #include <setjmp.h>
# include "oneDoperation.hpp"
# include "macros.hpp"
// jmp_buf env;


namespace legacy
{
template<typename valtype, typename indtype>
inline unsigned char LBiFind(
    indtype &ciLB, valtype **M, indtype ci_1LB, valtype &SR, indtype I, indtype &J, indtype *UB, bool useBinarySearch)
{
  if(ciLB < ci_1LB + 1) ciLB = ci_1LB + 1;
  valtype *v = M[0];
  SR += v[UB[I]];


  while(true)
  {
    if(UB[J] >= ciLB - (I - J)) break; // when I == J, certainly UB[I] >= ciLB and it breaks
    SR -= v[UB[J]];
    ++J;
  }


  // the 2nd force is that v[UB[J]]+v[UB[J]+1]+...+v[UB[J]+I-J] must >= SR
  while(true)
  {
    // what is the upper bound of I? UB[I], and the actual value is [UB[I]]
    if(J >= I)
    {
      if(v[UB[I]] < SR) return 0; // solution won't exist
      break;
    }


    // which column in M?: M[I-J]
    // what is the value we are looking for? v[UB[J]]+v[UB[J]+1]+...+v[UB[J]+I-J]
    // which row in that column?: M[I-J][UB[J]]
    if(M[I - J][UB[J]] < SR)
    {
      SR -= v[UB[J]];
      ++J;
    }
    else break;
  }


  // the free point J is now located, conduct binary search
  indtype I_J = I - J;
  if(useBinarySearch)
  {
    ciLB = lowerBoundBiMan<valtype, indtype> (&M[I_J][ciLB - I_J], &M[I_J][UB[J]] + 1, SR) - &M[I_J][0] + I_J;
  }
  else
  {
    ciLB = lowerBoundLr<valtype, indtype> (&M[I_J][ciLB - I_J], &M[I_J][UB[J]] + 1, SR) - &M[I_J][0] + I_J;
  }


  return 1;
}




template<typename valtype, typename indtype>
inline unsigned char UBiFind(
    indtype &ciUB, valtype **M, indtype ciP1UB, valtype &SR, indtype I, indtype &J, indtype *LB, bool useBinarySearch)
{
  if(ciUB > ciP1UB - 1) ciUB = ciP1UB - 1;
  valtype *v = M[0];
  SR += v[LB[I]];


  while(true)
  {
    if(LB[J] <= ciUB + (J - I)) break;
    SR -= v[LB[J]];
    --J;
  }


  // the 2nd force is that v[UB[J]]+v[UB[J]+1]+...+v[UB[J]+I-J] must >= SR
  while(true)
  {
    if(I == J)
    {
      if(v[LB[I]] > SR) return 0;
      break;
    }


    if(M[J - I][LB[J] - (J - I)] > SR)
    {
      SR -= v[LB[J]];
      --J;
    }
    else break;
  }


  // the free point J is now located, conduct binary search
  indtype J_I = J - I;
  if(useBinarySearch)
  {
    ciUB = upperBoundBiMan<valtype, indtype> (&M[J_I][LB[J] - J_I], &M[J_I][ciUB] + 1, SR) - 1 - &M[J_I][0];
  }
  else
  {
    ciUB = upperBoundLr<valtype, indtype> (&M[J_I][LB[J] - J_I], &M[J_I][ciUB] + 1, SR) - &M[J_I][0];
  }


  return 1;
}








//_________________________________________________________________________________________________
template<typename valtype, typename indtype>
inline indtype findBoundCpp(indtype len, valtype x, valtype ME,
                            indtype *LB, valtype &sumLB,
                            indtype *UB, valtype &sumUB,
                            valtype **M, bool useBinarySearch)
{
  valtype Min = x - ME, Max = x + ME;
  if(sumUB < Min or sumLB > Max) return 0;
  // if(std::abs(sumUB - sumLB) < eps) return 2;
  if(absDiff(sumLB, sumUB) < relaEps) return 2;


  bool boo = 0;
  valtype *v = M[0];


  while(true)
  {
    bool boundChanged = 0;


    indtype I = 0, J = 0;
    valtype SR = Min + v[UB[I]] - sumUB;


    //update first lower bound, 0;
    {
      indtype tmpLB = LB[I];
      if(useBinarySearch)
      {
        LB[I] = lowerBoundBiMan<valtype, indtype> (&v[LB[I]], &v[UB[I]] + 1, SR) - v;
      }
      else
      {
        LB[I] = lowerBoundLr<valtype, indtype> (&v[LB[I]], &v[UB[I]] + 1, SR) - v;
      }
      if(LB[I] > UB[I]) return 0;
      if(tmpLB != LB[I]) boundChanged = 1;
      sumLB = v[LB[I]];
    }


    ++I;


    //update lower bounds from 1 to the end
    {
      for(; I < len; ++I)
      {
        indtype LBI = LB[I];
        unsigned char b = LBiFind<valtype, indtype> (LB[I], M, LB[I - 1], SR, I, J, UB, useBinarySearch);
        if(b == 0) return 0;
        if(!boundChanged) boundChanged = (LBI != LB[I]);
        sumLB += v[LB[I]];
      }
    }


    // see if can jump out now
    {
      if(!boo) boo = 1;
      else
      {
        if(!boundChanged)
        {
          // if(std::abs(sumLB - sumUB) < eps) return 2;
          if(absDiff(sumLB, sumUB) < relaEps) return 2;
          break;
        }
      }
    }


    J = len - 1; I = J;
    SR = Max + v[LB[I]] - sumLB;
    boundChanged = 0;


    {
      indtype tmpUB = UB[I];
      if(useBinarySearch)
      {
        UB[I] = upperBoundBiMan<valtype, indtype> (&v[LB[I]], &v[UB[I]] + 1, SR) - v - 1;
      }
      else
      {
        UB[I] = upperBoundLr<valtype, indtype> (&v[LB[I]], &v[UB[I]] + 1, SR) - v;
      }
      if(LB[I] > UB[I]) return 0;
      if(tmpUB != UB[I]) boundChanged = 1;
      sumUB = v[UB[I]];
    }


    --I;


    // refresh the rest UB bounds
    {
      for(; I >= 0; --I)
      {
        indtype UBI = UB[I];
        unsigned char b = UBiFind<valtype, indtype> (UB[I], M, UB[I + 1], SR, I, J, LB, useBinarySearch);
        if(b == 0) return 0;
        if(!boundChanged) boundChanged = (UBI != UB[I]);
        sumUB += v[UB[I]];
      }
    }


    // see if keep updating LB is worth it
    {
      if(!boundChanged)
      {
        // if(std::abs(sumLB - sumUB) < eps) return 2;
        if(absDiff(sumLB, sumUB) < relaEps) return 2;
        break;
      }
    }

  }


  return 1;
}


}










