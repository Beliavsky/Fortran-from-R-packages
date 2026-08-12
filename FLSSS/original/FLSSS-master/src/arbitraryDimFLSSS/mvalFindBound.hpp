#pragma once
#include "triMat.hpp"
#include "makeKsumIndicatorTables.hpp"


template<typename ing>
inline ing LBiFind(ing &ciLB, uint64_t ***M, ing ci_1LB, uint64_t *SR,
                   ing d, ing I, ing &J, ing *UB)
{

  ciLB = std::max<ing> (ci_1LB + 1, ciLB);
  uint64_t **v = M[0];
  mvalPlus<ing> (SR, SR, v[UB[I]], d);


  // std::cout << "In LBiFind, SR = ";
  // for(int k = 0, kend = d; k < kend; ++k) std::cout << SR[k] << ", ";
  // std::cout << std::endl;
  // std::cout << "Is SR negative = " <<
  //   bool(isNegative<ing> (SR, d)) << std::endl;


  // when I == J, certainly UB[I] >= ciLB and it breaks
  while ( UB[J] < ciLB - (I - J) )
  {
    mvalMinus<ing> (SR, SR, v[UB[J]], d);
    ++J;
  }


  // The 2nd force is that v[UB[J]]+v[UB[J]+1]+...+v[UB[J]+I-J] must >= SR
  while(true)
  {
    bool SRisNegative = isNegative<ing> (SR, d);
    // std::cout << "In LBiFind first while(), SRisNegative = " << SRisNegative << ", ";


    // What is the upper bound of I? UB[I], and the actual value is [UB[I]]
    if(J >= I)
    {
      if ( !SRisNegative and !greaterEqual<ing> (v[UB[I]], SR, d) ) return 0;
      break;
    }


    // Which column in M?: M[I-J]
    // What is the value we are looking for?:
    //   v[UB[J]] + v[UB[J] + 1] + ... + v[UB[J] + I - J]
    // Which row in that column?: M[I-J][UB[J]]
    if ( !SRisNegative and
         !greaterEqual<ing> (M[I - J][UB[J]], SR, d) )
    {
      mvalMinus<ing> (SR, SR, v[UB[J]], d);
      ++J;
    }
    else break;
  }


  // std::cout << "\nisNegative<ing> (SR, d) = " << isNegative<ing> (SR, d) << std::endl;
  if ( !isNegative<ing> (SR, d) )
  {
    ing I_J = I - J;
    // std::cout << "I_J = " << I_J << ", ";
    // std::cout << "ciLB - I_J = " << ciLB - I_J << ", ";
    // std::cout << "UB[J] = " << UB[J] << ", ";
    // std::cout << "&M[I_J][UB[J]] + 1 - &M[I_J][ciLB - I_J] = " <<
    //   &M[I_J][UB[J]] + 1 - &M[I_J][ciLB - I_J] << ", ";
    // std::cout << "&M[I_J][UB[J]] + 1 - &M[I_J][0] = " <<
    //   &M[I_J][UB[J]] + 1 - &M[I_J][0] << ", ";


    ciLB = mvalLowerBoundLr<ing> (
      &M[I_J][ciLB - I_J], &M[I_J][UB[J]] + 1, SR, d) - &M[I_J][0] + I_J;
    // std::cout << "ciLB = " << ciLB << std::endl;
  }


  return 1;
}




template<typename ing>
inline ing UBiFind(ing &ciUB, uint64_t ***M, ing ciP1UB, uint64_t *SR,
                   ing d, ing I, ing &J, ing *LB)
{
  ciUB = std::min<ing> (ciUB, ciP1UB - 1);
  uint64_t **v = M[0];
  mvalPlus(SR, SR, v[LB[I]], d);


  while( LB[J] > ciUB + (J - I) )
  {
    mvalMinus<ing> (SR, SR, v[LB[J]], d);
    --J;
  }


  // The 2nd force is that v[UB[J]]+v[UB[J]+1]+...+v[UB[J]+I-J] must >= SR
  while(true)
  {
    if(I == J)
    {
      if(!lessEqual<ing> ( v[LB[I]], SR, d) ) return 0;
      break;
    }


    if( !lessEqual<ing> ( M[J - I][LB[J] - (J - I)], SR, d) )
    {
      mvalMinus<ing> (SR, SR, v[LB[J]], d);
      --J;
    }
    else break;
  }


  ing J_I = J - I;
  ciUB = mvalUpperBoundLr<ing> (
      &M[J_I][LB[J] - J_I], &M[J_I][ciUB] + 1, SR, d) - &M[J_I][0];


  return 1;
}




template<typename ing>
inline ing findBoundCpp(
    ing len, ing d, uint64_t *targetS, ing *LB, uint64_t *sumLB,
    ing *UB, uint64_t *sumUB, uint64_t ***M, uint64_t *SR,
    KsumLookUpTable<ing> *look)
  // SRVcntr is a buffer of size at least d.
{


  if (!look->query(targetS, len, d)) return 0;


  if (!greaterEqual(sumUB, targetS, d) or !lessEqual(sumLB, targetS, d))
    return 0;


  bool boo = 0;
  uint64_t **v = M[0];


  unsigned LBsum = 0, UBsum = 0;
  while (true)
  {
    bool boundChanged = false;


    ing I = 0, J = 0;
    mvalPlusMinus(SR, targetS, v[UB[I]], sumUB, d);


    // Update first lower bound, 0;
    LBsum = 0;
    {
      ing tmpLB = LB[I];
      if ( !isNegative<ing> (SR, d) )
        LB[I] = mvalLowerBoundLr<ing> ( &v[LB[I]], &v[UB[I]] + 1, SR, d ) - v;
      if ( LB[I] > UB[I] ) return 0;
      boundChanged |= tmpLB != LB[I];
      std::memcpy(sumLB, v[LB[I]], sizeof(uint64_t) * d);
    }
    LBsum += LB[I];
    ++I;


    // Update lower bounds from 1 to the end
    {
      for(; I < len; ++I)
      {
        ing LBI = LB[I];
        ing b = LBiFind( LB[I], M, LB[I - 1], SR, d, I, J, UB );
        if ( b == 0 ) return 0;
        boundChanged |= LBI != LB[I];
        mvalPlus(sumLB, sumLB, v[LB[I]], d);
        LBsum += LB[I];
      }
    }


    // Can I jump out now ?
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


    boundChanged = false;
    J = len - 1; I = J;
    mvalPlusMinus(SR, targetS, v[LB[I]], sumLB, d);


    UBsum = 0;
    {
      ing tmpUB = UB[I];
      UB[I] = mvalUpperBoundLr<ing> (&v[LB[I]], &v[UB[I]] + 1, SR, d) - v;
      if(LB[I] > UB[I]) return 0;
      boundChanged |= tmpUB != UB[I];
      std::memcpy(sumUB, v[UB[I]], sizeof(uint64_t) * d);
      UBsum += UB[I];
    }


    --I;


    // Refresh the rest UB bounds
    {
      for(; I >= 0; --I)
      {
        ing UBI = UB[I];
        ing b = UBiFind<ing> (UB[I], M, UB[I + 1], SR, d, I, J, LB);
        if(b == 0) return 0;
        boundChanged |= UBI != UB[I];
        mvalPlus(sumUB, sumUB, v[UB[I]], d);
        UBsum += UB[I];
      }
    }


    // Can I jump out now ?
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




















