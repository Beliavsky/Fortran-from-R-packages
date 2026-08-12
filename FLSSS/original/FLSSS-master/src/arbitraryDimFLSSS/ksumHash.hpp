#pragma once
#include <Rcpp.h>
#include "arithmetic.hpp"
#include "tiktok.hpp"
#include "charlieThreadPool.hpp"


#define XXH_STATIC_LINKING_ONLY   /* access advanced declarations */
#define XXH_IMPLEMENTATION   /* access definitions */
#include "xxhash.h"


#ifndef vec
#define vec std::vector
#endif


struct ComputeComboRecur
{
  char B[8];
  int d, N, k;
  uint64_t prime;
  int *lb, *ub; // lb and ub should be of size k, the subset size.
  std::atomic<unsigned char> *Htable;
  uint64_t **v;
  vec<uint64_t> csumV;
  vec<uint64_t*> csumPtr;
  uint64_t **csum;
  // vec<int> initialLB;
  vec<int> lbvar;
  int *initialLB;


  ComputeComboRecur(){}


  void reset(int d, int N, int k, std::atomic<unsigned char> *Htable,
             uint64_t prime, uint64_t **v)
  {
    this->d = d;
    this->N = N;
    this->k = k;
    this->prime = prime;
    this->Htable = Htable;
    this->v = v;
    csumV.assign(uint64_t(d) * (k + 1), 0);
    csumPtr.resize(k + 1);
    for (int i = 0, iend = csumPtr.size(); i < iend; ++i)
      csumPtr[i] = &csumV[0] + i * uint64_t(d);
    csum = &csumPtr[1];


    unsigned char a = 1;
    for (int i = 0; i < 8; ++i) B[i] = a << i;
    // B[0] = 1;  B[1] = 2;  B[2] = 4;  B[3] = 8;
    // B[4] = 16; B[5] = 32; B[6] = 64; B[7] = 128;
  }


  void reset(int *lb, int *ub)
  {
    // this->lb = lb;
    // this->ub = ub;
    // initialLB.assign(lb, lb + k);
    this->initialLB = lb;
    this->ub = ub;
    lbvar.assign(lb, lb + k);
    this->lb = &lbvar[0];


    std::fill(csumV.begin(), csumV.end(), 0);
    for (int i = 0; i < k; ++i) mvalPlus<int> (
      csum[i], csum[i - 1], v[lb[i]], d);
  }


  ComputeComboRecur(int d, int N, int k, int *lb, int *ub,
                    std::atomic<unsigned char> *Htable,
                    // unsigned char *Htable,
                    uint64_t prime, uint64_t **v)
  {
    reset(d, N, k, Htable, prime, v);
    reset(lb, ub);
  }


  // Don't worry baby. You have done checking the unique hash values and they
  // ARE pretty unique in numeric experiments.
  void setBitFun(uint64_t *s)
  {
    uint64_t hashval = XXH64(s, sizeof(uint64_t) * d, 42); // Rcpp::Rcout << hashval << ", ";
    uint64_t whichBit = hashval % prime;
    uint64_t whichByte = whichBit >> 3; // whichBit / 8;
    uint64_t whichBitInByte = whichBit & 7;
    Htable[whichByte].fetch_or(B[whichBitInByte], std::memory_order_relaxed);


    // if (s[0] == 6300766154302661ull)
    // {
    //   Rcpp::Rcout << "\nFor 6300766154302661, whichBit = " << whichBit <<
    //     ", whichByte = " << whichByte << ", whichBitInByte = " <<
    //       whichBitInByte << "\n\n";
    //   Rcpp::Rcout << "(Htable[whichByte] & B[whichBitInByte]) = " <<
    //     (Htable[whichByte] & B[whichBitInByte]) << "\n\n";
    // }
  }


  // Compute sum and set bits.
  uint64_t operator()()
  {
    uint64_t n = 0;
    while (true) // k is the subset size. Will not change.
    {
      if (lb[k - 1] >= ub[k - 1])
      {
        mvalPlus<int> (csum[k - 1], csum[k - 2], v[lb[k - 1]], d);


        // if (csum[k - 1][0] == 6300766154302661ull) Rcpp::Rcout << "Hey that is it!\n";


        setBitFun(csum[k - 1]); // Register the sum
        ++n;
        int i = k - 1;
        for (; i >= 0 and lb[i] >= ub[i]; --i);
        if (i < 0) break;
        lb[i] += 1;
        mvalPlusMinus<int> (csum[i], csum[i], v[lb[i]], v[lb[i] - 1], d);
        i += 1;
        for (int iend = k - 1; i < iend; ++i)
        {
          lb[i] = std::max<int> (initialLB[i], lb[i - 1] + 1);
          mvalPlus<int> (csum[i], csum[i - 1], v[lb[i]], d);
        }
        lb[i] = std::max<int> (initialLB[i], lb[i - 1] + 1);
      }
      else
      {
        mvalPlus<int> (csum[k - 1], csum[k - 2], v[lb[k - 1]], d);


        // if (csum[k - 1][0] == 6300766154302661ull) Rcpp::Rcout << "Hey that is it!\n";


        setBitFun(csum[k - 1]); // Register the sum
        ++n;
        lb[k - 1] += 1;
      }
    }
    return n;
  }


};




// Print all the permutations of the first k elements.
void allCombo(int *initialLB, int *ub, int k, int bsize,
              vec<vec<int> > &rst)
{
  vec<int> lb(initialLB, initialLB + bsize);


  auto liftBounds = [](vec<int> &x)->void
  {
    for (int i = 1, iend = x.size(); i < iend; ++i)
      x[i] = std::max<int> (x[i - 1] + 1, x[i]);
  };


  while (true) // k is the subset size. Will not change.
  {
    if (lb[k - 1] >= ub[k - 1])
    {
      rst.push_back(lb);
      liftBounds(rst.back());
      int i = k - 1;
      for (; i >= 0 and lb[i] >= ub[i]; --i);
      if (i < 0) break;
      lb[i] += 1;
      for (++i; i < k; ++i)
        lb[i] = std::max<int> (initialLB[i], lb[i - 1] + 1);
    }
    else
    {
      rst.push_back(lb);
      liftBounds(rst.back());
      lb[k - 1] += 1;
    }
  }


  // Rcpp::Rcout << "All combo = \n";
  // for (int i = 0, iend = rst.size(); i < iend; ++i)
  // {
  //   for (int j = 0, jend = rst[i].size(); j < jend; ++j)
  //   {
  //     Rcpp::Rcout << rst[i][j] << ", ";
  //   }
  //   Rcpp::Rcout << "\n";
  // }


}


// Compute the maximum number of sums given bounds.
struct NofSums
{
  vec<vec<uint64_t> > x;
  uint64_t operator() (int *lb, int *ub, int len)
  {

    // Rcpp::Rcout << "LB = ";
    // for (int i = 0; i < len; ++i) Rcpp::Rcout << lb[i] << ", ";
    // Rcpp::Rcout << "\n";
    // Rcpp::Rcout << "UB = ";
    // for (int i = 0; i < len; ++i) Rcpp::Rcout << ub[i] << ", ";
    // Rcpp::Rcout << "\n";


    x.resize(len);
    x[0].resize(ub[0] - lb[0] + 1);
    std::iota(x[0].begin(), x[0].end(), 1);
    for (int i = 1; i < len; ++i)
    {
      // Rcpp::Rcout << "ub[i] - lb[i] + 1 = " << ub[i] - lb[i] + 1 << "\n";
      x[i].resize(ub[i] - lb[i] + 1);
      x[i].resize(0);
      int m = std::min<int> (lb[i] - lb[i - 1], x[i - 1].size());
      // Rcpp::Rcout << "m - 1 = " << m - 1 << ", ";
      // Rcpp::Rcout << "x[i - 1].size() = " << x[i - 1].size() << "\n";
      // Rcpp::Rcout << "x[i - 1][m - 1] = " << x[i - 1][m - 1] << "\n";
      x[i].push_back(x[i - 1][m - 1]);
      for (int j = lb[i] + 1, jend = ub[i] + 1; j < jend; ++j)
      {
        m = std::min<int> (j - lb[i - 1], x[i - 1].size());
        x[i].push_back(x[i - 1][m - 1] + x[i].back());
      }
    }


    // Rcpp::Rcout << "NofSums.x = \n";
    // for (int i = 0, iend = x.size(); i < iend; ++i)
    // {
    //   for (int j = 0, jend = x[i].size(); j < jend; ++j)
    //     Rcpp::Rcout << x[i][j] << ", ";
    //   Rcpp::Rcout << "\n";
    // }
    // Rcpp::Rcout << "\n";


    return x.back().back();
  }
};


struct Ksum
{
  int N, d, k;
  vec<uint64_t*> v;
  uint64_t Hsize, modPrime;
  NofSums NS;
  std::atomic<unsigned char> *H;
  int *lb, *ub;
  CharlieThreadPool *tp;


  void setHsizePrime(int upscale)
  {
    uint64_t primes[62] = {
      5ull, 11ull, 23ull, 47ull, 97ull, 199ull, 409ull, 823ull, 1741ull,
      3469ull, 6949ull, 14033ull, 28411ull, 57557ull, 116731ull, 236897ull,
      480881ull, 976369ull, 1982627ull, 4026031ull, 8175383ull, 16601593ull,
      33712729ull, 68460391ull, 139022417ull, 282312799ull, 573292817ull,
      1164186217ull, 2364114217ull, 4294967291ull, 8589934583ull, 17179869143ull,
      34359738337ull, 68719476731ull, 137438953447ull, 274877906899ull,
      549755813881ull, 1099511627689ull, 2199023255531ull, 4398046511093ull,
      8796093022151ull, 17592186044399ull, 35184372088777ull, 70368744177643ull,
      140737488355213ull, 281474976710597ull, 562949953421231ull,
      1125899906842597ull, 2251799813685119ull, 4503599627370449ull,
      9007199254740881ull, 18014398509481951ull, 36028797018963913ull,
      72057594037927931ull, 144115188075855859ull, 288230376151711717ull,
      576460752303423433ull, 1152921504606846883ull, 2305843009213693951ull,
      4611686018427387847ull, 9223372036854775783ull, 18446744073709551557ull};


    Hsize = NS(lb, ub, k) * upscale; // Container in NS has been filled.
    uint64_t which = std::lower_bound(
      &primes[0], &primes[0] + 62, Hsize) - &primes[0];
    modPrime = primes[which];
    Hsize = (modPrime + 7) / 8; // Hsize is the number of bytes in H.
  }


  // lb and ub will not be changed.
  void reset(uint64_t *x, int N, int d, int k, int *lb, int *ub,
             int upscale, CharlieThreadPool &tp)
  {
    this->tp = &tp;
    v.resize(N);
    for (int i = 0; i < N; ++i) v[i] = x + i * uint64_t(d);
    this->N = N;
    this->d = d;
    this->k = k;
    this->lb = lb;
    this->ub = ub;
    setHsizePrime(upscale); // Set Hsize.
    // Rcpp::Rcout << "modPrime = " << modPrime << "\n";
    // Rcpp::Rcout << "Hsize = " << Hsize << "\n";
  }


  Rcpp::RawVector operator() (bool verbose)
  {
    H = new std::atomic<unsigned char> [Hsize];
    auto empfun = [](std::size_t t){ return false; };
    tp->parFor(0, Hsize, [&](std::size_t i, std::size_t t)
    {
      H[i].store(0, std::memory_order_relaxed);
      return false;
    }, 1000, empfun, empfun);


    uint64_t Njobs = 50 * tp->maxCore;
    int i = 0, iend = NS.x.size();
    for (; i < iend and NS.x[i].back() < Njobs; ++i);
    if (verbose) Rcpp::Rcout << "Bounded " << k << "-sum computes & hashes ";
    if (i >= iend) // The number of sums is so small that it can be easily done
      // using just 1 core.
    {
      // std::atomic<uint64_t> ic(0);
      if (verbose) Rcpp::Rcout << NS.x.back().back() << " sums with 1 thread: ";
      uint64_t Nsums = ComputeComboRecur(d, N, k, lb, ub, H, modPrime, &v[0])();
      if (verbose) Rcpp::Rcout << Nsums << "\n\n";
    }
    else
    {
      Njobs = NS.x[i].back();
      vec<vec<int> > allLowerBounds;
      allCombo(lb, ub, i + 1, k, allLowerBounds);


      vec<ComputeComboRecur> computeComboObjs(tp->maxCore);
      for ( unsigned j = 0; j < tp->maxCore; ++j )
        computeComboObjs[j].reset(d, N, k, H, modPrime, &v[0]);
      vec<vec<int> > upperBounds(tp->maxCore, vec<int> (ub, ub + k));


      if (verbose) Rcpp::Rcout << NS.x.back().back() << " sums with " <<
        tp->maxCore << " threads: ";


      vec<uint64_t> icv(tp->maxCore, 0);
      tp->parFor(0, allLowerBounds.size(), [&](std::size_t j, std::size_t t)
      {
        if (verbose and t == 0)
          Rcpp::Rcout << std::accumulate(icv.begin(), icv.end(), 0) << ", ";


        std::copy(&allLowerBounds[j][0], &allLowerBounds[j][i + 1],
                  &upperBounds[t][0]);


        computeComboObjs[t].reset(&allLowerBounds[j][0], &upperBounds[t][0]);
        icv[t] += computeComboObjs[t]();
        return false;
      }, 1, empfun, empfun);


      if (verbose) Rcpp::Rcout <<
        std::accumulate(icv.begin(), icv.end(), 0) << ".\n";
    }


    Rcpp::RawVector rst(Hsize);
    unsigned char *rstPtr = &rst[0];
    for (uint64_t i = 0; i < Hsize; ++i)
      rstPtr[i] = H[i].load(std::memory_order_relaxed);


    // Rcpp::Rcout << "int(rstPtr[54844] & 32) = " << int(rstPtr[54844] & 32) << "\n";


    if (verbose)
    {
      // unsigned char eight[8] = {1, 2, 4, 8, 16, 32, 64, 128};
      unsigned char eight[8], a = 1;
      for (int i = 0; i < 8; ++i) eight[i] = a << i;
      vec<uint64_t> onesv(tp->maxCore, 0);
      tp->parFor(0, Hsize, [&](std::size_t i, std::size_t t)
      {
        for (int k = 0; k < 8; ++k)
          onesv[t] += (rstPtr[i] & eight[k]) != 0;
        return false;
      }, 1000, empfun, empfun);
      uint64_t ones = std::accumulate(onesv.begin(), onesv.end(), 0);
      Rcpp::Rcout << "Indicator table takes " <<
        rst.size() / float(1e9) << " GB. ";
      Rcpp::Rcout << (ones + 0.0) / (uint64_t(8) * rst.size()) * 100.0 <<
        "% of all bits are 1. Modulo prime = " <<
          std::to_string(modPrime) << ".\n\n";
    }


    delete [] H;
    return rst;
  }
};













