#pragma once
#include "mvalFindBound.hpp"
#include "ksumHash.hpp"
// #include "mflsssOBJ.hpp"


struct PtrPrime
{
  unsigned char *p;
  uint64_t prime;
  PtrPrime() { p = nullptr; prime = 0; }
  void reset(unsigned char *p, uint64_t prime)
  {
    this->p = p; this->prime = prime;
  }
};


template<typename ing>
struct KsumLookUpTable
{
  unsigned char B[8];
  vec<PtrPrime> Q;


  // Always construct an object with proper size.
  // Do not add a constructor with no arguments.
  KsumLookUpTable(ing fullSubsetSize)
  {
    unsigned char a = 1;
    for (int i = 0; i < 8; ++i) B[i] = a << i;
    // B[0] = a << ;  B[1] = 2;  B[2] = 4;  B[3] = 8;
    // B[4] = 16; B[5] = 32; B[6] = 64; B[7] = 128;
    Q.assign(fullSubsetSize + 1, PtrPrime());
  }


  // Results will be stored in lb and ub.
  // void gatherBounds(vec<mflsssOBJ<ing> > &mflsssObjList,
  bool gatherBounds(vec<ing*> &existingLBs,
                    vec<ing*> &existingUBs,
                    vec<ing>  &existingBsizes,
                    ing supersetSize,
                    ing *lb, ing *ub, ing len)
  {
    std::fill(lb, lb + len, supersetSize);
    std::fill(ub, ub + len, 0);
    for (int i = 0, iend = existingBsizes.size(); i < iend; ++i)
    {
      ing l = existingBsizes[i];
      if (l >= len)
      {
        for (ing j = 0; j < len; ++j)
        {
          lb[j] = std::min<ing> (lb[j], existingLBs[i][j]);
          ub[j] = std::max<ing> (ub[j], existingUBs[i][l - len + j]);
        }
      }
    }
    // It is possible that all mflsss objects have been exhausted.
    return std::accumulate(ub, ub + len, 0ll) != 0;
  }


  // Will not consider 2-sum, 1-sum.
  // Make the lookup tables from information in mflsssObjList.
  Rcpp::List make( // vec<mflsssobj> &mflsssObjList,
      vec<ing*> &existingLBs,
      vec<ing*> &existingUBs,
      vec<ing>  &existingBsizes,
      ing targetSubsetSize, ing supersetSize, ing d,
      uint64_t *v, ing maxK, int upscale, CharlieThreadPool &tp, bool verbose)
  {
    maxK = std::min<ing> (targetSubsetSize, std::max<ing> (3, maxK));
    vec<ing> bv(maxK * 2);
    ing *lb = &bv[0], *ub = lb + maxK;
    Ksum KS;
    Rcpp::List rst(std::max<int> (0, maxK - 3 + 1));
    Rcpp::StringVector enames(rst.size());
    Q.assign(targetSubsetSize + 1, PtrPrime()); // f.subsetSize >= maxK;
    vec<int> lbubint(maxK * 2);
    Rcpp::StringVector eleName(2);
    eleName[0] = "prime"; eleName[1] = "table";


    for (ing k = 3; k <= maxK; ++k)
    {
      bool tableCanBeBuilt = gatherBounds(
        existingLBs, existingUBs, existingBsizes, supersetSize, lb, ub, k);

      if (!tableCanBeBuilt)
      {
        if (verbose) Rcpp::Rcout <<
          "All solutions have been found during decomposition. \
No k-sum lookup tables to build.\n\n";
        return Rcpp::List::create();
      }
      std::copy(lb, lb + k, &lbubint[0]);
      std::copy(ub, ub + k, &lbubint[k]);


      if (verbose)
      {
        Rcpp::Rcout << "Gather index bounds from flsss objects..\n";
        Rcpp::Rcout << "Index lower bounds for " << int(k) << "-sum: ";
        for (int u = 0, uend = k; u < uend; ++u) Rcpp::Rcout << lbubint[u] << ", ";
        Rcpp::Rcout << "\n";
        Rcpp::Rcout << "Index Upper bounds for " << int(k) << "-sum: ";
        for (int u = 0, uend = k; u < uend; ++u) Rcpp::Rcout << lbubint[k + u] << ", ";
        Rcpp::Rcout << "\n";
      }


      KS.reset(v, supersetSize, d, k, &lbubint[0],
               &lbubint[k], upscale, tp);
      Rcpp::RawVector tmp = KS(verbose);


      Rcpp::RawVector tmpprime = copy2rRaw<uint64_t> (KS.modPrime);
      auto elelist = Rcpp::List::create(tmpprime, tmp);
      elelist.names() = eleName;
      rst[k - 3] = elelist;
      Q[k].reset(&tmp[0], KS.modPrime);
      enames[k - 3] = std::to_string(k) + "-sum";
    }
    rst.names() = enames;


    return rst;
  }


  void read(Rcpp::List tables, ing fullSubsetSize)
  {
    ing qsize = std::max<ing> (tables.size() + 3, fullSubsetSize);
    Q.assign(qsize, PtrPrime());
    for (int k = 0, kend = tables.size(); k < kend; ++k)
    {
      Rcpp::List tmp = tables[k];
      Rcpp::RawVector prime = tmp["prime"];
      Rcpp::RawVector table = tmp["table"];
      Q[k + 3].reset(&table[0], copyRraw<uint64_t> (prime));
      // Rcpp::Rcout << Q[k + 3].prime << "\n";
      // Rcpp::Rcout << uint64_t(Q[k + 3].p) << "\n\n";
    }
    // Rcpp::Rcout << "table this = " << uint64_t(this) << "\n";
  }


  // Return true if sum could be reached.
  unsigned char query(uint64_t *sum, ing len, ing d)
  {

    // if (sum[0] == 6300766154302661ull)
    // {
    //   Rcpp::Rcout << "What is going on?\n\n";
    //   Rcpp::Rcout << "len = " << int(len) << "\n";
    //   Rcpp::Rcout << "(Q[len][54844] & 32) = " << int(Q[len].p[54844] & 32) << "\n";
    // }



    // Rcpp::Rcout << "Q.size() = " << Q.size() << ", ";
    // Rcpp::Rcout << "len = " << int(len) << "\n";
    auto &q = Q[len];
    if (q.p == nullptr) return 1;
    uint64_t whichBit = XXH64(sum, sizeof(uint64_t) * d, 42) % q.prime;


    // if (sum[0] == 6300766154302661ull)
    // {
    //   Rcpp::Rcout << "\nDuring query, for 6300766154302661, whichBit = " << whichBit <<
    //     ", whichByte = " << (whichBit >> 3) << ", whichBitInByte = " <<
    //       (whichBit & 7) << "\n\n";
    //   Rcpp::Rcout << "q.p[whichBit >> 3] & B[whichBit & 7] = " <<
    //     int(q.p[whichBit >> 3] & B[whichBit & 7]) << "\n";
    //   Rcpp::Rcout << "B[whichBit & 7] = " << int(B[whichBit & 7]) << "\n";
    //   Rcpp::Rcout << "B = ";
    //   for (int u = 0; u < 8; ++u) Rcpp::Rcout << int(B[u]) << ", ";
    //   Rcpp::Rcout << "\n";
    // }


    return q.p[whichBit >> 3] & B[whichBit & 7];
  }
};














