# pragma once
# include <fstream>
# include <RcppParallel.h>
# include <mutex>
# include <atomic>
# include <chrono>
# include <string>
# include "macros.hpp"
# include "mPATclass.hpp"
# include "dnyTasking.hpp"


// ==========================================================================={
// Move std::vector from source to destination, an invalid vector
// (header points to somewhere d-allocating which will raise undefined
// behavior) by moving the bits in the headers. This differs from
// std::swap().
// template<typename T>
// void vhmove(vec<T> &des, vec<T> &source)
// {
//   std::memcpy((char*)(&des), &source, sizeof(vec<T>));
//   std::fill((char*)(&source), (char*)(&source) + sizeof(vec<T>), 0);
// }
template<typename T>
Rcpp::RawVector copy2rRaw(T &x)
{
  Rcpp::RawVector rst(sizeof(T));
  std::memcpy(&rst[0], &x, sizeof(T));
  return rst;
}
template<typename T>
Rcpp::RawVector copyVec2rRaw(vec<T> &x)
{
  Rcpp::RawVector rst(x.size() * sizeof(T));
  std::memcpy(&rst[0], &x[0], rst.size());
  return rst;
}
template<typename T>
void copyRraw(T &x, Rcpp::RawVector v)
{
  std::memcpy((char*)(&x), (char*)(&v[0]), sizeof(T));
  // while(true){};
}
template<typename T>
void copyRraw2vec(vec<T> &x, Rcpp::RawVector v)
{
  x.resize(v.size() / sizeof(T));
  std::copy(v.begin(), v.end(), (char*)(&x[0]));
}
template<typename T>
void zeroVecH(vec<T> &x)
{
  std::fill((char*)(&x[0]), (char*)(&x[0]) + sizeof(x), 0);
}
// ===========================================================================}




template<typename valtype, typename indtype>
struct shared
{
  // bool useBiSrch;
  indtype subsetSize;
  indtype N, d, dlst, dl, dust, du;
  int sizeNeed;
  // tbb::atomic<int> totalSize;
  std::atomic<int> totalSize;
  std::chrono::time_point<std::chrono::steady_clock> endTime;
  valtype ***M;
  INT *mask; // could be nullptr


  shared() { totalSize = 0; }


  shared(indtype subsetSize,
         indtype N, indtype d, indtype dlst, indtype dl,
         indtype dust, indtype du, int sizeNeed,
         std::chrono::time_point<std::chrono::steady_clock> endTime,
         valtype ***M, INT *mask):
    subsetSize(subsetSize),
    N(N), d(d), dlst(dlst), dl(dl),
    dust(dust), du(du), sizeNeed(sizeNeed),
    endTime(endTime), M(M), mask(mask)
  {
    totalSize = 0;
  }


  // the following 4 members work for knapsack()
  // valtype *profitVec;
  double *profitVec;
  indtype *optimalSolution;
  // valtype optimalProfit;
  double optimalProfit;


  shared(indtype subsetSize,
         indtype N, indtype d, indtype dlst, indtype dl,
         indtype dust, indtype du,
         std::chrono::time_point<std::chrono::steady_clock> endTime,
         valtype ***M, INT *mask,
         double *profitVec,
         indtype *optimalSolution):
    subsetSize(subsetSize),
    N(N), d(d), dlst(dlst), dl(dl),
    dust(dust), du(du), endTime(endTime), M(M), mask(mask),
    profitVec(profitVec), optimalSolution(optimalSolution)
  {
    optimalProfit = 0;
    std::fill(optimalSolution, optimalSolution + subsetSize, 0);
    sizeNeed = 0;
    totalSize = 0;
  }


  Rcpp::RawVector save() { return copy2rRaw(*this); }


  void read(Rcpp::RawVector x)
  {
    std::memcpy((char*)(this), &x[0], x.size());
  }


  void read(valtype ***M, INT *mask,
            double *profitVec, indtype *optimalSolution)
  {
    this->M = M; // Auxiliary triangle matrix.
    this->mask = mask;
    this->profitVec = profitVec;
    this->optimalSolution = optimalSolution;
  }
};


// Return true if time is up.
struct timer
{
  std::size_t i, delta;
  std::chrono::time_point<std::chrono::steady_clock> endtime;
  timer(std::chrono::time_point<std::chrono::steady_clock> endtime,
        std::size_t d = 64) { i = 0; delta = d; this->endtime = endtime; }
  bool operator() ()
  {
    ++i;
    if(i % delta != 0) return false;
    return std::chrono::steady_clock::now() > endtime;
  }
};


template<typename valtype, typename indtype, bool mk, bool useBiSearch>
struct dummyContainers
{
  vec<indtype> hopeV;
  vec<mPAT<valtype, indtype, mk, useBiSearch> > SKvec;
  vec<indtype> indvec;
  vec<valtype> valvec;
  vec<valtype> SRVcntr;
  void swap(vec<indtype> &ahopeV,
            vec<mPAT<valtype, indtype, mk, useBiSearch> > &aSKvec,
            vec<indtype> &aindvec,
            vec<valtype> &avalvec, vec<valtype> &aSRVcntr)
  {
    hopeV.swap(ahopeV);
    SKvec.swap(aSKvec);
    indvec.swap(aindvec);
    valvec.swap(avalvec);
    SRVcntr.swap(aSRVcntr);
  }
};




template<typename valtype, typename indtype, bool mk, bool useBiSearch>
struct mflsssOBJ
{
  indtype *hope; // '*hope' points to the first element to write in 'hopeV'
  shared<valtype, indtype> *f;
  vec<indtype> hopeV;


  double existingProfitSum;
  // Memorize the subset sum in the extra dimension (profit dimension).
  // It is only used in knapsack problem.


  vec<mPAT<valtype, indtype, mk, useBiSearch> > SKvec;
  mPAT<valtype, indtype, mk, useBiSearch> *SKback;
  vec<indtype> indvec;
  vec<valtype> valvec;
  vec<valtype> SRVcntr;
  vec<vec<indtype> > result;


  void print()
  {
    Rcpp::Rcout << "hopeV = ";
    for(int i = 0, iend = hope - &hopeV[0]; i < iend; ++i)
      Rcpp::Rcout << (int)hopeV[i] << ", ";
    Rcpp::Rcout << "\n";
    Rcpp::Rcout << "-------------------------\n";


    for(int i = 0, iend = SKback - &SKvec[0]; i < iend; ++i)
    {
      mPAT<valtype, indtype, mk, useBiSearch> &X = SKvec[i];


      Rcpp::Rcout << "beenUpdated = " << (int)X.beenUpdated << "\n";
      Rcpp::Rcout << "position = " << (int)X.position << "\n";
      Rcpp::Rcout << "len = " << (int)X.len << "\n";
      Rcpp::Rcout << "Nzeroed = " << (int)X.Nzeroed << "\n";
      Rcpp::Rcout << "beenUpdated = " << (int)X.beenUpdated << "\n";


      Rcpp::Rcout << "LB = ";
      for(int k = 0, kend = X.len; k < kend; ++k)
        Rcpp::Rcout << (int)X.LB[k] << ", ";
      Rcpp::Rcout << "\n";


      Rcpp::Rcout << "UB = ";
      for(int k = 0, kend = X.len; k < kend; ++k)
        Rcpp::Rcout << (int)X.UB[k] << ", ";
      Rcpp::Rcout << "\n";


      Rcpp::Rcout << "Bresv = ";
      for(int i = 0, iend = X.position + 1; i < iend; ++i)
        Rcpp::Rcout << (int)X.Bresv[i] <<", ";
      Rcpp::Rcout << "\n";


      Rcpp::Rcout << "MIN = ";
      for(int i = 0, iend = f->dl; i < iend; ++i) Rcpp::Rcout << X.MIN[i] << ", ";
      Rcpp::Rcout << "\n";


      Rcpp::Rcout << "MAX = ";
      for(int i = 0, iend = f->du; i < iend; ++i) Rcpp::Rcout << X.MAX[i] << ", ";
      Rcpp::Rcout << "\n";


      Rcpp::Rcout << "sumLB = ";
      for(int i = 0, iend = f->d; i < iend; ++i) Rcpp::Rcout << X.sumLB[i] << ", ";
      Rcpp::Rcout << "\n";


      Rcpp::Rcout << "sumUB = ";
      for(int i = 0, iend = f->d; i < iend; ++i) Rcpp::Rcout << X.sumUB[i] << ", ";
      Rcpp::Rcout << "\n";


      Rcpp::Rcout << "sumBresv = ";
      for(int i = 0, iend = f->d; i < iend; ++i) Rcpp::Rcout << X.sumBresv[i] << ", ";
      Rcpp::Rcout << "\n";
      Rcpp::Rcout << "-------------------------\n";
    }
    Rcpp::Rcout << "================================\n";

  }


  void swap(mflsssOBJ &X)
  {
    std::swap(X.hope, hope);
    std::swap(X.f, f);
    std::swap(X.hopeV, hopeV);
    std::swap(X.existingProfitSum, existingProfitSum);
    std::swap(X.SKvec, SKvec);
    std::swap(X.SKback, SKback);
    std::swap(X.indvec, indvec);
    std::swap(X.valvec, valvec);
    std::swap(X.SRVcntr, SRVcntr);
    std::swap(X.result, result);
  }


  void setAnSK(mPAT<valtype, indtype, mk, useBiSearch> *sk,
               indtype *ind, valtype *val, indtype len)
  {
    sk->beenUpdated = 1;
    sk->MIN = val;
    sk->MAX = sk->MIN + f->dl;
    sk->sumLB = sk->MAX + f->du;
    sk->sumUB = sk->sumLB + f->d;
    sk->sumBresv = sk->sumUB + f->d;
    sk->LB = ind;
    sk->UB = sk->LB + len;
    sk->Bresv = sk->UB + len;
    sk->len = len;
  }


  void initialize(shared<valtype, indtype> *fixedInfo,
                  valtype *target, valtype *ME,
                  indtype *LB, indtype *UB,
                  dummyContainers<valtype, indtype, mk, useBiSearch> *dummies = nullptr)
  {
    f = fixedInfo;
    std::size_t stackLen = (unsigned)f->subsetSize + 2;
    unsigned biscaleFactor = (unsigned)std::log2(f->N + 1.0 - f->subsetSize) + 3;
    // about '+3': once, there was a failed test case given '+0'
    if(dummies != nullptr)
    {
      dummies->swap(hopeV, SKvec, indvec, valvec, SRVcntr);
    }
    indvec.assign(stackLen * (stackLen + 1) / 2 * 3 * biscaleFactor, 0);
    valvec.assign((3 * (std::size_t)f->d + (std::size_t)f->dl +
      (std::size_t)f->du) * stackLen * biscaleFactor, 0);
    SKvec.resize((unsigned)f->subsetSize * biscaleFactor);
    hopeV.assign(f->subsetSize, 0);
    SRVcntr.assign(f->d, 0);
    hope = &hopeV[0];
    // extraDimSum = 0; // No need to initialize this for subset sum
    mPAT<valtype, indtype, mk, useBiSearch> *SKbegin = &SKvec.front();


    setAnSK(SKbegin, &indvec[0], &valvec[0], f->subsetSize);
    SKback = SKbegin + 1;


    for(indtype i = 0; i < SKbegin->len; ++i)
    {
      SKbegin->LB[i] = LB[i];
      SKbegin->UB[i] = UB[i];
    }


    // assign MIN and MAX
    {
      for(indtype i = f->dlst, iend = f->dlst + f->dl; i < iend; ++i)
      {
        SKbegin->MIN[i - f->dlst] = target[i] - ME[i];
      }
      for(indtype i = f->dust, iend = f->dust + f->du; i < iend; ++i)
      {
        SKbegin->MAX[i - f->dust] = target[i] + ME[i];
      }
    }


    iterSum<valtype, indtype> (SKbegin->sumLB, f->M[0], SKbegin->LB, SKbegin->len, f->d);
    iterSum<valtype, indtype> (SKbegin->sumUB, f->M[0], SKbegin->UB, SKbegin->len, f->d);


    result.resize(7);
    result.resize(0);
  }


  void finalize(dummyContainers<valtype, indtype, mk, useBiSearch> *dummies = nullptr)
  {
    if(dummies != nullptr)
      dummies->swap(hopeV, SKvec, indvec, valvec, SRVcntr);
  }


  int TTTstackRun( // std::ofstream *outfile = nullptr
  )
  {
    mPAT<valtype, indtype, mk, useBiSearch> *SK = &SKvec[0];
    int rstCurrentSize = result.size();


    // std::ofstream outfile("debug.csv", std::ofstream::app);
    // outfile << "output\n\n\n";
    // outfile << "(SKback - 1)->len = " << int((SKback - 1)->len) << "\n";


    // Rcpp::Rcout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n";
    // Rcpp::Rcout << "(SKback - 1)->len = " << int((SKback - 1)->len) << "\n";
    // Rcpp::Rcout << "hope - &hopeV[0] = " << int(hope - &hopeV[0]) << "\n";


    valtype **V = f->M[0];


    if((SKback - 1)->len == 1)
    {
      SK = SKback - 1;
      // for(indtype i = 0; i < f->N; ++i)
      for(indtype i = SK->LB[0], iend = SK->UB[0] + 1; i < iend; ++i)
      {
        bool allBetween = 1;
        for(indtype k = 0; k < f->dl and allBetween; ++k)
        {
          if(V[i][k + f->dlst] < SK->MIN[k]) allBetween = 0;
        }
        for(indtype k = 0; k < f->du and allBetween; ++k)
        {
          if(V[i][k + f->dust] > SK->MAX[k]) allBetween = 0;
        }
        if(allBetween)
        {
          *hope = i;
          result.push_back(hopeV);
        }
      }
      // update totalSize
      {
        int addSize = result.size() - rstCurrentSize;
        // if(addSize > 0) f->totalSize.fetch_add(addSize);
        if(addSize > 0) f->totalSize.fetch_add(addSize);
      }
      return SK - &SKvec[0];
    }


    timer atimer(f->endTime, 1);
    while(true)
    {
      // outfile << "rstCurrentSize = " << (int)rstCurrentSize << ", ";
      // (SKback - 1)->print(f->d, f->dl, f->du, outfile);
      // outfile << "parent printed ___________________________________\n\n";


      // Rcpp::Rcout << "Parent LB = \n";
      // for(int u = 0, uend = SKback[-1].len; u < uend; ++u)
      // {
      //   Rcpp::Rcout << (int)SKback[-1].LB[u] << ", ";
      // }
      // // Rcpp::Rcout << "\nParent UB = \n";
      // for(int u = 0, uend = SKback[-1].len; u < uend; ++u)
      // {
      //   Rcpp::Rcout << (int)SKback[-1].UB[u] << ", ";
      // }
      // Rcpp::Rcout << "\n";



      SKback->copyParentGene(*(SKback - 1), f->d, f->dl, f->du);


      // SKback->print(f->d, f->dl, f->du, outfile);
      // outfile << "parent copied ___________________________________\n\n";


      indtype boo = SKback->grow(f->M, f->d, f->dlst, f->dl, f->dust, f->du, hope,
                                 f->mask, SRVcntr //, &outfile
      );


      // outfile << "boo == " << (int)boo << "\n";
      // SKback->print(f->d, f->dl, f->du, outfile);
      // outfile << "child grown ___________________________________\n\n";


      // continue to give birth.
      if(boo == 1)
      {
        ++SKback;
        continue;
      }


      if(boo == 3) // if len in the child becomes 1
      {
        indtype i = SKback->LB[0], iend = SKback->UB[0] + 1;
        for(; i < iend; ++i)
        {
          hopeV.back() = i;
          result.push_back(hopeV);

          // Rcpp::Rcout << "boo = 3, ";
          // Rcpp::Rcout << "i = " << i << "\n";
          // for(int u = 0, uend = hopeV.size(); u < uend; ++u)
            // Rcpp::Rcout << (int)hopeV[u] << ", ";
          // Rcpp::Rcout << "\n\n";
        }
      }
      else if(boo == 2) // if lower bounds and upper bounds overlap
      {


        // Rcpp::Rcout << "boo = 2\n";
        // Rcpp::Rcout << "hope - &hopeV[0] = " << int(hope - &hopeV[0]) << "\n";


        // for(int u = 0, uend = hopeV.size(); u < uend; ++u)
        //   Rcpp::Rcout << (int)hopeV[u] << ", ";
        // Rcpp::Rcout << "\n";


        std::copy(SKback->UB, SKback->UB + SKback->len, hope);
        result.push_back(hopeV);


        // Rcpp::Rcout << "UB = ";
        // for(int u = 0, uend = SKback->len; u < uend; ++u)
        //   Rcpp::Rcout << (int)SKback->UB[u] << ", ";
        // Rcpp::Rcout << "\n\n";
      }


      while(true)
      {
        bool updateBool = (SKback - 1)->update(f->M, f->d, f->dlst, f->dl, f->dust, f->du);


        // (SKback - 1)->print(f->d, f->dl, f->du, outfile);
        // outfile << "parent updated ___________________________________\n\n";


        if(updateBool != 0) break;


        hope -= (SKback - 1)->Nzeroed;
        --SKback;


        if(SKback - SK <= 1)
        {
          // update totalSize
          int addSize = result.size() - rstCurrentSize;
          if(addSize > 0) f->totalSize.fetch_add(addSize);
          // Rcpp::Rcout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n";
          return 0; // all the combinations have been tried
        }
      }


      // update totalSize
      {
        int addSize = result.size() - rstCurrentSize;
        if(addSize > 0) f->totalSize.fetch_add(addSize);
        rstCurrentSize += addSize;
      }


      if(f->totalSize >= f->sizeNeed or atimer()) break;
    }


    // Rcpp::Rcout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n";
    return SKback - SK;
  }


  void initializeForKnapsack(
      shared<valtype, indtype>*, valtype*, valtype*, indtype*, indtype*);


  // int TTTstackRunForKnapsack(tbb::spin_mutex*, bool);
  // int TTTstackRunForKnapsack(tbb::spin_mutex*);
  int TTTstackRunForKnapsack(std::mutex*);


  Rcpp::List save()
  {
    // Rcpp::RawVector self = copy2rRaw(*this);
    Rcpp::RawVector hopeV_ = copyVec2rRaw(hopeV);
    int hopeOffset = hope - &hopeV[0];


    // mPAT<valtype, indtype, mk, useBiSearch> *mpat = &SKvec[0];
    // Rcpp::RawVector SKvec_0ptr = copy2rRaw(mpat);
    int SKbackOffset = SKback - &SKvec[0];


    Rcpp::RawVector SKvec_ = copyVec2rRaw(SKvec);


    indtype *indvec_0ptr_ = &indvec[0];
    Rcpp::RawVector indvec_0ptr = copy2rRaw(indvec_0ptr_);
    Rcpp::RawVector indvec_ = copyVec2rRaw(indvec);


    valtype *valvec_0ptr_ = &valvec[0];
    Rcpp::RawVector valvec_0ptr = copy2rRaw(valvec_0ptr_);
    Rcpp::RawVector valvec_ = copyVec2rRaw(valvec);


    Rcpp::RawVector SRVcntr_ = copyVec2rRaw(SRVcntr);


    Rcpp::List result_(result.size());
    for(int i = 0, iend = result.size(); i < iend; ++i)
      result_[i] = copyVec2rRaw(result[i]);


    return Rcpp::List::create(Rcpp::Named("existingProfitSum") = existingProfitSum,
                              Rcpp::Named("hopeOffset") = hopeOffset,
                              Rcpp::Named("hopeV") = hopeV_,
                              Rcpp::Named("SKbackOffset") = SKbackOffset,
                              Rcpp::Named("SKvec") = SKvec_,
                              Rcpp::Named("indvec_0ptr") = indvec_0ptr,
                              Rcpp::Named("indvec") = indvec_,
                              Rcpp::Named("valvec_0ptr") = valvec_0ptr,
                              Rcpp::Named("valvec") = valvec_,
                              Rcpp::Named("SRVcntr") = SRVcntr_,
                              Rcpp::Named("result") = result_);
  }


  void read(Rcpp::List X, shared<valtype, indtype> *f)
  {
    this->f = f;
    // Rcpp::Rcout << "1.15\n";

    // Rcpp::RawVector x0 = X["self"];
    existingProfitSum = X["existingProfitSum"];
    // Rcpp::Rcout << "self size = " << sizeof(*this) << "\n";
    // Rcpp::Rcout << "x0 size = " << x0.size() << "\n";
    // return;
    // copyRraw(*this, x0);
    // std::memcpy((char*)(this), (char*)(&x0[0]), x0.size());
    // Rcpp::Rcout << "1.16\n";
    // zeroVecH(hopeV);
    // Rcpp::Rcout << "1.17\n";
    // zeroVecH(SKvec);
    // Rcpp::Rcout << "1.18\n";
    // zeroVecH(indvec);
    // Rcpp::Rcout << "1.19\n";
    // zeroVecH(valvec);
    // Rcpp::Rcout << "1.191\n";
    // zeroVecH(SRVcntr);
    // Rcpp::Rcout << "1.192\n";
    // zeroVecH(result);
    // Rcpp::Rcout << "1.193\n";
    // return;


    // Rcpp::Rcout << "1.2\n";

    Rcpp::RawVector x1 = X["hopeV"]; copyRraw2vec(hopeV, x1);

    // Rcpp::Rcout << "1.3\n";

    Rcpp::RawVector x2 = X["SKvec"]; copyRraw2vec(SKvec, x2);

    // Rcpp::Rcout << "1.31\n";

    Rcpp::RawVector x3 = X["indvec"]; copyRraw2vec(indvec, x3);

    // Rcpp::Rcout << "1.32\n";

    Rcpp::RawVector x4 = X["valvec"]; copyRraw2vec(valvec, x4);

    // Rcpp::Rcout << "1.33\n";

    Rcpp::RawVector x5 = X["SRVcntr"]; copyRraw2vec(SRVcntr, x5);

    // Rcpp::Rcout << "1.34\n";

    Rcpp::List x6 = X["result"];
    result.resize(x6.size());
    for(int i = 0, iend = result.size(); i < iend; ++i)
    {
      Rcpp::RawVector v = x6[i];
      copyRraw2vec(result[i], v);
    }
    // Rcpp::Rcout << "1.35\n";


    indtype *indvec_0ptr;
    Rcpp::RawVector x7 = X["indvec_0ptr"]; copyRraw(indvec_0ptr, x7);
    valtype *valvec_0ptr;
    Rcpp::RawVector x8 = X["valvec_0ptr"]; copyRraw(valvec_0ptr, x8);
    // mPAT<valtype, indtype, mk, useBiSearch> *SKvec_0ptr;
    // Rcpp::RawVector x9 = X["SKvec_0ptr"]; copyRraw(SKvec_0ptr, x9);


    int hopeOffset = X["hopeOffset"];
    hope = &hopeV[0] + hopeOffset;
    // SKback = SKback - SKvec_0ptr + &SKvec[0] + int(X["SKbackOffset"]);
    int SKbackOffset = X["SKbackOffset"];
    SKback = &SKvec[0] + SKbackOffset;
    // for(int i = 0, iend = SKvec.size(); i < iend; ++i)
    for(int i = 0; i < SKbackOffset; ++i)
    {
      mPAT<valtype, indtype, mk, useBiSearch> &x = SKvec[i];


      // std::cout << std::size_t(x.LB) << ", " << std::size_t(indvec_0ptr) << "\n";
      // std::cout << std::size_t(x.LB) - std::size_t(indvec_0ptr) << "\n";
      // std::cout << "\n";


      x.LB = x.LB - indvec_0ptr + &indvec[0];
      x.UB = x.UB - indvec_0ptr + &indvec[0];
      // for(int k = 0, kend = x.len; k < kend; ++k)
      // {
      //   Rcpp::Rcout << (int)x.LB[k] << ", ";
      // }
      // Rcpp::Rcout << "\n";
      x.Bresv = x.Bresv - indvec_0ptr + &indvec[0];
      x.MIN = x.MIN - valvec_0ptr + &valvec[0];
      x.MAX = x.MAX - valvec_0ptr + &valvec[0];
      x.sumLB = x.sumLB - valvec_0ptr + &valvec[0];
      // for(int k = 0, kend = x.len; k < kend; ++k)
      // {
      //   Rcpp::Rcout << x.sumLB[k] << ", ";
      // }
      // Rcpp::Rcout << "\n\n";
      x.sumUB = x.sumUB - valvec_0ptr + &valvec[0];
      x.sumBresv = x.sumBresv - valvec_0ptr + &valvec[0];
    }
  }
};




template<typename valtype, typename indtype, bool mk, bool useBiSearch>
indtype growTwin(mflsssOBJ<valtype, indtype, mk, useBiSearch> &Xmflsss,
                 mflsssOBJ<valtype, indtype, mk, useBiSearch> &Ymflsss,
                 std::ofstream *outfile = nullptr) // outfile prints log
{
  shared<valtype, indtype> &f = *Xmflsss.f;
  mPAT<valtype, indtype, mk, useBiSearch> &X = *Xmflsss.SKback;
  X.copyParentGene(*(Xmflsss.SKback - 1), f.d, f.dl, f.du);


  indtype boo = findBoundCpp<valtype, indtype, mk, useBiSearch> (
    X.len, f.d, f.dlst, f.dl, f.dust, f.du, X.MIN, X.MAX,
    X.LB, X.sumLB, X.UB, X.sumUB, f.M, f.mask, Xmflsss.SRVcntr);


  if(outfile != nullptr)
  {
    X.print(f.d, f.dl, f.du, *outfile);
    *outfile << "\n\nBounds found ___________________________________, boo = " << (int)boo << "\n";
  }


  if(boo == 0) return 0;
  if(X.len == 1) return 3;
  if(boo == 2) return 2;


  // find the slot that has the least gap
  // Rcpp::Rcout << "1.131\n";
  X.position = 0;
  indtype nonzeroMin = -1;


  vec<indtype> acntr(X.len);
  indtype *overlapPosition = &*acntr.begin(), *olpend = overlapPosition;


  // Rcpp::Rcout << "1.132\n";
  for(indtype i = 0; i < X.len; ++i)
  {
    indtype tmp = X.UB[i] - X.LB[i];
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


  // Rcpp::Rcout << "1.133\n";
  // erase all positions where LB and UB meet.
  indtype &Nzeroed = X.Nzeroed;
  Nzeroed = olpend - overlapPosition;
  if(Nzeroed > 0)
  {
    vec<valtype> cntrS(f.d); valtype *S = &*cntrS.begin();
    *olpend = X.len;
    for(indtype i = 0; i < Nzeroed; ++i)
    {
      indtype &st = overlapPosition[i], &end = overlapPosition[i + 1];
      mvalPlus(S, S, f.M[0][X.UB[st]], f.d);
      std::copy(X.LB + st + 1, X.LB + end, X.LB + st - i);
      std::copy(X.UB + st + 1, X.UB + end, X.UB + st - i);
    }
    X.len -= Nzeroed;


    mvalMinus(X.MIN, X.MIN, S + f.dlst, f.dl); // target changes, so MIN and MAX change
    mvalMinus(X.MAX, X.MAX, S + f.dust, f.du);
    mvalMinus(X.sumLB, X.sumLB, S, f.d);
    mvalMinus(X.sumUB, X.sumUB, S, f.d);


    // after erasion, position may change. Adjust position
    {
      indtype tmp = 0;
      for(indtype *i = overlapPosition; i < olpend; ++i)
      {
        if(X.position > *i) ++tmp;
        else break;
      }
      X.position -= tmp;
    }
  }


  // Rcpp::Rcout << "1.134\n";
  // x, pass wisdom to your twin!
  {
    Ymflsss.f = &f;
    Ymflsss.hopeV.assign(f.subsetSize, 0);
    // Rcpp::Rcout << "1.1341\n";
    std::copy(&Xmflsss.hopeV[0], Xmflsss.hope, &Ymflsss.hopeV[0]);
    Ymflsss.hope = &Ymflsss.hopeV[0] + (Xmflsss.hope - &Xmflsss.hopeV[0]);
    Ymflsss.SKvec.resize(Xmflsss.SKvec.size() - 1);
    // Rcpp::Rcout << "1.1342\n";
    Ymflsss.SKback = &Ymflsss.SKvec[1];
    Ymflsss.indvec.assign(Xmflsss.indvec.size() - 2 * (int)Xmflsss.SKvec[0].len, 0);
    // Rcpp::Rcout << "Xmflsss.indvec.size() - 2 * (int)Xmflsss.SKvec[0].len = "
                // << Xmflsss.indvec.size() - 2 * (int)Xmflsss.SKvec[0].len << "\n";
    // Rcpp::Rcout << "1.1343\n";
    Ymflsss.valvec.assign(Xmflsss.valvec.size() - 2 * ((int)f.d + f.dl + f.du), 0);
    // Rcpp::Rcout << "Xmflsss.valvec.size() - 2 * ((int)f.d + f.dl + f.du) = "
                // << Xmflsss.valvec.size() - 2 * ((int)f.d + f.dl + f.du) << "\n";
    // Rcpp::Rcout << "1.1344\n";
    Ymflsss.SRVcntr.assign(f.d, 0);
    // Rcpp::Rcout << "1.1345\n";
    // Ymflsss.result.reserve(f.sizeNeed);
    Ymflsss.result.reserve(7);
  }


  // Rcpp::Rcout << "1.135\n";
  // initialize Ymflsss.SKvec, the first element
  Ymflsss.setAnSK(&Ymflsss.SKvec[0], &Ymflsss.indvec[0],
                  &Ymflsss.valvec[0], Xmflsss.SKback->len);


  // X takes the lower half, Y takes the upper half
  mPAT<valtype, indtype, mk, useBiSearch> &Y = *(Ymflsss.SKback - 1);


  X.beenUpdated = 1;
  Y.beenUpdated = 1;
  std::copy(X.MIN, X.MIN + f.dl + f.du, Y.MIN);


  std::copy(X.sumUB, X.sumUB + f.d, Y.sumUB);
  std::copy(X.UB, X.UB + X.len, Y.UB);
  indtype cap = (X.UB[X.position] + X.LB[X.position]) / 2;
  indtype capResv = cap;
  indtype i = X.position;
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


  ++Xmflsss.SKback;
  return 1;
}




template<typename valtype, typename indtype, bool mk, bool useBiSearch>
void mitosis(vec<mflsssOBJ<valtype, indtype, mk, useBiSearch> > &descendants,
             shared<valtype, indtype> &f,
             vec<vec<indtype> > &rstCollection,
             indtype *LB, indtype *UB, valtype *target, valtype *ME,
             int threads, int avgThreadLoad)
{
  int Ndescendants = 1;
  if(threads > 1)
  {
    Ndescendants = 1 << ((int)std::log2(threads * avgThreadLoad + 0.0) + 1);
  }
  // Rcpp::Rcout << " Ndescendants = " << Ndescendants << "\n";
  descendants.resize(Ndescendants);
  descendants[0].initialize(&f, target, ME, LB, UB);
  vec<unsigned char> acntr(Ndescendants, 0);
  unsigned char *dead = &*acntr.begin();
  int j = 1;
  // Rcpp::Rcout << "1.11\n";


  while(j < Ndescendants)
  {
    int iend = j;
    for(int i = 0; i < iend; ++i, ++j)
    {
      if(f.totalSize >= f.sizeNeed) return;


      if(dead[i])
      {
        dead[j] = 1;
        continue;
      }


      // Rcpp::Rcout << "i = " << i << ", j = " << j << ", ";
      // Rcpp::Rcout << "descendants.size() = " << descendants.size() << "\n";
      indtype boo = growTwin(descendants[i], descendants[j] //, &mfile
      );
      // Rcpp::Rcout << "1.13\n";


      // std::cout << "grown\n";
      mPAT<valtype, indtype, mk, useBiSearch> &tmp = *descendants[i].SKback;
      if(boo == 0)
      {
        dead[i] = 1;
        dead[j] = 1;
      }
      else if(boo == 3)
      {
        indtype k = tmp.LB[0], kend = tmp.UB[0] + 1;
        for(; k < kend; ++k)
        {
          // std::cout << "k = " + std::to_string(k) + "\n";
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


  // Rcpp::Rcout << "1.2\n";
  if(f.totalSize >= f.sizeNeed) return;
  // Rcpp::Rcout << "1.3\n";


  // cleansing descendants[]
  int validDescendents = Ndescendants - std::accumulate(dead, dead + Ndescendants, (int)0);
  if(validDescendents == Ndescendants) return;
  // Rcpp::Rcout << "1.4\n";


  vec<mflsssOBJ<valtype, indtype, mk, useBiSearch> >
    descendantsRemain(validDescendents);
  for(int i = 0, k = 0; i < Ndescendants; ++i)
  {
    if(dead[i]) continue;
    descendants[i].swap(descendantsRemain[k]);
    ++k;
  }
  descendants.swap(descendantsRemain);
  // Rcpp::Rcout << "1.5\n";
}




// Grow once.
template<typename valtype, typename indtype, bool mk, bool useBiSearch>
void mitosis(
    vec<mflsssOBJ<valtype, indtype, mk, useBiSearch> > &descendants,
    mflsssOBJ<valtype, indtype, mk, useBiSearch> &root,
    shared<valtype, indtype> &f,
    vec<vec<indtype> > &rstCollection)
{
  // Rcpp::Rcout << "\n========================================mitosis\n";
  int Ndescendants = 2;
  descendants.resize(Ndescendants);
  // descendants[0].initialize(&f, target, ME, LB, UB);
  descendants[0].swap(root);
  vec<unsigned char> acntr(Ndescendants, 0);
  unsigned char *dead = &*acntr.begin();
  int j = 1;


  // Rcpp::Rcout << "\nwhile begin\n";
  while(j < Ndescendants)
  {
    int iend = j;
    for(int i = 0; i < iend; ++i, ++j)
    {
      if(f.totalSize >= f.sizeNeed) return;


      if(dead[i])
      {
        dead[j] = 1;
        continue;
      }


      // Rcpp::Rcout << "\ngrow twin begin\n";
      indtype boo = growTwin(descendants[i], descendants[j] //, &mfile
      );
      // Rcpp::Rcout << "\ngrow twin end\n";

      // std::cout << "grown\n";
      mPAT<valtype, indtype, mk, useBiSearch> &tmp = *descendants[i].SKback;
      if(boo == 0)
      {
        // Rcpp::Rcout << "boo == 0\n";
        dead[i] = 1;
        dead[j] = 1;
      }
      else if(boo == 3)
      {
        // Rcpp::Rcout << "boo == 3\n";
        indtype k = tmp.LB[0], kend = tmp.UB[0] + 1;
        for(; k < kend; ++k)
        {
          descendants[i].hopeV.back() = k;
          // Rcpp::Rcout << "rstCollection.push_back begin\n";
          rstCollection.push_back(descendants[i].hopeV);
          // Rcpp::Rcout << "rstCollection.size() = " << rstCollection.size() << "\n";
          ++f.totalSize;
        }
        dead[i] = 1;
        dead[j] = 1;
      }
      else if(boo == 2)
      {
        // Rcpp::Rcout << "boo == 2\n";
        std::copy(tmp.UB, tmp.UB + tmp.len, descendants[i].hope);
        rstCollection.push_back(descendants[i].hopeV);
        ++f.totalSize;
        dead[i] = 1;
        dead[j] = 1;
      }
    }
  }
  // Rcpp::Rcout << "\nwhile end\n";


  if(f.totalSize >= f.sizeNeed) return;


  // cleansing descendants[]
  int validDescendents = Ndescendants - std::accumulate(dead, dead + Ndescendants, (int)0);
  if(validDescendents == Ndescendants) return;


  vec<mflsssOBJ<valtype, indtype, mk, useBiSearch> >
    descendantsRemain(validDescendents);
  for(int i = 0, k = 0; i < Ndescendants; ++i)
  {
    if(dead[i]) continue;
    descendants[i].swap(descendantsRemain[k]);
    ++k;
  }
  descendants.swap(descendantsRemain);


  // Rcpp::Rcout << "\n========================================mitosis\n";
}




template<typename valtype, typename indtype, bool mk, bool useBiSearch>
struct parMflsssOBJ: public RcppParallel::Worker // works for mflsssComoPar() or fully initialized mflsssOBJs.
{
  vec<mflsssOBJ<valtype, indtype, mk, useBiSearch> > &mflsssOBJvec;
  shared<valtype, indtype> *f;
  dynamicTasking *dT;
  void operator() (std::size_t st, std::size_t end)
  {
    for(;;)
    {
      std::size_t objI = 0;
      if(!dT->nextTaskID(objI)) break;
      mflsssOBJvec[objI].TTTstackRun();
      if(f->totalSize >= f->sizeNeed or
           std::chrono::steady_clock::now() > f->endTime) break;
    }
  }


  parMflsssOBJ(vec<mflsssOBJ<valtype, indtype, mk, useBiSearch> >
                 &mflsssOBJvec, std::size_t maxCore): mflsssOBJvec(mflsssOBJvec)
  {
    dynamicTasking dt(maxCore, mflsssOBJvec.size());
    dT = &dt;
    f = mflsssOBJvec.begin()->f;
    RcppParallel::parallelFor(0, dT->NofCore, *this);
  }
};




template<typename valtype, typename indtype, bool mk, bool useBiSearch>
struct parMflsssOBJbyCore: public RcppParallel::Worker // works for mflsssPar()
{
  vec<mflsssOBJ<valtype, indtype, mk, useBiSearch> > &mflsssOBJvec; // mflsssOBJvec is only of size maxCore
  shared<valtype, indtype> *f;
  valtype *target;
  dummyContainers<valtype, indtype, mk, useBiSearch> *dummyCs;
  valtype *ME;
  indtype *commonLB, *commonUB;
  vec<vec<vec<indtype> > > &rst; // rst is of size
  dynamicTasking *dT;
  void operator() (std::size_t st, std::size_t end)
  {
    for(;;)
    {
      std::size_t objI = 0;
      if(!dT->nextTaskID(objI)) break;
      mflsssOBJvec[st].initialize(
          f, target + objI * f->d, ME, commonLB, commonUB, dummyCs + st);
      mflsssOBJvec[st].TTTstackRun();
      mflsssOBJvec[st].finalize(dummyCs + st);


      // harvest
      {
        for(int i = 0, iend = mflsssOBJvec[st].result.size(); i < iend; ++i)
        {
          rst[st].resize(rst[st].size() + 1);
          rst[st].back().swap(mflsssOBJvec[st].result[i]);
        }
      }


      if(f->totalSize >= f->sizeNeed or
           std::chrono::steady_clock::now() > f->endTime) break;
    }
  }


  parMflsssOBJbyCore(vec<mflsssOBJ<valtype, indtype, mk, useBiSearch> > &mflsssOBJvec,
                     shared<valtype, indtype> *f, valtype *target,
                     dummyContainers<valtype, indtype, mk, useBiSearch> *dummyCs,
                     valtype *ME, indtype *commonLB, indtype *commonUB,
                     vec<vec<vec<indtype> > > &rst,
                     unsigned keyTargetSize, std::size_t maxCore):
    mflsssOBJvec(mflsssOBJvec), f(f), target(target),
    dummyCs(dummyCs), ME(ME), commonLB(commonLB), commonUB(commonUB), rst(rst)
  {
    mflsssOBJvec.resize(maxCore);
    dynamicTasking dt(maxCore, keyTargetSize); dT = &dt;
    RcppParallel::parallelFor(0, dT->NofCore, *this);
  }
};







// ================================================================================================
// Knapsack problems
// ================================================================================================
template<typename valtype, typename indtype, bool mk, bool useBiSearch>
inline void mflsssOBJ<valtype, indtype, mk, useBiSearch>::
  initializeForKnapsack(
    shared<valtype, indtype> *fixedInfo, valtype *target,
    valtype *ME, indtype *LB, indtype *UB)
  {
    f = fixedInfo;
    std::size_t stackLen = (unsigned)f->subsetSize + 2;
    unsigned biscaleFactor = (unsigned)std::log2(f->N + 0.0 - f->subsetSize) + 3;
    indvec.assign(stackLen * (stackLen + 1) / 2 * 3 * biscaleFactor, 0);
    valvec.assign((3 * (std::size_t)f->d + (std::size_t)f->dl +
      (std::size_t)f->du) * stackLen * biscaleFactor, 0);
    SKvec.resize((unsigned)f->subsetSize * biscaleFactor);
    SRVcntr.assign(f->d, 0);


    hopeV.assign(f->subsetSize, 0);
    hope = &hopeV[0];
    existingProfitSum = 0;


    mPAT<valtype, indtype, mk, useBiSearch> *SKbegin = &SKvec.front();


    setAnSK(SKbegin, &indvec[0], &valvec[0], f->subsetSize);
    SKback = SKbegin + 1;


    for(indtype i = 0; i < SKbegin->len; ++i)
    {
      SKbegin->LB[i] = LB[i];
      SKbegin->UB[i] = UB[i];
    }


    // assign MIN and MAX
    {
      for(indtype i = f->dlst, iend = f->dlst + f->dl; i < iend; ++i)
      {
        SKbegin->MIN[i - f->dlst] = target[i] - ME[i];
      }
      for(indtype i = f->dust, iend = f->dust + f->du; i < iend; ++i)
      {
        SKbegin->MAX[i - f->dust] = target[i] + ME[i];
      }
    }


    iterSum<valtype, indtype> (SKbegin->sumLB, f->M[0], SKbegin->LB, SKbegin->len, f->d);
    iterSum<valtype, indtype> (SKbegin->sumUB, f->M[0], SKbegin->UB, SKbegin->len, f->d);
  }




template<typename valtype, typename indtype, bool mk, bool useBiSearch>
inline int mflsssOBJ<valtype, indtype, mk, useBiSearch>::
  TTTstackRunForKnapsack(
    // tbb::spin_mutex *mx // , bool verbose //, std::ofstream *outfile = nullptr
    std::mutex *mx
  )
  {
    mPAT<valtype, indtype, mk, useBiSearch> *SK = &SKvec[0];


    // std::ofstream outfile("proboutput.csv", std::ofstream::out|std::ofstream::app);
    // outfile << "output\n\n\n";


    // valtype **V = f->M[0];
    while(true)
    {
      // outfile << "rstCurrentSize = " << rstCurrentSize << ", ";
      // (SKback - 1)->print(f->d, f->dl, f->du, outfile);
      // outfile << "parent printed ___________________________________\n\n";


      SKback->copyParentGene(*(SKback - 1), f->d, f->dl, f->du);


      // SKback->print(f->d, f->dl, f->du, outfile);
      // outfile << "parent copied ___________________________________\n\n";


      indtype boo = SKback->growForKnapsack(
        f->M, f->d, f->dlst, f->dl, f->dust, f->du, hope,
        f->mask, f->profitVec, existingProfitSum, f->optimalProfit, SRVcntr);


      // outfile << "boo == " << (int)boo << "\n";
      // SKback->print(f->d, f->dl, f->du, outfile);
      // outfile << "child grown ___________________________________\n\n";


      // continue to give birth.
      if(boo == 1)
      {
        ++SKback;
        continue;
      }


      // if(boo == 3 or boo == 2)
      if(boo != 0)
      {
        std::copy(SKback->UB, SKback->UB + SKback->len, hope);
        valtype tmpProfit = 0;
        for(indtype i = 0; i < f->subsetSize; ++i)
        {
          tmpProfit += f->profitVec[hopeV[i]];
        }
        mx->lock();
        {
          if(tmpProfit > f->optimalProfit)
          {
            f->optimalProfit = tmpProfit;
            std::copy(hopeV.begin(), hopeV.end(), f->optimalSolution);
            // if(verbose) std::cout << "Updated profit = " << tmpProfit << "\n";
            // In Linux environment, the above still smashes C stack.
          }
        }
        mx->unlock();
      }


      while(true)
      {
        bool updateBool = (SKback - 1)->update(f->M, f->d, f->dlst, f->dl, f->dust, f->du);


        // (SKback - 1)->print(f->d, f->dl, f->du, outfile);
        // outfile << "parent updated ___________________________________\n\n";


        if(updateBool != 0) break;


        // Recover hope and the profit sum. Still struggling to decide if using an array to store previous
        // 'existingProfitSum' would be more efficient.
        {
          indtype *i = hope - 1;
          --SKback;
          hope -= SKback->Nzeroed;
          for(; i >= hope; --i)
          {
            existingProfitSum -= f->profitVec[*i];
          }
        }


        if(SKback - SK <= 1)
        {
          // std::cout << " All combinations have been tried\n";
          return 0; // all the combinations have been tried
        }
      }


      if(std::chrono::steady_clock::now() > f->endTime)
      {
        return -1;
      }
    }


    return SKback - SK;
  }




// This growTwin() function is made for knapsack
template<typename valtype, typename indtype, bool mk, bool useBiSearch>
inline indtype growTwinForKnapsack(
    mflsssOBJ<valtype, indtype, mk, useBiSearch> &Xmflsss,
    mflsssOBJ<valtype, indtype, mk, useBiSearch> &Ymflsss,
    std::ofstream *outfile = nullptr) // outfile prints log
{
  shared<valtype, indtype> &f = *Xmflsss.f;
  mPAT<valtype, indtype, mk, useBiSearch> &X = *Xmflsss.SKback;
  X.copyParentGene(*(Xmflsss.SKback - 1), f.d, f.dl, f.du);


  // X.print(f.d, f.dl, f.du, *outfile);
  // *outfile << "\n\n";


  indtype boo = findBoundCpp<valtype, indtype, mk, useBiSearch> (
    X.len, f.d, f.dlst, f.dl, f.dust, f.du, X.MIN, X.MAX,
    X.LB, X.sumLB, X.UB, X.sumUB, f.M, f.mask, Xmflsss.SRVcntr);


  // See if the sum of upper bounds is less than the current maximal profit
  if(boo != 0)
  {
    double S = Xmflsss.existingProfitSum;
    for(indtype i = 0; i < X.len; ++i)
    {
      S += f.profitVec[X.UB[i]];
    }
    if(S <= f.optimalProfit) return 0;
  }


  if(outfile != nullptr)
  {
    X.print(f.d, f.dl, f.du, *outfile);
    *outfile << "\n\nBounds found ___________________________________, boo = "
             << (int)boo << "\n";
  }


  if(boo == 0) return 0;
  if(X.len == 1) return 3;
  if(boo == 2) return 2;


  // find the slot that has the least gap
  X.position = 0;
  indtype nonzeroMin = -1;


  vec<indtype> acntr(X.len);
  indtype *overlapPosition = &*acntr.begin(), *olpend = overlapPosition;


  for(indtype i = 0; i < X.len; ++i)
  {
    indtype tmp = X.UB[i] - X.LB[i];
    if(tmp == 0)
    {
      *Xmflsss.hope = X.UB[i];
      // std::cout << "in growTwinForKnapsack(), we get a hope\n";
      Xmflsss.existingProfitSum += f.profitVec[*Xmflsss.hope];
      // std::cout << "Xmflsss.existingProfitSum becomes = " << Xmflsss.existingProfitSum << "\n";
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


  // erase all positions where LB and UB meet.
  indtype &Nzeroed = X.Nzeroed;
  Nzeroed = olpend - overlapPosition;
  if(Nzeroed > 0)
  {
    vec<valtype> cntrS(f.d);
    valtype *S = &*cntrS.begin();
    *olpend = X.len;
    for(indtype i = 0; i < Nzeroed; ++i)
    {
      indtype &st = overlapPosition[i], &end = overlapPosition[i + 1];
      mvalPlus(S, S, f.M[0][X.UB[st]], f.d);
      std::copy(X.LB + st + 1, X.LB + end, X.LB + st - i);
      std::copy(X.UB + st + 1, X.UB + end, X.UB + st - i);
    }
    X.len -= Nzeroed;


    mvalMinus(X.MIN, X.MIN, S + f.dlst, f.dl); // target changes, so MIN and MAX change
    mvalMinus(X.MAX, X.MAX, S + f.dust, f.du);
    mvalMinus(X.sumLB, X.sumLB, S, f.d);
    mvalMinus(X.sumUB, X.sumUB, S, f.d);


    // after erasion, position may change. Adjust position
    {
      indtype tmp = 0;
      for(indtype *i = overlapPosition; i < olpend; ++i)
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


    Ymflsss.existingProfitSum = Xmflsss.existingProfitSum;


    Ymflsss.SKvec.resize(Xmflsss.SKvec.size() - 1);
    Ymflsss.SKback = &Ymflsss.SKvec[1];
    Ymflsss.indvec.assign(Xmflsss.indvec.size() - 2 * (int)Xmflsss.SKvec[0].len, 0);
    Ymflsss.valvec.assign(Xmflsss.valvec.size() - 2 * ((int)f.d + f.dl + f.du), 0);
    Ymflsss.SRVcntr.assign(f.d, 0);
    Ymflsss.result.reserve(7);
  }


  // initialize Ymflsss.SKvec, the first element
  Ymflsss.setAnSK(&Ymflsss.SKvec[0], &Ymflsss.indvec[0],
                  &Ymflsss.valvec[0], Xmflsss.SKback->len);


  // X takes the lower half, Y takes the upper half
  mPAT<valtype, indtype, mk, useBiSearch> &Y = *(Ymflsss.SKback - 1);


  X.beenUpdated = 1;
  Y.beenUpdated = 1;
  std::copy(X.MIN, X.MIN + f.dl + f.du, Y.MIN);


  std::copy(X.sumUB, X.sumUB + f.d, Y.sumUB);
  std::copy(X.UB, X.UB + X.len, Y.UB);
  indtype cap = (X.UB[X.position] + X.LB[X.position]) / 2;
  indtype capResv = cap;
  indtype i = X.position;
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


  ++Xmflsss.SKback;
  return 1;
}




template<typename valtype, typename indtype, bool mk, bool useBiSearch>
inline void mitosisForKnapsack(
    vec<mflsssOBJ<valtype, indtype, mk, useBiSearch> > &descendants,
    shared<valtype, indtype> &f,
    indtype *LB, indtype *UB, valtype *target, valtype *ME,
    int threads, int avgThreadLoad, bool verbose)
{
  int Ndescendants = 1;
  if(threads > 1)
  {
    Ndescendants = 1 << ((int)std::log2(threads * avgThreadLoad) + 1);
  }


  descendants.resize(Ndescendants);
  descendants[0].initializeForKnapsack(&f, target, ME, LB, UB);
  vec<unsigned char> acntr(Ndescendants, 0);
  unsigned char *dead = &*acntr.begin();
  int j = 1;


  while(j < Ndescendants)
  {
    int iend = j;
    for(int i = 0; i < iend; ++i, ++j)
    {
      if(dead[i])
      {
        dead[j] = 1;
        continue;
      }


      indtype boo = growTwinForKnapsack(descendants[i], descendants[j] //, &outfile
      );


      mPAT<valtype, indtype, mk, useBiSearch> &tmp = *descendants[i].SKback;
      if(boo == 0)
      {
        dead[i] = 1;
        dead[j] = 1;
      }
      else if(boo != 1)
      {
        std::copy(tmp.UB, tmp.UB + tmp.len, descendants[i].hope);
        valtype tmpProfit = 0;
        for(indtype k = 0; k < f.subsetSize; ++k)
        {
          tmpProfit += f.profitVec[descendants[i].hopeV[k]];
        }
        if(tmpProfit > f.optimalProfit)
        {
          f.optimalProfit = tmpProfit;
          std::copy(descendants[i].hopeV.begin(),
                    descendants[i].hopeV.end(), f.optimalSolution);
          if(verbose) Rcpp::Rcout << "Updated profit = " << tmpProfit << "\n";
        }
        dead[i] = 1;
        dead[j] = 1;
      }
    }
  }


  // cleansing descendants[]
  int validDescendents = Ndescendants - std::accumulate(dead, dead + Ndescendants, (int)0);
  // std::cout << "validDescendents = " << validDescendents << "\n";
  if(validDescendents == Ndescendants) return;


  vec<mflsssOBJ<valtype, indtype, mk, useBiSearch> > descendantsRemain(validDescendents);
  // std::cout << "validDescendents = " << validDescendents << "\n";
  for(int i = 0, k = 0; i < Ndescendants; ++i)
  {
    if(dead[i]) continue;
    descendants[i].swap(descendantsRemain[k]);
    ++k;
  }
  descendants.swap(descendantsRemain);
  return;
}




template<typename valtype, typename indtype, bool mk, bool useBiSearch>
struct parMflsssOBJforKnapsack: public RcppParallel::Worker // works for mflsssComoPar()
{
  // bool verbose;
  // tbb::spin_mutex *mx;
  std::mutex *mx;
  vec<mflsssOBJ<valtype, indtype, mk, useBiSearch> > &mflsssOBJvec;
  dynamicTasking *dT;
  void operator() (std::size_t st, std::size_t end)
  {
    for(;;)
    {
      std::size_t objI = 0;
      if(!dT->nextTaskID(objI)) break;
      int tmp = mflsssOBJvec[objI].TTTstackRunForKnapsack(mx // , verbose
                                                            );
      if(tmp == -1) break;
    }
  }


  parMflsssOBJforKnapsack(
    vec<mflsssOBJ<valtype, indtype, mk, useBiSearch> > &mflsssOBJvec, std::size_t maxCore):
    mflsssOBJvec(mflsssOBJvec)
  {
    dynamicTasking dt(maxCore, mflsssOBJvec.size());
    dT = &dt;
    // tbb::spin_mutex mtx;
    std::mutex mtx;
    mx = &mtx;
    RcppParallel::parallelFor(0, dT->NofCore, *this);
  }
};
// ============================================================================






