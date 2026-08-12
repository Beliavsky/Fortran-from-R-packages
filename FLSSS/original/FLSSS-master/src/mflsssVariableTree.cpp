// [[Rcpp::depends(RcppParallel)]]
# include <Rcpp.h>
# include <RcppParallel.h>
# include "header/dnyTasking.hpp"
// # include "header/multiDIO.hpp"
# include "legacy/multiDstack.hpp"
using namespace Rcpp;
using namespace RcppParallel;




/*
template<typename valtype, typename indtype>
struct parMflsss: public Worker
{
  bool useBisearchInFindBounds;
  indtype len, N;
  indtype d; // d is the dimensionality after padding the key column
  indtype dlst, dl, dust, du, keyInd;
  int sizeNeeded;
  // int &firstRunningTask; // record the last + 1 keyTarget right before all threads end
  valtype *originalTarget;
  valtype *keyTarget;
  valtype *scaleFactor;
  valtype ***M, *ME;
  indtype *commonLB, *commonUB;
  std::atomic<int> *totalSize; // initialized as 0
  vec<vec<vec<indtype> > > &result; // result.size() equals the number of threads
  double endTime;
  vec<vec<indtype> > &intCtnrGroup; // each thread takes up a container
  vec<vec<valtype> > &valCtnrGroup;
  vec<vec<mPAT<valtype, indtype> > > &SKgroup;
  dynamicTasking *dT;


  void operator()(std::size_t st, std::size_t end)
  {
    for(;;)
    {
      std::size_t objI = 0;
      if(!dT->nextTaskID(objI)) break;
      // std::cout << (int)objI << ", ";


      // __________________________________________________________________________________________
      // thread function:
      {
        vec<indtype> &intCtnr = intCtnrGroup[st];
        vec<valtype> &valCtnr = valCtnrGroup[st];
        vec<mPAT<valtype, indtype> > &SK = SKgroup[st];
        mPAT<valtype, indtype> *SKbegin = &SK.front();


        // fill the first SK
        {
          SKbegin->beenUpdated = 0;
          SKbegin->MIN = &valCtnr[0];
          SKbegin->MAX = SKbegin->MIN + dl;
          SKbegin->sumLB = SKbegin->MAX + du;
          SKbegin->sumUB = SKbegin->sumLB + d;
          // SKbegin->sumBresv = SKbegin->sumUB + d;


          SKbegin->LB = &intCtnr[0];
          SKbegin->UB = SKbegin->LB + len;
          SKbegin->Bresv = SKbegin->UB + len;
          SKbegin->len = len;


          for(indtype i = 0; i < len; ++i)
          {
            SKbegin->LB[i] = commonLB[i];
            SKbegin->UB[i] = commonUB[i];
          }


          // compute MIN and MAX from target and ME
          {
            valtype target[d];
            // target[0] = keyTarget[objI];
            // for(indtype i = 1; i < d; ++i)
            // {
            //   target[i] = keyTarget[objI] * scaleFactor[i - 1] + originalTarget[i - 1];
            // }


            for(indtype k = 0; k < d; ++k)
            {
              target[k] = keyTarget[objI] * scaleFactor[k] + originalTarget[k];
            }


            // Derive MIN and MAX
            valtype *vst = target + dlst;
            valtype *me = ME + dlst;
            for(indtype i = 0; i < dl; ++i)
            {
              SKbegin->MIN[i] = vst[i] - me[i];
              // std::cout << SKbegin->MIN[i] << ", ";
            }
            // std::cout << "\n";
            vst = target + dust;
            me = ME + dust;
            for(indtype i = 0; i < du; ++i)
            {
              SKbegin->MAX[i] = vst[i] + me[i];
              // std::cout << SKbegin->MAX[i] << ", ";
            }
            // std::cout << "------\n";
          }


          // SKbegin->target[0] = keyTarget[objI];
          // for(indtype i = 1; i < d; ++i)
          // {
          //   SKbegin->target[i] = keyTarget[objI] * scaleFactor[i - 1] + originalTarget[i - 1];
          // }


          std::fill(SKbegin->sumLB, SKbegin->sumLB + d, 0);
          iterSum <valtype, indtype> (SKbegin->sumLB, M[0], SKbegin->LB, len, d);
          std::fill(SKbegin->sumUB, SKbegin->sumUB + d, 0);
          iterSum <valtype, indtype> (SKbegin->sumUB, M[0], SKbegin->UB, len, d);
        }


        mPAT<valtype, indtype> *SKback = SKbegin + 1;
        vec<vec<indtype> > &rst = result[st];


        // TTTstack<valtype, indtype> (
        //     len, N, d, M, ME, rst, sizeNeeded,
        //     SKbegin, SKback, useBisearchInFindBounds, *totalSize, endTime);


        TTTstack<valtype, indtype> (
            len, N, d, dlst, dl, dust, du, M, rst, sizeNeeded,
            SKbegin, SKback, useBisearchInFindBounds, *totalSize, endTime);


        // std::cout << "(int)*totalSize = " << (int)*totalSize << ", ";
        // std::cout << "sizeNeeded = " << sizeNeeded << ", ";
        // std::cout << "rst.size() = " << rst.size() << "\n";
        if(*totalSize >= sizeNeeded or (double)std::clock() > endTime) break;
      }
      // __________________________________________________________________________________________
    }
  }




  parMflsss(bool useBisearchInFindBounds, indtype len, indtype N,
            indtype d, indtype dlst, indtype dl, indtype dust, indtype du, indtype keyInd,
            int sizeNeeded, // int &firstRunningTask,
            valtype *originalTarget, valtype *keyTarget,
            valtype *scaleFactor,
            valtype ***M, valtype *ME, indtype *commonLB, indtype *commonUB,
            vec<vec<vec<indtype> > > &result, double endTime,
            vec<vec<indtype> > &intCtnrGroup,
            vec<vec<valtype> > &valCtnrGroup,
            vec<vec<mPAT<valtype, indtype> > > &SKgroup,
            int maxThreads, int keyTargetSize):
    useBisearchInFindBounds(useBisearchInFindBounds),
    len(len), N(N),
    d(d), dlst(dlst), dl(dl), dust(dust), du(du), keyInd(keyInd),
    sizeNeeded(sizeNeeded), // firstRunningTask(firstRunningTask),
    originalTarget(originalTarget), keyTarget(keyTarget), scaleFactor(scaleFactor),
    M(M), ME(ME), commonLB(commonLB), commonUB(commonUB),
    result(result), endTime(endTime),
    intCtnrGroup(intCtnrGroup), valCtnrGroup(valCtnrGroup),
    SKgroup(SKgroup)
  {
    std::atomic<int> tsize = 0;
    totalSize = &tsize;
    dynamicTasking dtask(maxThreads, keyTargetSize);
    dT = &dtask;
      parallelFor(0, dT->NofCore, *this);
  }
};
*/




// ________________________________________________________________________________________________
template<typename valtype, typename indtype>
struct parMflsss: public Worker
{
  bool useBisearchInFindBounds;
  indtype len, N;
  indtype d; // d is the dimensionality after padding the key column
  indtype dlst, dl, dust, du, keyInd;
  int sizeNeeded;
  // int &firstRunningTask; // record the last + 1 keyTarget right before all threads end
  valtype *originalTarget;
  valtype *keyTarget;
  valtype *scaleFactor;
  valtype ***M, *ME;
  indtype *commonLB, *commonUB;
  std::atomic<int> *totalSize; // initialized as 0
  vec<vec<vec<indtype> > > &result; // result.size() equals the number of threads
  double endTime;
  vec<vec<indtype> > &intCtnrGroup; // each thread takes up a container
  vec<vec<valtype> > &valCtnrGroup;
  vec<vec<mPAT<valtype, indtype> > > &SKgroup;
  dynamicTasking *dT;


  void operator()(std::size_t st, std::size_t end)
  {
    for(;;)
    {
      std::size_t objI = 0;
      if(!dT->nextTaskID(objI)) break;
      // std::cout << (int)objI << ", ";


      // __________________________________________________________________________________________
      // thread function:
      {
        vec<indtype> &intCtnr = intCtnrGroup[st];
        vec<valtype> &valCtnr = valCtnrGroup[st];
        vec<mPAT<valtype, indtype> > &SK = SKgroup[st];
        mPAT<valtype, indtype> *SKbegin = &SK.front();


        // fill the first SK
        {
          // SKbegin->beenUpdated = 0;
          SKbegin->MIN = &valCtnr[0];
          SKbegin->MAX = SKbegin->MIN + dl;
          SKbegin->sumLB = SKbegin->MAX + du;
          SKbegin->sumUB = SKbegin->sumLB + d;
          // SKbegin->sumBresv = SKbegin->sumUB + d;


          SKbegin->LB = &intCtnr[0];
          SKbegin->UB = SKbegin->LB + len;
          SKbegin->Bresv = SKbegin->UB + len;
          SKbegin->len = len;


          for(indtype i = 0; i < len; ++i)
          {
            SKbegin->LB[i] = commonLB[i];
            SKbegin->UB[i] = commonUB[i];
          }


          // compute MIN and MAX from target and ME
          {
            // valtype target[d];
            vec<valtype> targetContainer(d);
            valtype *target = &targetContainer[0];
            for(indtype k = 0; k < d; ++k)
            {
              target[k] = keyTarget[objI] * scaleFactor[k] + originalTarget[k];
            }


            // Derive MIN and MAX
            valtype *vst = target + dlst;
            valtype *me = ME + dlst;
            for(indtype i = 0; i < dl; ++i)
            {
              SKbegin->MIN[i] = vst[i] - me[i];
            }
            vst = target + dust;
            me = ME + dust;
            for(indtype i = 0; i < du; ++i)
            {
              SKbegin->MAX[i] = vst[i] + me[i];
            }
          }


          std::fill(SKbegin->sumLB, SKbegin->sumLB + d, 0);
          iterSum <valtype, indtype> (SKbegin->sumLB, M[0], SKbegin->LB, len, d);
          std::fill(SKbegin->sumUB, SKbegin->sumUB + d, 0);
          iterSum <valtype, indtype> (SKbegin->sumUB, M[0], SKbegin->UB, len, d);
        }


        mPAT<valtype, indtype> *SKback = SKbegin + 1;
        vec<vec<indtype> > &rst = result[st];


        TTTstack<valtype, indtype> (
            len, N, d, dlst, dl, dust, du, M, rst, sizeNeeded,
            SKbegin, SKback, useBisearchInFindBounds, *totalSize, endTime);


        if(*totalSize >= sizeNeeded or (double)std::clock() > endTime) break;
      }
      // __________________________________________________________________________________________
    }
  }




  parMflsss(bool useBisearchInFindBounds, indtype len, indtype N,
            indtype d, indtype dlst, indtype dl, indtype dust, indtype du, indtype keyInd,
            int sizeNeeded, // int &firstRunningTask,
            valtype *originalTarget, valtype *keyTarget,
            valtype *scaleFactor,
            valtype ***M, valtype *ME, indtype *commonLB, indtype *commonUB,
            vec<vec<vec<indtype> > > &result, double endTime,
            vec<vec<indtype> > &intCtnrGroup,
            vec<vec<valtype> > &valCtnrGroup,
            vec<vec<mPAT<valtype, indtype> > > &SKgroup,
            std::size_t maxThreads, int keyTargetSize):
    useBisearchInFindBounds(useBisearchInFindBounds),
    len(len), N(N),
    d(d), dlst(dlst), dl(dl), dust(dust), du(du), keyInd(keyInd),
    sizeNeeded(sizeNeeded),
    originalTarget(originalTarget), keyTarget(keyTarget), scaleFactor(scaleFactor),
    M(M), ME(ME), commonLB(commonLB), commonUB(commonUB),
    result(result), endTime(endTime),
    intCtnrGroup(intCtnrGroup), valCtnrGroup(valCtnrGroup),
    SKgroup(SKgroup)
  {
    // std::atomic<int> tsize = 0;
    std::atomic<int> tsize(0);
    totalSize = &tsize;
    dynamicTasking dtask(maxThreads, keyTargetSize);
    dT = &dtask;
      parallelFor(0, dT->NofCore, *this);
  }
};








template<typename valtype, typename indtype>
// vr has N rows and _d columns
List mFLSSScpp(
    int len, List vr,
    int _d, int dlst, int dl, int dust, int du, int keyInd, int N,
    NumericVector scaleFactorr,
    NumericVector originalTargetr,
    NumericVector keyTargetr,
    NumericVector MEr,
    IntegerVector LBr, IntegerVector UBr,
    int sizeNeeded, double endTime, int maxCore,
    bool useBiSearchInFindBounds)
{
  triM<valtype, indtype> mat;
  std::size_t byteSize = mat.containerByteSize(_d, N, len);
  vec<unsigned char> matContainer(byteSize);
  mat.alloc(&matContainer[0], _d, len, N);
  mat.make(&matContainer[0], len, vr);
  valtype ***M = mat.mat;


  valtype *ME = (valtype*)&*MEr.begin();
  vec<valtype> MEcontain;
  if(sizeof(valtype) < 8)
  {
    MEcontain.assign(MEr.begin(), MEr.end());
    ME = &*MEcontain.begin();
  }


  valtype *originalTarget = (valtype*)&*originalTargetr.begin();
  vec<valtype> originalTargetContain;
  if(sizeof(valtype) < 8)
  {
    originalTargetContain.assign(originalTargetr.begin(), originalTargetr.end());
    originalTarget = &*originalTargetContain.begin();
  }


  valtype *keyTarget = (valtype*)&*keyTargetr.begin();
  vec<valtype> keyTargetContain;
  if(sizeof(valtype) < 8)
  {
    keyTargetContain.assign(keyTargetr.begin(), keyTargetr.end());
    keyTarget = &*keyTargetContain.begin();
  }


  valtype *scaleFactor = (valtype*)&*scaleFactorr.begin();
  vec<valtype> scaleFactorContain;
  if(sizeof(valtype) < 8)
  {
    scaleFactorContain.assign(scaleFactorr.begin(), scaleFactorr.end());
    scaleFactor = &*scaleFactorContain.begin();
  }


  vec<indtype> iniContain(int(2) * len);
  indtype *commonLB = &iniContain[0], *commonUB = commonLB + len;
  // indtype commonLB[len], commonUB[len];
  for(indtype i = 0; i < len; ++i)
  {
    commonLB[i] = LBr[i] - 1;
    commonUB[i] = UBr[i] - 1;
  }


  vec<vec<vec<indtype> > > result(maxCore);
  vec<vec<indtype> > intCtnrGroup(maxCore);
  vec<vec<valtype> > valCtnrGroup(maxCore);
  vec<vec<mPAT<valtype, indtype> > > SKgroup(maxCore);
  {
    std::size_t stackLen = (unsigned)len + 3;
    // unsigned biscaleFactor = (unsigned)std::log2(N + 0.0 - (unsigned)len) + 1;
    unsigned biscaleFactor = 1;
    for(int i = 0; i < maxCore; ++i)
    {
      intCtnrGroup[i].resize(stackLen * (stackLen + 1) / 2 * 3 * biscaleFactor);
      valCtnrGroup[i].resize(  stackLen * ((std::size_t)_d * 3 +
        (std::size_t)dl + (std::size_t)du) * biscaleFactor  );
      // contain 4 values: MIN, MAX, sumLB, sumUB. MIN and MAX can be less than _d
      SKgroup[i].resize(stackLen * biscaleFactor);
      result[i].reserve(sizeNeeded);
    }
  }


  // int firstRunningTask = 0;
  parMflsss<valtype, indtype> (
      useBiSearchInFindBounds, len, N,
      _d, dlst, dl, dust, du, keyInd,
      sizeNeeded, // firstRunningTask,
      originalTarget, keyTarget, scaleFactor,
      M, ME, commonLB, commonUB, result,
      endTime, intCtnrGroup, valCtnrGroup, SKgroup,
      maxCore, keyTargetr.size());


  int solutionN = 0;
  {
    for(int i = 0, iend = result.size(); i < iend; ++i)
    {
      solutionN += result[i].size();
    }
  }


  List lis(solutionN);
  for(int k = 0, i = 0, iend = result.size(); i < iend; ++i)
  {
    for(int j = 0, jend = result[i].size(); j < jend; ++j)
    {
      IntegerVector tmp(result[i][j].size());
      for(int u = 0, uend = result[i][j].size(); u < uend; ++u)
      {
        tmp[u] = result[i][j][u] + 1;
      }
      lis[k] = tmp;
      ++k;
    }
  }


  // List memImage = exportAll <valtype, indtype> (
  //   _d, len, N, firstRunningTask, ME, matContainer, keyTargetr, scaleFactorr,
  //   originalTargetr, lis, useBiSearchInFindBounds);


  return lis;
}




// [[Rcpp::export]]
List z_mFLSSSvariableTree(int maxCore, int len, List vr,
              int d, int dlst, int dl, int dust, int du, int keyInd,
              NumericVector originalTarget, NumericVector keyTarget,
              NumericVector scaleFactor, NumericVector MEr, IntegerVector LBr, IntegerVector UBr,
              int sizeNeed, double duration, bool useFloat, bool useBisearchInFindBounds = 0)
{
  NumericVector tmp = vr[0];
  int vlen = tmp.size();
  double endTime = (double)std::clock() + duration * CLOCKS_PER_SEC;
  List result;


  if(std::max(vlen, d) < 127)
  {
    if(useFloat) result = mFLSSScpp <float, signed char> (len, vr, d, dlst, dl, dust, du, keyInd, vlen, scaleFactor, originalTarget, keyTarget, MEr, LBr, UBr, sizeNeed, endTime, maxCore, useBisearchInFindBounds);
    else result = mFLSSScpp <double, signed char> (len, vr, d, dlst, dl, dust, du, keyInd, vlen, scaleFactor, originalTarget, keyTarget, MEr, LBr, UBr, sizeNeed, endTime, maxCore, useBisearchInFindBounds);
  }
  else if(std::max(vlen, d) < 32767)
  {
    if(useFloat) result = mFLSSScpp <float, short> (len, vr, d, dlst, dl, dust, du, keyInd, vlen, scaleFactor, originalTarget, keyTarget, MEr, LBr, UBr, sizeNeed, endTime, maxCore, useBisearchInFindBounds);
    else result = mFLSSScpp <double, short> (len, vr, d, dlst, dl, dust, du, keyInd, vlen, scaleFactor, originalTarget, keyTarget, MEr, LBr, UBr, sizeNeed, endTime, maxCore, useBisearchInFindBounds);
  }
  else
  {
    if(useFloat) result = mFLSSScpp <float, int> (len, vr, d, dlst, dl, dust, du, keyInd, vlen, scaleFactor, originalTarget, keyTarget, MEr, LBr, UBr, sizeNeed, endTime, maxCore, useBisearchInFindBounds);
    else result = mFLSSScpp <double, int> (len, vr, d, dlst, dl, dust, du, keyInd, vlen, scaleFactor, originalTarget, keyTarget, MEr, LBr, UBr, sizeNeed, endTime, maxCore, useBisearchInFindBounds);
  }
  return result;
}




/*
template<typename valtype, typename indtype>
// vr has N rows and _d columns
List mFLSSSimportCpp(
    List memImage,
    int len, List vr, int _d, int N,
    NumericVector scaleFactorr,
    NumericVector originalTargetr,
    NumericVector keyTargetr,
    NumericVector MEr,
    IntegerVector LBr, IntegerVector UBr,
    int sizeNeeded, double endTime, int maxCore,
    bool useBisearchInFindBounds)
{
  triM<valtype, indtype> mat;
  std::size_t byteSize = mat.containerByteSize(_d, N, len);
  vec<unsigned char> matContainer(byteSize);
  mat.alloc(&matContainer[0], _d, len, N);
  mat.make(&matContainer[0], len, vr);
  valtype ***M = mat.mat;


  valtype *ME = (valtype*)&*MEr.begin();
  vec<valtype> MEcontain;
  if(sizeof(valtype) < 8)
  {
    MEcontain.assign(MEr.begin(), MEr.end());
    ME = &*MEcontain.begin();
  }


  valtype *originalTarget = (valtype*)&*originalTargetr.begin();
  vec<valtype> originalTargetContain;
  if(sizeof(valtype) < 8)
  {
    originalTargetContain.assign(originalTargetr.begin(), originalTargetr.end());
    originalTarget = &*originalTargetContain.begin();
  }


  valtype *keyTarget = (valtype*)&*keyTargetr.begin();
  vec<valtype> keyTargetContain;
  if(sizeof(valtype) < 8)
  {
    keyTargetContain.assign(keyTargetr.begin(), keyTargetr.end());
    keyTarget = &*keyTargetContain.begin();
  }


  valtype *scaleFactor = (valtype*)&*scaleFactorr.begin();
  vec<valtype> scaleFactorContain;
  if(sizeof(valtype) < 8)
  {
    scaleFactorContain.assign(scaleFactorr.begin(), scaleFactorr.end());
    scaleFactor = &*scaleFactorContain.begin();
  }


  indtype commonLB[len], commonUB[len];
  for(indtype i = 0; i < len; ++i)
  {
    commonLB[i] = LBr[i] - 1;
    commonUB[i] = UBr[i] - 1;
  }


  vec<vec<vec<indtype> > > result(maxCore);
  vec<vec<indtype> > intCtnrGroup(maxCore);
  vec<vec<valtype> > valCtnrGroup(maxCore);
  vec<vec<mPAT<valtype, indtype> > > SKgroup(maxCore);
  vec<std::size_t> SKbackGroup(maxCore, 1);
  {
    std::size_t stackLen = len + 2, valIntSize = stackLen * (stackLen + 1) / 2 * 3;
    for(int i = 0; i < maxCore; ++i)
    {
      intCtnrGroup[i].resize(valIntSize);
      valCtnrGroup[i].resize(valIntSize);
      SKgroup[i].resize(len + 6);
      result[i].reserve(sizeNeeded);
    }
  }


  parMflsss<valtype, indtype> (
      useBisearchInFindBounds, len, N, _d, sizeNeeded,
      originalTarget, keyTarget, scaleFactor,
      M, ME, commonLB, commonUB, result,
      endTime, intCtnrGroup, valCtnrGroup, SKgroup, SKbackGroup,
      maxCore, keyTargetr.size());


  int solutionN = 0;
  {
    for(int i = 0, iend = result.size(); i < iend; ++i)
    {
      solutionN += result[i].size();
    }
  }


  if(solutionN == 0) return List::create();


  Rcpp::List memImage = exportAll <valtype, indtype> (
    _d, len, N, SKbackGroup, ME, matContainer, intCtnrGroup,
    valCtnrGroup, SKgroup, useBisearchInFindBounds);


  List lis(solutionN);
  for(int k = 0, i = 0, iend = result.size(); i < iend; ++i)
  {
    for(int j = 0, jend = result[i].size(); j < jend; ++j)
    {
      IntegerVector tmp(result[i][j].size());
      for(int u = 0, uend = result[i][j].size(); u < uend; ++u)
      {
        tmp[u] = result[i][j][u] + 1;
      }
      lis[k] = tmp;
      ++k;
    }
  }


  return List::create(Named("solution") = lis, Named("memoryImage") = memImage);
}
*/










































