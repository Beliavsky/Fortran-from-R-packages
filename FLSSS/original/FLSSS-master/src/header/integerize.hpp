# include "macros.hpp"


/*
// integerize a vector
template<typename valtype, typename indtype>
inline void integerize(INT *result, valtype *v, indtype len, indtype N,
                valtype target, valtype me,
                int &targetInt, int &meINT //, unsigned &digits
                , int precisionLevel)
{
  // std::cout << "\n\n\n";
  // std::cout << me << "\n";
  double r = 0.5 / me;
  target = target * r;
  double sh = (INT)target + 1 - target, avgsh = sh / len;
  target = (INT)target + 1;
  double *rst = (double*)result;
  for(indtype i = 0; i < N; ++i)
  {
    rst[i] = v[i] * r + avgsh;
    rst[i] -= (INT)rst[i];
  }
  // vec<valtype> vsort(v, v + N);
  // std::sort(vsort.begin(), vsort.end());
  std::nth_element(rst, rst + N - len, rst + N);
  valtype tmp = std::accumulate(rst + N - len, rst + N, 0.0);
  int mulfactor = 1 << (unsigned(std::log2(tmp)) + 1 + 1 + precisionLevel);
  // unsigned(std::log2(tmp)) + 1 is the digit that decimals would reach
  // unsigned(std::log2(tmp)) + 1 + 1 is the digit that decimals could affect
  // unsigned(std::log2(tmp)) + 1 + 1 + 1 is even safer


  targetInt = target * mulfactor;
  meINT = mulfactor / 2;
  for(indtype i = 0; i < N; ++i)
  {
    // std::cout << v[i] << ", ";
    rst[i] = (v[i] * r + avgsh) * mulfactor;
    // std::cout << rst[i] << ", ";
    result[i] = std::round(rst[i]);
    // std::cout << result[i] << "        ";
  }
}
*/


template<typename valtype, typename indtype>
inline void integerize(INT *result, valtype *v, indtype len, indtype N,
                       valtype target, valtype me,
                       int &targetInt, int &meINT, int precisionLevel)
{
  double r = 1 / me;
  target = target * r;
  // double sh = (INT)target + 1 - target, avgsh = sh / len;
  // target = (INT)target + 1;
  double *vnew = (double*)result;
  for(indtype i = 0; i < N; ++i)
  {
    // vnew[i] = v[i] * r + avgsh;
    vnew[i] = v[i] * r;
  }


  if(precisionLevel == 0)
  {
    int ratio = 1;
    valtype theMax = *std::max_element(vnew, vnew + N);
    int N8 = N * 8;
    while(theMax * ratio < N8) ratio *= 2; // Every 8 has an element on average
    for(indtype i = 0; i < N; ++i)
    {
      result[i] = std::round(vnew[i] * ratio);
    }
    // targetInt = target * ratio;
    targetInt = std::round(target * ratio);
    meINT = ratio;
    return;
  }


  if(precisionLevel == -1)
  {
    // valtype vsort[N];
    vec<valtype> acntr(N); valtype *vsort = &*acntr.begin();
    std::copy(vnew, vnew + N, vsort);
    std::sort(vsort, vsort + N);
    indtype leastDiffi = 0;
    double currentMax = 1e308;
    for(indtype i = 1; i < N; ++i)
    {
      valtype tmp = vsort[i] - vsort[i - 1];
      if(tmp > 1e-10 and tmp < currentMax)
      {
        leastDiffi = i;
        currentMax = tmp;
      }
    }
    valtype tmp = vsort[leastDiffi] - vsort[leastDiffi - 1];
    int ratio = 1;
    while(INT(std::round(tmp * ratio)) < 1) ratio *= 2;


    for(indtype i = 0; i < N; ++i)
    {
      result[i] = std::round(vnew[i] * ratio);
    }
    targetInt = std::round(target * ratio);
    meINT = ratio;
    return;
  }


  int ratio = 1;
  valtype theMax = *std::max_element(vnew, vnew + N);
  while(theMax * ratio < precisionLevel) ratio *= 2;
  for(indtype i = 0; i < N; ++i)
  {
    result[i] = std::round(vnew[i] * ratio);
  }
  targetInt = std::round(target * ratio);
  meINT = ratio;
}




// columnMajorV has been shifted. Each column's minimum is 0
// integerize a numeric matrix
template<typename valtype, typename indtype>
inline void integerize(INT *rst, int *targetInt, int *meInt,
                valtype *columnMajorV,
                indtype len, indtype N, indtype d,
                valtype *target, valtype *me, int *precisionLevel)
{
  for(int i = 0, iend = d; i < iend; ++i)
  {
    integerize(rst + i * N, columnMajorV + i * N, len, N,
               target[i], me[i], targetInt[i],
               meInt[i], precisionLevel[i]);
  }
}



