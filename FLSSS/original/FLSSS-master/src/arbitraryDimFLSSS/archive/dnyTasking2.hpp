#pragma once
// [[Rcpp::depends(RcppParallel)]]
# include <RcppParallel.h>
# include <atomic>
# include <thread>


struct dynamicTasking
{
  std::size_t jobEnd, grainSize;
  std::atomic<std::size_t> counter;


  // Return false if no job remains.
  bool nextTaskIDrange(std::size_t &begin, std::size_t &end)
  {
    begin = counter.fetch_add(grainSize, std::memory_order_relaxed);
    if(begin >= jobEnd) return false;
    end = std::min(jobEnd, begin + grainSize);
    return true;
  }


  void reset(std::size_t jobBegin, std::size_t jobEnd, std::size_t grainSize)
  {
    counter = jobBegin;
    this->jobEnd = jobEnd;
    this->grainSize = grainSize;
  }


  dynamicTasking(std::size_t jobBegin, std::size_t jobEnd, std::size_t grainSize)
  {
    reset(jobBegin, jobEnd, grainSize);
  }
};


const std::size_t TOTALNTHREAD = std::thread::hardware_concurrency();


// For lambda function or functor style.
// [&](std::size_t i, std::size_t t){}: WHAT SHOULD BE DONE FOR THE i_th JOB
// AND THE THREAD (t_th) THAT RUNS IT. THIS FUNCTION DOES NOT IMPLY THE i_th JOB
// WILL BE OR SHOULD BE RUN BY THE t_th THREAD, BUT RATHER, IT IS
template<typename RunFun, typename BeforeRunFun, typename AfterRunFun>
struct ParaFor: public RcppParallel::Worker
{
  dynamicTasking *dT;
  RunFun run;
  BeforeRunFun beforeRun;
  AfterRunFun afterRun;


  void operator() (std::size_t st, std::size_t end)
  {
    beforeRun(st);
    for(;;)
    {
      std::size_t I, Iend;
      if(!dT->nextTaskIDrange(I, Iend)) break;
      for(; I < Iend; ++I) run(I, st);
    }
    afterRun(st);
  }
  // ###########################################################################
  // begin: index of the first job.
  // end - 1: index of the last job.
  // begin AND end ARE NOT NECESSARILY THE BOUNDS OF SOME CONTAINER!!!
  // RunFun is a lambda function in the form of
  //   [&](std::size_t i, std::size_t t) { define and run the i_th job. }
  // Or, RunFun is a functor with overloaded operator() (std::size_t i, std::size_t t).
  // Similarly, BeforeRunFun is a lambda function in the form of [&](std::size_t t){}
  // where t is the thread ID. This function will be executed by the thread
  // before running the jobs.
  // AfterRunFun is a lambda function in the form of [&](std::size_t t){}
  // where t is the thread ID. This function will be executed by the thread
  // After running the jobs.
  // ###########################################################################
  ParaFor(std::size_t begin, std::size_t end, RunFun run,
          std::size_t maxCore, std::size_t grainSize,
          BeforeRunFun beforeRun, AfterRunFun afterRun):
    run(run), beforeRun(beforeRun), afterRun(afterRun)
  {
    maxCore = std::min<std::size_t>(
      end - begin, std::min<std::size_t> (maxCore, TOTALNTHREAD));
    dynamicTasking dt(begin, end, grainSize); dT = &dt;
    parallelFor(0, maxCore, *this);
  }
};


// DO NOT DELETE!!!
// Example use:
// // [[Rcpp::export]]
// double paraSummation(NumericVector x, int maxCore = 15)
// {
//   NumericVector S(maxCore);
//   ParaFor(0, x.size(), [&](std::size_t i, std::size_t t)
//   {
//     S[t] += x[i];
//   }, maxCore, 128, [](std::size_t t){}, [](std::size_t t){});
//   return std::accumulate(S.begin(), S.end(), 0.0);
// }




// run, beforeRun, afterRun return true if early return is demanded.
template<typename RunFun, typename BeforeRunFun, typename AfterRunFun>
struct ParaForEarlyReturn: public RcppParallel::Worker
{
  dynamicTasking *dT;
  RunFun run;
  BeforeRunFun beforeRun;
  AfterRunFun afterRun;


  void operator() (std::size_t st, std::size_t end)
  {
    bool earlyReturn = false;
    earlyReturn = beforeRun(st);
    if (earlyReturn) return;
    for (; !earlyReturn; )
    {
      std::size_t I, Iend;
      if (!dT->nextTaskIDrange(I, Iend)) break;
      for (; I < Iend and !earlyReturn; ++I) earlyReturn = run(I, st);
    }
    if (!earlyReturn) afterRun(st);
  }
  // ###########################################################################
  // begin: index of the first job.
  // end - 1: index of the last job.
  // begin AND end ARE NOT NECESSARILY THE BOUNDS OF SOME CONTAINER!!!
  // RunFun is a lambda function in the form of
  //   [&](std::size_t i, std::size_t t) { define and run the i_th job. }
  // Or, RunFun is a functor with overloaded operator() (std::size_t i, std::size_t t).
  // Similarly, BeforeRunFun is a lambda function in the form of [&](std::size_t t){}
  // where t is the thread ID. This function will be executed by the thread
  // before running the jobs.
  // AfterRunFun is a lambda function in the form of [&](std::size_t t){}
  // where t is the thread ID. This function will be executed by the thread
  // After running the jobs.
  // ###########################################################################
  ParaForEarlyReturn(std::size_t begin, std::size_t end, RunFun run,
                     std::size_t maxCore, std::size_t grainSize,
                     BeforeRunFun beforeRun, AfterRunFun afterRun):
    run(run), beforeRun(beforeRun), afterRun(afterRun)
  {
    maxCore = std::min<std::size_t>(
      end - begin, std::min<std::size_t> (maxCore, TOTALNTHREAD));
    dynamicTasking dt(begin, end, grainSize); dT = &dt;
    parallelFor(0, maxCore, *this);
  }
};





























