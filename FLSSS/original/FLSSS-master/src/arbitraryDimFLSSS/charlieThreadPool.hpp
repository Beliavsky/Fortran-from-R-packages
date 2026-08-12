#pragma once
#include <thread>
#include <atomic>


struct DTask
{
  std::size_t jobEnd, grainSize;
  std::atomic<std::size_t> counter;


  // Return false if no job remains.
  bool nextTaskIDrange(std::size_t &begin, std::size_t &end)
  {
    begin = counter.fetch_add(grainSize, std::memory_order_relaxed);
    if (begin >= jobEnd) return false;
    end = std::min(jobEnd, begin + grainSize);
    return true;
  }


  void reset(std::size_t jobBegin, std::size_t jobEnd, std::size_t grainSize)
  {
    counter = jobBegin;
    this->jobEnd = jobEnd;
    this->grainSize = grainSize;
  }


  DTask(){}
  DTask(std::size_t jobBegin, std::size_t jobEnd, std::size_t grainSize)
  {
    reset(jobBegin, jobEnd, grainSize);
  }
};


struct CharlieThreadPool
{
  unsigned maxCore;
  volatile bool *haveFood; // haveFood[maxCore] will be appEnd indicator.
  std::thread *tp;
  std::function<bool(std::size_t, std::size_t)> run;
  std::function<bool(std::size_t)> beforeRun;
  std::function<bool(std::size_t)> afterRun;
  DTask dT; // Will be set by ParaFor object.


  void runJobs(std::size_t threadID) // threadID = 0 is the main thread.
  {
    bool earlyReturn = false;
    earlyReturn = beforeRun(threadID);
    if (earlyReturn) return;
    while (!earlyReturn)
    {
      std::size_t I, Iend;
      if (!dT.nextTaskIDrange(I, Iend)) break;
      for (; I < Iend and !earlyReturn; ++I) earlyReturn = run(I, threadID);
    }
    if (!earlyReturn) afterRun(threadID);
  }


  void live(int threadID)
  {
    while (true)
    {
      while ( !haveFood[threadID] )
      {
        if (haveFood[maxCore]) return; // if (appEnd) return;
      }
      runJobs(threadID);
      haveFood[threadID] = false;
    }
  }


  void initialize(unsigned maxCore)
  {
    maxCore = std::min<unsigned> (
      std::thread::hardware_concurrency(), maxCore);
    this->maxCore = maxCore;
    tp = new std::thread [maxCore];
    haveFood = new volatile bool [maxCore + 1];
    std::fill(haveFood, haveFood + maxCore + 1, false); // appEnd = false;
    for (unsigned i = 1; i < maxCore; ++i) // Fire up all the worker threads.
      tp[i] = std::thread(&CharlieThreadPool::live, this, i);
  }


  void destroy()
  {
    haveFood[maxCore] = true; // appEnd = true;
    for (unsigned i = 1; i < maxCore; ++i) tp[i].join();
    delete [] tp;
    tp = nullptr;
    delete [] haveFood;
    haveFood = nullptr;
  }


  void reset(unsigned maxCore)
  {
    maxCore = std::min<unsigned> (
      std::thread::hardware_concurrency(), maxCore);
    if (maxCore != this->maxCore)
    {
      destroy();
      initialize(maxCore);
    }
  }


  CharlieThreadPool(unsigned maxCore) { initialize(maxCore); }


  ~CharlieThreadPool() { if (haveFood != nullptr) destroy(); }


  void parFor(std::size_t begin, std::size_t end,
              std::function<bool(std::size_t, std::size_t)> run,
              std::size_t grainSize,
              std::function<bool(std::size_t)> beforeRun,
              std::function<bool(std::size_t)> afterRun)
  {
    this->run = run;
    this->beforeRun = beforeRun;
    this->afterRun = afterRun;
    this->dT.reset(begin, end, grainSize);
    std::fill(this->haveFood, this->haveFood + this->maxCore, true); // Kick off job runs.
    this->runJobs(0); // Main thread also runs jobs.
    bool allfinished = false;
    while (!allfinished)
    {
      allfinished = true;
      for (unsigned i = 1; i < this->maxCore; ++i)
        allfinished &= !this->haveFood[i];
    }
  }
};




// DO NOT DELETE!!!
// Example use:
// // [[Rcpp::export]]
// double paraSummation(NumericVector x, int maxCore = 15)
// {
//   NumericVector S(maxCore);
//   CharlieThreadPool ctp(maxCore);
//   ParaFor(0, x.size(), [&](std::size_t i, std::size_t t)
//   {
//     S[t] += x[i];
//   }, maxCore, 128, [](std::size_t t){}, [](std::size_t t){});
//   return std::accumulate(S.begin(), S.end(), 0.0);
// }













