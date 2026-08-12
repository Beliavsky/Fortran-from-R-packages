# pragma once
// // [[Rcpp::depends(RcppParallel)]]
// # include <Rcpp.h>
#include <atomic>
// #include <RcppParallel.h>
// using namespace Rcpp;
// using namespace RcppParallel;


struct dynamicTasking
{
  std::size_t NofCore;
  std::size_t NofAtom;
  // tbb::atomic<std::size_t> counter;
  std::atomic<std::size_t> counter;


  void reset(std::size_t &NofCPU, std::size_t NofTask)
  {
    NofCore = NofCPU;
    if(NofCore > NofTask) NofCore = NofTask;
    NofAtom = NofTask;
    counter = 0;
  }


  dynamicTasking(std::size_t &NofCPU, std::size_t NofTask)
  {
    reset(NofCPU, NofTask);
  }


  bool nextTaskID(std::size_t &taskID)
  {
    // taskID = counter.fetch_and_increment();
	taskID = counter.fetch_add(1);
    return taskID < NofAtom;
  }


  bool nextTaskID(std::size_t &taskID, std::size_t increment)
  {
    // taskID = counter.fetch_and_add(increment);
	taskID = counter.fetch_add(increment);
    return taskID < NofAtom;
  }
};









