#ifndef maxflow_solver_h
#define maxflow_solver_h

#include "performance.h"
#include "dynamic.h"
#include "logs.h"

class maxflow_solver{
public:
  //dynamic::fixed_array1<double> info;
  struct tinfo{
    int nV;
    int nE;
    int nBE;
    std::string name;
    double flow;
    int sweeps;
    debug::PerformanceCounter construct_t;
    debug::PerformanceCounter solve_t;
    debug::PerformanceCounter msg_t;
    debug::PerformanceCounter augment_t;
    debug::PerformanceCounter relabel_t;
    debug::PerformanceCounter gap_t;
    debug::PerformanceCounter diskr_t;
    debug::PerformanceCounter diskw_t;
    size_t mem_shared;
    size_t mem_region;
    long long diskr_b;
    long long diskw_b;
  }info;
  //
  struct tparams{
    int n_threads;
  }params;
 protected:
  int S,T;
 public:
  typedef float tcap;
  typedef double tflow;
 public:
  maxflow_solver();
  virtual ~maxflow_solver(){};
  virtual tflow maxflow()=0;
  void print_info();//const char * file, bool append = false);
  virtual void save_cut_to_file(const std::string & filename)=0;
  virtual void save_cut()=0;
  virtual void get_cut(int * S){};
  virtual tflow cut_cost(int * S){return 0;};
  virtual tflow cut_cost(){return 0;};
};

#endif
