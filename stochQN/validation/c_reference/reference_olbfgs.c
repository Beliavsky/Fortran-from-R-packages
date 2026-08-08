#include <stdio.h>
#include "stochqn.h"
int main(void){
  double x[2]={2.0,-1.0}, grad[2]={0,0}, *req=NULL;
  task_enum task; info_enum info;
  workspace_oLBFGS *w=initialize_oLBFGS(2,4,1.0,0.0,0.0,1,1);
  run_oLBFGS(0.1,x,grad,&req,&task,w,&info);
  while(w->niter<5 || task!=calc_grad){
    grad[0]=req[0]; grad[1]=req[1];
    run_oLBFGS(0.1,x,grad,&req,&task,w,&info);
  }
  printf("%.17g %.17g %zu\n",x[0],x[1],w->bfgs_memory->mem_used);
  dealloc_oLBFGS(w);
  return 0;
}
