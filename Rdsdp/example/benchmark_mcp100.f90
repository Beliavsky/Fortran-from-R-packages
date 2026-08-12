program benchmark_mcp100
   use rdsdp
   implicit none
   type(dsdp_problem) :: p
   type(dsdp_solution) :: ss,sd
   type(dsdp_control) :: cs,cd
   real(dp) :: t0,t1,ts,td

   call read_sdpa('data/mcp100.dat-s',p)
   cs=dsdp_control(); cs%use_sparse_data=.true.
   call cpu_time(t0); call dsdp_solve(p,ss,cs); call cpu_time(t1); ts=t1-t0
   write(*,'(a,f10.4,a,es16.8)') 'sparse data time: ',ts,'  objective: ',ss%dobj
   write(*,'(a,i0,a,i0)') 'stored SDP nnz: ',ss%sdp_data_nnz,'  sparse pair evaluations: ',ss%sparse_pair_evals

   cd=cs; cd%use_sparse_data=.false.
   call cpu_time(t0); call dsdp_solve(p,sd,cd); call cpu_time(t1); td=t1-t0
   write(*,'(a,f10.4,a,es16.8)') 'dense reference time: ',td,'  objective: ',sd%dobj
   if (ts>0.0_dp) write(*,'(a,f10.2)') 'dense/sparse time ratio: ',td/ts
   write(*,'(a,es12.4)') 'objective difference: ',abs(ss%dobj-sd%dobj)
end program benchmark_mcp100
