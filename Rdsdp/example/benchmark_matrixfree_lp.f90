program benchmark_matrixfree_lp
   use rdsdp, only : dp, dsdp_problem, dsdp_control, dsdp_solution, dsdp_lp_block, dsdp_solve
   implicit none
   integer, parameter :: n=300
   type(dsdp_problem) :: p
   type(dsdp_control) :: direct_ctrl,mf_ctrl
   type(dsdp_solution) :: direct_sol,mf_sol
   real(dp) :: pi,scale,t0,t1,t_direct,t_mf
   integer :: i,j

   pi=acos(-1.0_dp); scale=sqrt(2.0_dp/real(n+1,dp))
   p%m=n; allocate(p%b(n),p%block(1))
   p%block(1)%category=dsdp_lp_block; p%block(1)%n=n
   allocate(p%block(1)%cdiag(n),p%block(1)%adiag(n,n))
   p%block(1)%cdiag=2.0_dp
   do j=1,n
      do i=1,n
         p%block(1)%adiag(i,j)=scale*sin(real(i*j,dp)*pi/real(n+1,dp))
      end do
   end do
   p%b=0.5_dp*matmul(transpose(p%block(1)%adiag),spread(1.0_dp,1,n))

   direct_ctrl=dsdp_control()
   direct_ctrl%use_sparse_schur_factor=.false.
   direct_ctrl%gaptol=1.0e-6_dp; direct_ctrl%pinfeastol=1.0e-6_dp
   mf_ctrl=direct_ctrl
   mf_ctrl%use_cg=.true.; mf_ctrl%cg_matrix_free=.true.; mf_ctrl%cg_threshold=1
   mf_ctrl%cg_tol=1.0e-9_dp; mf_ctrl%cg_maxiter=100

   call cpu_time(t0); call dsdp_solve(p,direct_sol,direct_ctrl); call cpu_time(t1); t_direct=t1-t0
   call cpu_time(t0); call dsdp_solve(p,mf_sol,mf_ctrl); call cpu_time(t1); t_mf=t1-t0

   write(*,'("assembled direct: ",f9.4," s  objective=",f14.6)') t_direct,direct_sol%dobj
   write(*,'("matrix-free PCG: ",f9.4," s  objective=",f14.6,"  matvecs=",i0)') &
      t_mf,mf_sol%dobj,mf_sol%matrix_free_matvecs
   if (t_mf>0.0_dp) write(*,'("time ratio direct/matrix-free: ",f8.2)') t_direct/t_mf
end program benchmark_matrixfree_lp
