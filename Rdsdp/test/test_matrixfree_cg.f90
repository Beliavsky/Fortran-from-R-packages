program test_matrixfree_cg
   use rdsdp, only : dp, dsdp_problem, dsdp_control, dsdp_solution, dsdp_lp_block, dsdp_solve
   implicit none
   integer, parameter :: n=48
   type(dsdp_problem) :: p
   type(dsdp_control) :: ctrl
   type(dsdp_solution) :: sol
   real(dp) :: pi,scale
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

   ctrl=dsdp_control()
   ctrl%use_cg=.true.; ctrl%cg_matrix_free=.true.; ctrl%cg_threshold=1
   ctrl%cg_tol=1.0e-10_dp; ctrl%cg_maxiter=100
   ctrl%gaptol=1.0e-6_dp; ctrl%pinfeastol=1.0e-6_dp
   call dsdp_solve(p,sol,ctrl)

   if (sol%status/=0) error stop 'matrix-free PCG LP did not converge'
   if (abs(sol%dobj-real(n,dp))>2.0e-4_dp) error stop 'matrix-free PCG objective mismatch'
   if (sol%matrix_free_cg_solves<=0 .or. sol%matrix_free_matvecs<=0) &
      error stop 'matrix-free PCG path was not exercised'
   if (sol%direct_schur_solves/=0) error stop 'matrix-free PCG unexpectedly fell back to dense solve'
   print *, 'test_matrixfree_cg: PASS'
end program test_matrixfree_cg
