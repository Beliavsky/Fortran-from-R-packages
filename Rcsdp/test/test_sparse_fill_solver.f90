program test_sparse_fill_solver
   use rcsdp
   implicit none
   integer, parameter :: n=120
   type(csdp_problem) :: p
   type(csdp_solution) :: s
   type(csdp_control) :: ctrl
   real(dp), allocatable :: cm(:,:)
   integer :: ii(1), jj(1), k
   real(dp) :: vv(1)

   ! A deliberately sparse-fill SDP.  Maximize -trace(X) subject to
   ! X(1,1)=1 and X >= 0.  The optimum is -1 and the makefill pattern is
   ! only the diagonal (120 of 14400 matrix positions).
   call init_problem(p,[csdp_matrix],[n],1)
   allocate(cm(n,n)); cm=0.0_dp
   do k=1,n
      cm(k,k)=-1.0_dp
   end do
   call set_c_matrix_block(p,1,cm)
   ii=[1]; jj=[1]; vv=[1.0_dp]
   call set_sparse_a_block(p,1,1,ii,jj,vv)
   p%b=[1.0_dp]

   ctrl%printlevel=0
   ctrl%maxiter=80
   call csdp(p,s,ctrl)

   if (s%status/=csdp_success) error stop 'sparse-fill solver did not converge'
   if (abs(s%pobj+1.0_dp)>2.0e-7_dp) error stop 'sparse-fill objective mismatch'
   if (s%fill_nnz/=n .or. s%fill_full_entries/=n*n) error stop 'sparse-fill pattern mismatch'
   if (s%fill_sparse_products<=0) error stop 'integrated sparse fill kernels were not exercised'
   if (s%fill_dense_products/=0) error stop 'unexpected dense fill-product fallback'
   print *, 'test_sparse_fill_solver: PASS'
end program test_sparse_fill_solver
