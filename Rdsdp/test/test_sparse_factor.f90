program test_sparse_factor
   use rdsdp, only : dp, dsdp_problem, dsdp_control, dsdp_solution, dsdp_lp_block, dsdp_solve
   use rdsdp_sparse_ldlt, only : sparse_ldlt_cache, sparse_ldlt_solve
   implicit none

   call test_factor_reuse
   call test_solver_sparse_factor
   print *, 'test_sparse_factor: PASS'

contains

   subroutine test_factor_reuse
      integer, parameter :: n=250
      type(sparse_ldlt_cache) :: cache
      real(dp), allocatable :: a(:,:),b(:),x(:)
      real(dp) :: res
      integer :: i
      logical :: ok

      allocate(a(n,n),b(n),x(n)); a=0.0_dp; b=1.0_dp
      do i=1,n
         a(i,i)=4.0_dp
      end do
      do i=1,n-1
         a(i,i+1)=-1.0_dp
         a(i+1,i)=-1.0_dp
      end do
      call sparse_ldlt_solve(a,b,x,1.0e-14_dp,0.0_dp,cache,ok)
      if (.not.ok) error stop 'first sparse LDL solve failed'
      res=maxval(abs(matmul(a,x)-b))
      if (res>2.0e-12_dp) error stop 'first sparse LDL residual too large'
      if (cache%symbolic_analyses/=1) error stop 'expected one symbolic analysis'

      do i=1,n
         a(i,i)=a(i,i)+0.25_dp
      end do
      call sparse_ldlt_solve(a,b,x,1.0e-14_dp,0.0_dp,cache,ok)
      if (.not.ok) error stop 'second sparse LDL solve failed'
      res=maxval(abs(matmul(a,x)-b))
      if (res>2.0e-12_dp) error stop 'second sparse LDL residual too large'
      if (cache%symbolic_analyses/=1) error stop 'symbolic pattern was not reused'
      if (cache%numeric_factorizations/=2) error stop 'numeric factorization count mismatch'
      if (cache%factor_nnz>2*n) error stop 'unexpected fill in tridiagonal factor'
   end subroutine test_factor_reuse

   subroutine test_solver_sparse_factor
      integer, parameter :: n=120
      type(dsdp_problem) :: p
      type(dsdp_control) :: ctrl
      type(dsdp_solution) :: sol
      integer :: i

      p%m=n; allocate(p%b(n),p%block(1)); p%b=0.5_dp
      p%block(1)%category=dsdp_lp_block; p%block(1)%n=n
      allocate(p%block(1)%cdiag(n),p%block(1)%adiag(n,n))
      p%block(1)%cdiag=1.0_dp; p%block(1)%adiag=0.0_dp
      do i=1,n
         p%block(1)%adiag(i,i)=1.0_dp
      end do
      ctrl=dsdp_control()
      ctrl%use_sparse_schur_factor=.true.
      ctrl%sparse_schur_threshold=20
      ctrl%sparse_schur_density_limit=0.10_dp
      call dsdp_solve(p,sol,ctrl)
      if (sol%status/=0) error stop 'sparse-factor LP did not converge'
      if (abs(sol%dobj-60.0_dp)>1.0e-4_dp) error stop 'sparse-factor LP objective mismatch'
      if (sol%sparse_factor_solves<=0) error stop 'sparse-factor path was not exercised'
      if (sol%sparse_symbolic_analyses/=1) error stop 'solver did not reuse symbolic analysis'
      if (sol%direct_schur_solves/=0) error stop 'unexpected dense fallback in sparse-factor LP'
      if (sol%schur_factor_nnz>=sol%schur_matrix_nnz) error stop 'sparse factor did not reduce storage'
   end subroutine test_solver_sparse_factor

end program test_sparse_factor
