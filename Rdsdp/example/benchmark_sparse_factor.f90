program benchmark_sparse_factor
   use rdsdp_kinds, only : dp
   use rdsdp_linalg, only : solve_spd
   use rdsdp_sparse_ldlt, only : sparse_ldlt_cache, sparse_ldlt_solve
   implicit none
   integer, parameter :: n=1000
   real(dp) :: a(n,n),b(n),xs(n),xd(n),t0,t1,ts,td,rs,rd
   type(sparse_ldlt_cache) :: cache
   integer :: i
   logical :: ok

   a=0.0_dp; b=1.0_dp
   do i=1,n
      a(i,i)=4.0_dp
   end do
   do i=1,n-1
      a(i,i+1)=-1.0_dp; a(i+1,i)=-1.0_dp
   end do
   call cpu_time(t0); call sparse_ldlt_solve(a,b,xs,1.0e-14_dp,0.0_dp,cache,ok); call cpu_time(t1)
   if (.not.ok) error stop 'sparse factorization failed'
   ts=t1-t0; rs=maxval(abs(matmul(a,xs)-b))
   call cpu_time(t0); call solve_spd(a,b,xd,1.0e-14_dp,ok); call cpu_time(t1)
   if (.not.ok) error stop 'dense Cholesky failed'
   td=t1-t0; rd=maxval(abs(matmul(a,xd)-b))
   write(*,'("sparse RCM+LDL: ",f9.5," s  factor nnz=",i0," residual=",es10.2)') &
      ts,cache%factor_nnz,rs
   write(*,'("dense Cholesky:  ",f9.5," s  dense entries=",i0," residual=",es10.2)') td,n*n,rd
   if (ts>0.0_dp) write(*,'("time ratio dense/sparse: ",f8.2)') td/ts
end program benchmark_sparse_factor
