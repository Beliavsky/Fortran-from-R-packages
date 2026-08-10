program test_lanczos_linesearch
   use rcsdp_kinds, only : dp
   use rcsdp_types, only : csdp_block_matrix, csdp_matrix
   use rcsdp_block_ops, only : allocate_block_matrix, line_search_pd
   implicit none
   type(csdp_block_matrix) :: x, dx
   real(dp) :: a_exact, a_lanczos, lambda
   integer :: i, n
   logical :: used

   n=220
   call allocate_block_matrix([csdp_matrix],[n],x)
   call allocate_block_matrix([csdp_matrix],[n],dx)
   x%block(1)%mat=0.0_dp
   dx%block(1)%mat=0.0_dp
   do i=1,n
      x%block(1)%mat(i,i)=1.0_dp
      lambda=0.1_dp+1.9_dp*real(i-1,dp)/real(n-1,dp)
      dx%block(1)%mat(i,i)=-lambda
   end do

   a_exact=line_search_pd(dx,x,0.97_dp,1.0_dp,.false.)
   a_lanczos=line_search_pd(dx,x,0.97_dp,1.0_dp,.true.,180,30,used)
   if (.not.used) error stop 'Lanczos path was not exercised'
   if (abs(a_exact-0.485_dp) > 1.0e-12_dp) error stop 'exact line search reference mismatch'
   if (a_lanczos > a_exact*(1.0_dp+1.0e-10_dp)) error stop 'Lanczos step was not conservative'
   if (abs(a_lanczos-a_exact) > 0.01_dp*a_exact) error stop 'Lanczos line search mismatch'
   print *, 'test_lanczos_linesearch: PASS'
end program test_lanczos_linesearch
