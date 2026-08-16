program basic
   use discrete_laplace
   use, intrinsic :: iso_fortran_env, only : int64
   implicit none
   integer(int64) :: x(12)
   integer :: xi(12),i,st
   real(dp) :: par(2)
   call set_rng_seed(123)
   call rdlaplace2(x,0.4_dp,0.3_dp)
   do i=1,size(x); xi(i)=int(x(i)); end do
   par=estdlaplace2(xi,'ML',status=st)
   print '(a,2f10.5)', 'ML p,q: ',par
   print '(a,i0)', 'status: ',st
end program
