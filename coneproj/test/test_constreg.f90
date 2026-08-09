program test_constreg
   use coneproj
   implicit none
   integer, parameter :: n=8
   real(dp) :: x(n,2), y(n), a(1,2)
   type(constreg_result) :: ans
   integer :: i
   do i=1,n
      x(i,1)=1.0_dp
      x(i,2)=real(i-1,dp)/real(n-1,dp)
   end do
   y = 1.5_dp + 2.0_dp*x(:,2)
   a = reshape([0.0_dp,1.0_dp],[1,2])
   call seed_rng(1234)
   call constreg_fit(y,x,a,ans,nsim_cov=20)
   if (ans%status /= coneproj_success) error stop 'constreg status'
   if (maxval(abs(ans%coefs-[1.5_dp,2.0_dp])) > 1.0e-7_dp) then
      print *, ans%coefs
      error stop 'constreg coefficients'
   end if
   print *, 'test_constreg: PASS'
end program test_constreg
