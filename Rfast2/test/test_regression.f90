program test_regression
   use rfast2
   implicit none
   type(regression_result) :: lr,pr,zr
   type(multinomial_result) :: mr
   real(dp) :: x(10,1),yb(10),yp(10)
   integer :: ym(10),yz(10),i

   do i=1,10
      x(i,1)=real(i-5,dp)/2.0_dp
   end do
   yb=[0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp]
   lr=logistic_reg(yb,x)
   if (lr%status /= 0 .or. size(lr%beta) /= 2 .or. lr%beta(2) <= 0.0_dp) error stop 1
   yp=[1.0_dp,1.0_dp,1.0_dp,2.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp,6.0_dp,8.0_dp]
   pr=poisson_reg(yp,x)
   if (pr%status /= 0 .or. pr%beta(2) <= 0.0_dp) error stop 2
   ym=[1,2,3,1,2,3,1,2,3,1]
   mr=multinom_reg(ym,x)
   if (mr%status /= 0) error stop 3
   if (.not. allocated(mr%beta)) error stop 31
   if (size(mr%beta,1) /= 2 .or. size(mr%beta,2) /= 2) error stop 32
   yz=[1,1,1,2,2,3,4,5,7,9]
   zr=ztp_reg(yz,x)
   if (zr%status /= 0 .or. any(zr%fitted <= 1.0_dp)) error stop 4
   print '(a)', 'test_regression: PASS'
end program test_regression
