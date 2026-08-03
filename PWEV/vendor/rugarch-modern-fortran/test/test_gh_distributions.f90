program test_gh_distributions
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use rugarch
   use rugarch_gh, only : bessel_k_log
   implicit none
   real(dp) :: x, p, q, expected
   integer :: i
   real(dp), parameter :: xs(5) = [0.01_dp,0.1_dp,1.0_dp,10.0_dp,100.0_dp]

   do i=1,size(xs)
      x=xs(i)
      expected=0.5_dp*log(acos(-1.0_dp)/(2.0_dp*x))-x
      if(abs(bessel_k_log(x,0.5_dp)-expected)>2.0e-9_dp) error stop 'Bessel K check failed'
   end do

   if(abs(distribution_pdf(0.1_dp,dist_ghyp,1.0_dp,0.2_dp)- &
      dsgh(0.1_dp,0.0_dp,1.0_dp,0.2_dp,1.0_dp,1.0_dp))>1.0e-11_dp) &
      error stop 'default generalized-hyperbolic lambda is not one'

   if(abs(dsnig(0.1_dp,0.0_dp,1.0_dp,0.2_dp,1.0_dp)- &
      dsgh(0.1_dp,0.0_dp,1.0_dp,0.2_dp,1.0_dp,-0.5_dp))>1.0e-11_dp) &
      error stop 'NIG is not the GH lambda=-1/2 special case'

   if(abs(dsghst(0.2_dp,0.0_dp,1.0_dp,0.0_dp,8.0_dp)- &
      dstd(0.2_dp,0.0_dp,1.0_dp,8.0_dp))>1.0e-11_dp) &
      error stop 'symmetric GH skew-t does not match standardized t'

   p=psnig(0.35_dp,0.0_dp,1.0_dp,0.25_dp,1.2_dp)
   q=qsnig(p,0.0_dp,1.0_dp,0.25_dp,1.2_dp)
   if(abs(q-0.35_dp)>2.0e-4_dp) error stop 'NIG CDF/quantile round trip failed'

   p=psgh(-0.25_dp,0.0_dp,1.0_dp,-0.15_dp,1.1_dp,0.7_dp)
   q=qsgh(p,0.0_dp,1.0_dp,-0.15_dp,1.1_dp,0.7_dp)
   if(abs(q+0.25_dp)>3.0e-4_dp) error stop 'GH CDF/quantile round trip failed'

   p=psghst(0.15_dp,0.0_dp,1.0_dp,0.5_dp,8.0_dp)
   q=qsghst(p,0.0_dp,1.0_dp,0.5_dp,8.0_dp)
   if(abs(q-0.15_dp)>3.0e-4_dp) error stop 'GH skew-t CDF/quantile round trip failed'

   call seed_rng(8347)
   do i=1,100
      if(.not.finite_value(rsnig(0.0_dp,1.0_dp,0.2_dp,1.0_dp))) error stop 'invalid NIG random value'
      if(.not.finite_value(rsgh(0.0_dp,1.0_dp,-0.1_dp,1.2_dp,0.5_dp))) error stop 'invalid GH random value'
      if(.not.finite_value(rsghst(0.0_dp,1.0_dp,0.4_dp,8.0_dp))) error stop 'invalid GHST random value'
   end do

   print '(a)', 'GH distribution tests passed'
contains
   pure elemental logical function finite_value(value)
      real(dp),intent(in)::value
      finite_value=ieee_is_finite(value)
   end function finite_value
end program test_gh_distributions
