program test_distribution
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   use argus, only : dp, dargus, pargus, qargus, &
      dargus_recycle, pargus_recycle, qargus_recycle
   implicit none

   real(dp), parameter :: tol = 2.0e-12_dp
   real(dp), parameter :: chis(7) = [0.001_dp, 0.01_dp, 0.1_dp, 0.6_dp, 1.0_dp, 2.0_dp, 6.0_dp]
   real(dp), parameter :: probs(5) = [0.01_dp, 0.1_dp, 0.5_dp, 0.9_dp, 0.99_dp]
   real(dp) :: x, p, q, lim, ans(4)
   real(dp) :: xa(2), ca(4), pa(2)
   integer :: i, j

   call assert_close(dargus(0.3_dp,1.0_dp), 0.7289121037966952_dp, 5.0e-15_dp, "density reference")
   call assert_close(pargus(0.3_dp,1.0_dp), 0.1094952963875544_dp, 5.0e-15_dp, "cdf reference")
   call assert_close(qargus(0.9_dp,2.0_dp), 0.9398552320663749_dp, 2.0e-14_dp, "quantile reference")
   call assert_close(pargus(0.9_dp,6.0_dp,lower=.false.), 0.9228238746473566_dp, 5.0e-15_dp, "upper-tail reference")

   call assert_close(exp(dargus(0.3_dp,1.0_dp,log_pdf=.true.)), dargus(0.3_dp,1.0_dp), tol, "log density")
   call assert_close(exp(pargus(0.9_dp,6.0_dp,lower=.false.,log_p=.true.)), &
      pargus(0.9_dp,6.0_dp,lower=.false.), tol, "log upper cdf")
   call assert_close(qargus(log(0.1_dp),2.0_dp,lower=.false.,log_p=.true.), &
      qargus(0.9_dp,2.0_dp), 2.0e-14_dp, "log upper quantile")

   do i = 1, size(chis)
      do j = 1, size(probs)
         p = probs(j)
         q = qargus(p,chis(i))
         call assert_close(pargus(q,chis(i)), p, 2.0e-11_dp, "cdf-quantile round trip")
         q = qargus(1.0_dp-p,chis(i),lower=.false.)
         call assert_close(pargus(q,chis(i)), p, 2.0e-11_dp, "upper quantile round trip")
      end do
   end do

   ! As chi -> 0, f(x) -> 3*x*sqrt(1-x^2) and
   ! F(x) -> 1-(1-x^2)^(3/2).
   x = 0.3_dp
   lim = 3.0_dp*x*sqrt(1.0_dp-x*x)
   call assert_close(dargus(x,1.0e-8_dp), lim, 2.0e-14_dp, "small-chi density limit")
   lim = 1.0_dp-(1.0_dp-x*x)**1.5_dp
   call assert_close(pargus(x,1.0e-8_dp), lim, 2.0e-14_dp, "small-chi cdf limit")

   call assert_close(dargus(-0.1_dp,2.0_dp), 0.0_dp, 0.0_dp, "density below support")
   call assert_close(dargus(1.1_dp,2.0_dp), 0.0_dp, 0.0_dp, "density above support")
   call assert_close(pargus(-0.1_dp,2.0_dp), 0.0_dp, 0.0_dp, "cdf below support")
   call assert_close(pargus(1.1_dp,2.0_dp), 1.0_dp, 0.0_dp, "cdf above support")
   call assert_close(qargus(0.0_dp,2.0_dp), 0.0_dp, 0.0_dp, "quantile zero")
   call assert_close(qargus(1.0_dp,2.0_dp), 1.0_dp, 0.0_dp, "quantile one")
   if (.not. ieee_is_nan(dargus(0.5_dp,-1.0_dp))) error stop "invalid chi should return NaN"

   ! Elemental scalar expansion.
   ca = [0.3_dp, 0.6_dp, 1.0_dp, 2.0_dp]
   ans = dargus(0.3_dp,ca)
   do i = 1, 4
      call assert_close(ans(i),dargus(0.3_dp,ca(i)),0.0_dp,"elemental expansion")
   end do

   ! Explicit R-style recycling helpers.
   xa = [0.2_dp,0.7_dp]
   call dargus_recycle(xa,ca,ans)
   do i = 1, 4
      call assert_close(ans(i),dargus(xa(1+mod(i-1,2)),ca(i)),0.0_dp,"density recycle")
   end do
   call pargus_recycle(xa,ca,ans)
   do i = 1, 4
      call assert_close(ans(i),pargus(xa(1+mod(i-1,2)),ca(i)),0.0_dp,"cdf recycle")
   end do
   pa = [0.2_dp,0.8_dp]
   call qargus_recycle(pa,ca,ans)
   do i = 1, 4
      call assert_close(ans(i),qargus(pa(1+mod(i-1,2)),ca(i)),0.0_dp,"quantile recycle")
   end do

   print '(a)', "test_distribution: PASS"

contains

   subroutine assert_close(actual, expected, atol, label)
      real(dp), intent(in) :: actual, expected, atol
      character(*), intent(in) :: label
      real(dp) :: err, scale

      err = abs(actual-expected)
      scale = max(1.0_dp,abs(expected))
      if (err > max(atol,20.0_dp*epsilon(1.0_dp)*scale)) then
         print '(a)', "FAIL: "//label
         print '(a,es24.16)', " actual   = ", actual
         print '(a,es24.16)', " expected = ", expected
         print '(a,es24.16)', " error    = ", err
         error stop 1
      end if
   end subroutine assert_close

end program test_distribution
