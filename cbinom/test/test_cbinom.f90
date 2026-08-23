program test_cbinom
   use cbinom, only : dp, dcbinom, pcbinom, qcbinom, rcbinom
   implicit none
   integer :: failures, k, n
   real(dp) :: p, f1, f2, q, x, s, d, meanx
   real(dp), allocatable :: r(:)

   failures = 0
   n = 20
   p = 0.2_dp

   ! For integer n, F_cbinom(k+1) equals ordinary Binomial(n,p) CDF at k.
   do k = 0, n
      f1 = pcbinom(real(k+1,dp), real(n,dp), p)
      f2 = binom_cdf(k, n, p)
      if (abs(f1-f2) > 2.0e-10_dp) then
         print *, 'integer-CDF parity failed', k, f1, f2
         failures = failures + 1
      end if
   end do

   ! CDF/quantile inversion for a non-integer size.
   do k = 1, 19
      q = real(k,dp)/20.0_dp
      x = qcbinom(q, 7.5_dp, 0.37_dp)
      f1 = pcbinom(x, 7.5_dp, 0.37_dp)
      if (abs(f1-q) > 2.0e-7_dp) then
         print *, 'quantile inversion failed', q, x, f1
         failures = failures + 1
      end if
   end do

   ! Density is nonnegative and integrates approximately to 1.
   s = 0.0_dp
   do k = 0, 2000
      x = 9.25_dp*real(k,dp)/2000.0_dp
      d = dcbinom(x, 8.25_dp, 0.61_dp)
      if (d < -1.0e-12_dp) failures = failures + 1
      if (k == 0 .or. k == 2000) then
         s = s + 0.5_dp*d
      else
         s = s + d
      end if
   end do
   s = s*9.25_dp/2000.0_dp
   if (abs(s-1.0_dp) > 3.0e-3_dp) then
      print *, 'density integral failed', s
      failures = failures + 1
   end if

   allocate(r(20000))
   call rcbinom(r, 10.0_dp, 0.3_dp)
   if (any(r < 0.0_dp) .or. any(r > 11.0_dp)) then
      print *, 'random support failed'
      failures = failures + 1
   end if
   meanx = sum(r)/real(size(r),dp)
   if (abs(meanx - 3.5_dp) > 0.08_dp) then
      print *, 'random mean failed', meanx
      failures = failures + 1
   end if

   if (failures == 0) then
      print *, 'test_cbinom: PASS'
   else
      print *, 'test_cbinom: FAIL', failures
      error stop 1
   end if

contains

   pure function binom_cdf(k, n, p) result(ans)
      integer, intent(in) :: k, n
      real(dp), intent(in) :: p
      real(dp) :: ans, term
      integer :: j
      if (k < 0) then
         ans = 0.0_dp
         return
      else if (k >= n) then
         ans = 1.0_dp
         return
      end if
      if (p == 0.0_dp) then
         ans = 1.0_dp
         return
      else if (p == 1.0_dp) then
         ans = 0.0_dp
         return
      end if
      term = (1.0_dp-p)**n
      ans = term
      do j = 0, k-1
         term = term * real(n-j,dp)/real(j+1,dp) * p/(1.0_dp-p)
         ans = ans + term
      end do
   end function binom_cdf

end program test_cbinom
