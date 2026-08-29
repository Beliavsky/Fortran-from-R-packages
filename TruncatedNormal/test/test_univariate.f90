program test_univariate
 use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_positive_inf, ieee_negative_inf
 use truncated_normal
 use r_compat, only: set_seed_int
 implicit none
 real(dp) :: mu1(1), sig1(1,1), lo1(1), hi1(1)
 real(dp) :: mu2(2), sig2(2,2), x(2,2), lo2(2), hi2(2), pinf, ninf
 real(dp) :: a(2), b(2), lp(2), d0(2), d1(2), t0(2), t1(2)
 type(prob_result) :: pr

 call set_seed_int(101)
 mu1 = 0.0_dp
 sig1 = reshape([1.0_dp], [1,1])
 lo1 = -1.0_dp
 hi1 = 1.0_dp
 pr = pmvnorm(mu1, sig1, lo1, hi1)
 call check_close('univariate normal probability', pr%prob, 0.6826894921370859_dp, 2.0e-13_dp)
 pr = pmvt(mu1, sig1, 5.0_dp, lo1, hi1)
 call check_close('univariate Student probability', pr%prob, 0.6367825323508771_dp, 2.0e-10_dp)

 a = [8.0_dp, -9.0_dp]
 b = [9.0_dp, -8.0_dp]
 lp = lnNpr(a, b)
 if (any(.not. (lp < 0.0_dp))) error stop 'lnNpr extreme-tail test failed'

 pinf = ieee_value(1.0_dp, ieee_positive_inf)
 ninf = ieee_value(1.0_dp, ieee_negative_inf)
 mu2 = [0.2_dp, -0.1_dp]
 sig2 = reshape([1.0_dp, 0.4_dp, 0.4_dp, 1.0_dp], [2,2])
 x(1,:) = [0.5_dp, 0.2_dp]
 x(2,:) = [-0.1_dp, 0.6_dp]
 lo2 = ninf
 hi2 = pinf
 d0 = dmvnorm(x, mu2, sig2)
 d1 = dtmvnorm(x, mu2, sig2, lo2, hi2)
 call check_vec_close('untruncated dtmvnorm', d1, d0, 2.0e-14_dp)
 call check_close('mvn density reference 1', d0(1), 0.16284017_dp, 1.0e-8_dp)
 call check_close('mvn density reference 2', d0(2), 0.11125410_dp, 1.0e-8_dp)

 t0 = dmvt(x, mu2, sig2, 5.0_dp)
 t1 = dtmvt(x, mu2, sig2, 5.0_dp, lo2, hi2)
 call check_vec_close('untruncated dtmvt', t1, t0, 2.0e-14_dp)
 t1 = dtmvt(x, mu2, sig2, 0.0_dp, lo2, hi2)
 call check_vec_close('df=0 normal density limit', t1, d0, 2.0e-14_dp)
 t1 = dtmvt(x, mu2, sig2, pinf, lo2, hi2)
 call check_vec_close('df=Inf normal density limit', t1, d0, 2.0e-14_dp)
 call check_close('mvt density reference 1', t0(1), 0.15888689_dp, 1.0e-8_dp)
 call check_close('mvt density reference 2', t0(2), 0.09784735_dp, 1.0e-8_dp)

 print '(a)', 'test_univariate: PASS'
contains
 subroutine check_close(name, got, want, tol)
  character(len=*), intent(in) :: name
  real(dp), intent(in) :: got, want, tol
  if (abs(got-want) > tol) then
   print '(a,2(1x,es24.16))', trim(name)//' failed:', got, want
   error stop 1
  end if
 end subroutine check_close
 subroutine check_vec_close(name, got, want, tol)
  character(len=*), intent(in) :: name
  real(dp), intent(in) :: got(:), want(:), tol
  if (maxval(abs(got-want)) > tol) then
   print '(a)', trim(name)//' failed'
   error stop 1
  end if
 end subroutine check_vec_close
end program test_univariate
