program test_parity
  use mnormt
  implicit none
  real(dp) :: mu4(4), s4(4,4), x4(4), p
  real(dp) :: mu2(2), s2(2,2), lo2(2), up2(2), p1, p2
  real(dp) :: mu(2), s(2,2), lo(2), up(2)
  real(dp), allocatable :: raw(:)
  integer :: i, fails, kappa(2)
  type(probability_result) :: pr
  type(trunc_moment_result) :: mr, cr
  type(mardia_result) :: md
  real(dp) :: data(8,2)

  fails = 0
  mu4 = 0.0_dp
  s4 = 0.0_dp
  do i = 1, 4
    s4(i,i) = 1.0_dp
  end do
  x4 = 0.0_dp
  p = pmnorm(x4, mu4, s4, 40000, 1.0e-9_dp, 0.0_dp)
  if (abs(p - 0.0625_dp) > 2.0e-6_dp) then
    print *, '4d mvn fail', p
    fails = fails + 1
  end if

  mu2 = 0.0_dp
  s2 = reshape([1.0_dp, 0.3_dp, 0.3_dp, 1.0_dp], [2,2])
  lo2 = [-1.0_dp, -2.0_dp]
  up2 = [0.5_dp, 1.5_dp]
  p1 = biv_nt_prob(5.0_dp, lo2, up2, mu2, s2)
  pr = sadmvt_prob(5.0_dp, lo2, up2, mu2, s2, 40000, 1.0e-9_dp, 0.0_dp)
  p2 = pr%value
  if (abs(p1 - p2) > 2.0e-6_dp) then
    print *, 'biv/sadmvt fail', p1, p2
    fails = fails + 1
  end if

  mu = [0.5_dp, -1.0_dp]
  s = reshape([3.0_dp, 3.0_dp, 3.0_dp, 6.0_dp], [2,2])
  lo = [-1.0_dp, -2.8_dp]
  up = [1.5_dp, 1.5_dp]
  kappa = [2,2]
  mr = mom_mtruncnorm(kappa, mu, s, lo, up)
  if (abs(mr%probability - 0.3826730991937587_dp) > 4.0e-8_dp) then
    print *, 'corr moment p fail', mr%probability
    fails = fails + 1
  end if
  if (maxval(abs(mr%mean - &
      [0.3438244816828486_dp, -0.8530756310097924_dp])) > 3.0e-8_dp) then
    print *, 'corr moment mean fail', mr%mean
    fails = fails + 1
  end if
  if (maxval(abs(mr%covariance - reshape([ &
      0.46067645_dp, 0.18450576_dp, 0.18450576_dp, 1.26780832_dp], &
      [2,2]))) > 5.0e-8_dp) then
    print *, 'corr moment cov fail', mr%covariance
    fails = fails + 1
  end if

  allocate(raw(5))
  raw = [1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 3.0_dp]
  cr = mom2cum(raw, [5])
  if (abs(cr%covariance(1,1) - 1.0_dp) > 1.0e-14_dp .or. &
      abs(cr%marginal_skewness(1)) > 1.0e-14_dp .or. &
      abs(cr%marginal_excess_kurtosis(1)) > 1.0e-14_dp) then
    print *, 'mom2cum fail'
    fails = fails + 1
  end if

  data = reshape([ &
      -2.0_dp, -1.0_dp, 0.0_dp, 1.0_dp, 2.0_dp, -1.5_dp, 0.5_dp, 1.5_dp, &
       0.0_dp,  1.0_dp,-1.0_dp, 2.0_dp,-2.0_dp,  0.5_dp,-0.5_dp, 1.5_dp], &
      [8,2])
  md = sample_mardia_measures(data)
  if (md%status /= 0 .or. md%b1 < 0.0_dp .or. md%b2 < 0.0_dp) then
    print *, 'mardia fail'
    fails = fails + 1
  end if

  if (fails == 0) then
    print *, 'test_parity: PASS'
  else
    print *, 'test_parity: FAIL', fails
    error stop 1
  end if
end program test_parity
