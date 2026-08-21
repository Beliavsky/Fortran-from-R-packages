program test_distributions
  use betafunctions
  implicit none
  real(dp) :: v, q, s
  integer :: k

  call check_close(beta4_pdf(0.4_dp, 0.25_dp, 0.75_dp, 5.0_dp, 3.0_dp), &
                   0.83349_dp, 2.0e-12_dp, 'beta4 pdf')
  call check_close(beta4_cdf(0.4_dp, 0.25_dp, 0.75_dp, 5.0_dp, 3.0_dp), &
                   0.0287955_dp, 2.0e-12_dp, 'beta4 cdf')
  q = beta4_quantile(0.7_dp, 0.25_dp, 0.75_dp, 5.0_dp, 3.0_dp)
  call check_close(q, 0.6118301530475466_dp, 2.0e-11_dp, 'beta4 quantile')

  call check_close(beta_binomial_pmf(5, 20, 0.0_dp, 1.0_dp, 2.5_dp, 4.0_dp), &
                   0.0851736132455381_dp, 2.0e-10_dp, 'beta-binomial pmf')
  call check_close(beta_binomial_cdf(5, 20, 0.0_dp, 1.0_dp, 2.5_dp, 4.0_dp), &
                   0.2461646473725213_dp, 3.0e-9_dp, 'beta-binomial P(X<q) convention')
  s = 0.0_dp
  do k = 0, 20
    s = s + beta_binomial_pmf(k, 20, 0.0_dp, 1.0_dp, 2.5_dp, 4.0_dp)
  end do
  call check_close(s, 1.0_dp, 2.0e-9_dp, 'beta-binomial normalization')

  call check_close(gamma_binomial_pdf(2.5_dp, 10.5_dp, 0.3_dp, .true.), &
                   0.2543061124405252_dp, 3.0e-9_dp, 'gamma-binomial normalized pdf')
  call check_close(gamma_binomial_cdf(5.0_dp, 10.5_dp, 0.3_dp), &
                   0.8861273115023546_dp, 3.0e-9_dp, 'gamma-binomial cdf')
  call check_close(beta_times_beta_tail(0.4_dp, 0.1_dp, 0.9_dp, 5.0_dp, 3.0_dp, &
                   20.0_dp, 0.5_dp, .true.), 0.8317407537106192_dp, 3.0e-10_dp, 'beta times beta lower tail')
  call check_close(beta_times_beta_tail(0.4_dp, 0.1_dp, 0.9_dp, 5.0_dp, 3.0_dp, &
                   20.0_dp, 0.5_dp, .false.), 0.18213425391877536_dp, 3.0e-10_dp, 'beta times beta upper tail')

  s = 0.0_dp
  do k = 0, 20
    s = s + compound_binomial_pmf(k, 20, 1.2_dp, 0.4_dp)
  end do
  call check_close(s, 1.0_dp, 2.0e-12_dp, 'compound-binomial normalization')

  v = beta_ms_cdf(0.4_dp, 0.5625_dp, 0.006510416666666667_dp, 0.25_dp, 0.75_dp)
  call check_close(v, 0.0287955_dp, 2.0e-11_dp, 'mean/variance beta cdf')

  print '(a)', 'test_distributions: PASS'

contains

  subroutine check_close(got, want, tol, label)
    real(dp), intent(in) :: got, want, tol
    character(*), intent(in) :: label
    if (abs(got - want) > tol) then
      write(*,'(a,2es24.15)') trim(label)//' mismatch: ', got, want
      error stop 1
    end if
  end subroutine check_close

end program test_distributions
