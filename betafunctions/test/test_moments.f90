program test_moments
  use betafunctions
  implicit none
  type(moment_set) :: m
  type(beta_params) :: p
  real(dp) :: empty(0)

  call beta_moments(5.0_dp, 3.0_dp, 0.25_dp, 0.75_dp, 4, m)
  call check_close(m%raw(1), 0.5625_dp, 2.0e-14_dp, 'beta mean')
  call check_close(m%central(2), 0.006510416666666667_dp, 2.0e-14_dp, 'beta variance')

  call binomial_moments(10, 0.3_dp, 4, m)
  call check_close(m%raw(1), 3.0_dp, 2.0e-13_dp, 'binomial mean')
  call check_close(m%central(2), 2.1_dp, 2.0e-12_dp, 'binomial variance')

  p = beta2_fit(empty, 0.6_dp, 0.02_dp, 0.0_dp, 1.0_dp, .false.)
  call check_close(p%alpha, 6.6_dp, 2.0e-13_dp, 'beta2 alpha')
  call check_close(p%beta, 4.4_dp, 2.0e-13_dp, 'beta2 beta')

  call check_close(descending_factorial(7.0_dp, 3), 210.0_dp, 0.0_dp, 'falling factorial')
  call check_close(ascending_factorial(5.0_dp, 3), 210.0_dp, 0.0_dp, 'rising factorial')

  print '(a)', 'test_moments: PASS'

contains
  subroutine check_close(got, want, tol, label)
    real(dp), intent(in) :: got, want, tol
    character(*), intent(in) :: label
    if (abs(got - want) > tol) then
      write(*,'(a,2es24.15)') trim(label)//' mismatch: ', got, want
      error stop 1
    end if
  end subroutine check_close
end program test_moments
