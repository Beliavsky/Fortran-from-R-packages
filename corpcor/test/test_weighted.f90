program test_weighted
  use corpcor, only : dp, moments_result, scale_result, wt_var, wt_moments, wt_scale
  implicit none
  real(dp), parameter :: tol = 1.0e-12_dp
  real(dp) :: x(6, 4), w(6), v
  type(moments_result) :: mom
  type(scale_result) :: sc

  x = reshape([ &
    1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp, &
    2.0_dp, 1.0_dp, 5.0_dp, 4.0_dp, 7.0_dp, 8.0_dp, &
    3.0_dp, 4.0_dp, 2.0_dp, 5.0_dp, 6.0_dp, 7.0_dp, &
    4.0_dp, 3.0_dp, 6.0_dp, 7.0_dp, 8.0_dp, 9.0_dp], [6, 4])
  w = [1.0_dp, 2.0_dp, 1.0_dp, 3.0_dp, 2.0_dp, 1.0_dp]
  w = w / sum(w)

  mom = wt_moments(x, w)
  call assert_close(maxval(abs(mom%mean - [3.6_dp, 4.3_dp, 4.7_dp, 6.2_dp])), 0.0_dp, tol)
  call assert_close(maxval(abs(mom%variance - [2.8_dp, 7.0125_dp, 2.5125_dp, 5.2_dp])), 0.0_dp, tol)
  v = wt_var(x(:, 1), w)
  call assert_close(v, 2.8_dp, tol)

  sc = wt_scale(x, w)
  call assert_close(maxval(abs(sum(spread(w, 2, 4) * sc%x, dim=1))), 0.0_dp, tol)
  call assert_close(maxval(abs(sum(spread(w, 2, 4) * sc%x**2, dim=1) / &
    (1.0_dp - sum(w*w)) - 1.0_dp)), 0.0_dp, tol)
  print '(a)', 'test_weighted: PASS'
contains
  subroutine assert_close(a, b, tolerance)
    real(dp), intent(in) :: a, b, tolerance
    if (abs(a - b) > tolerance) error stop 'assert_close failed'
  end subroutine assert_close
end program test_weighted
