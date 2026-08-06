program influence_demo
  use rpeif, only : dp, rpeif_options, influence_result, influence_series
  implicit none
  real(dp) :: returns(16)
  type(rpeif_options) :: options
  type(influence_result) :: result
  integer :: i

  returns = [-0.030_dp, 0.012_dp, 0.008_dp, -0.011_dp, 0.016_dp, 0.005_dp, &
    -0.018_dp, 0.021_dp, 0.009_dp, -0.006_dp, 0.013_dp, 0.004_dp, &
    -0.090_dp, 0.018_dp, 0.011_dp, 0.007_dp]

  options%alpha = 0.25_dp
  options%clean_outliers = .false.
  options%efficiency = 0.99_dp
  options%source_compatibility = .false.

  call influence_series('ES', returns, result, options)
  print '(a)', 'Expected-shortfall influence series'
  print '(a)', ' index       cleaned return        influence'
  do i = 1, size(result%values)
    print '(i6,2f22.10)', i, result%x(i), result%values(i)
  end do
  print '(a,i0)', 'status = ', result%status
end program influence_demo
