program midas_weights
  use rumidas
  implicit none
  real(dp), allocatable :: beta(:), almon(:)
  integer :: status, i

  beta = beta_weights(12, 1.0_dp, 5.0_dp, status)
  almon = exponential_almon_weights(12, 0.1_dp, -0.05_dp, status)
  print '(a)', ' lag       beta       almon'
  do i = 1, 12
    print '(i4,2f12.6)', i, beta(i), almon(i)
  end do
end program midas_weights
