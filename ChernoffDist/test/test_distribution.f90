program test_distribution
  use chernoffdist, only: dp, dchern, pchern
  implicit none

  real(dp), parameter :: x(5) = [0.0_dp, 0.5_dp, 1.0_dp, 1.5_dp, 2.0_dp]
  real(dp), parameter :: dref(5) = [ &
    0.7583445580539164_dp, 0.4912075008596545_dp, 0.1208801617403622_dp, &
    0.008618552068725742_dp, 0.0001212017887020576_dp ]
  real(dp), parameter :: pref(5) = [ &
    0.5_dp, 0.8311165703196535_dp, 0.9752206565664467_dp, &
    0.9988642217326905_dp, 0.9999891453806221_dp ]
  real(dp) :: d, p
  integer :: i, failures

  failures = 0
  do i = 1, size(x)
    d = dchern(x(i))
    p = pchern(x(i))
    if (abs(d - dref(i)) > 2.0e-9_dp * max(1.0_dp, abs(dref(i)))) failures = failures + 1
    if (abs(p - pref(i)) > 3.0e-9_dp) failures = failures + 1
  end do

  if (failures /= 0) then
    print '(a,i0)', 'test_distribution: FAIL ', failures
    error stop 1
  end if
  print '(a)', 'test_distribution: PASS'
end program test_distribution
