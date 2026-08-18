! SPDX-License-Identifier: MIT
program test_distribution
  use ewens, only : dp, dewens, dewens_counts, dewens_k, ewens_k_exact
  implicit none
  integer :: labels(6), one_class(4), mixed(4), mvec(3), k
  real(dp) :: p, total, mean_k, target

  labels = [1, 1, 1, 2, 2, 3]
  mvec = [1, 1, 1]
  call assert_close(dewens(labels, 2.0_dp), 4.0_dp / 21.0_dp, 2.0e-14_dp)
  call assert_close(dewens_counts(mvec, 2.0_dp), 4.0_dp / 21.0_dp, 2.0e-14_dp)

  one_class = 7
  mixed = [1, 1, 2, 2]
  call assert_close(dewens(one_class, 0.0_dp), 1.0_dp, 0.0_dp)
  call assert_close(dewens(mixed, 0.0_dp), 0.0_dp, 0.0_dp)

  call assert_close(dewens_k(1, 20, 1.0_dp), 0.05_dp, 2.0e-15_dp)
  call assert_close(ewens_k_exact(20, 1.0_dp), 3.597739657143682_dp, 2.0e-15_dp)

  total = 0.0_dp
  mean_k = 0.0_dp
  do k = 1, 50
    p = dewens_k(k, 50, 2.3_dp)
    total = total + p
    mean_k = mean_k + real(k, dp) * p
  end do
  target = ewens_k_exact(50, 2.3_dp)
  call assert_close(total, 1.0_dp, 3.0e-13_dp)
  call assert_close(mean_k, target, 2.0e-12_dp)

  print '(a)', 'test_distribution: PASS'
contains
  subroutine assert_close(x, y, tol)
    real(dp), intent(in) :: x, y, tol
    if (abs(x - y) > tol) then
      print '(a,3es24.15)', 'assert_close failed: ', x, y, abs(x - y)
      error stop 1
    end if
  end subroutine assert_close
end program test_distribution
