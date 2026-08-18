! SPDX-License-Identifier: GPL-3.0-or-later
program test_compositions_special
  use copula
  implicit none
  type(copula_model) :: models(2), base
  real(dp) :: weights(2), u(2), expected, total
  integer :: k
  u = [0.4_dp,0.8_dp]
  models(1) = independence_copula(2)
  models(2) = clayton_copula(2.0_dp)
  weights = [0.25_dp,0.75_dp]
  expected = 0.25_dp*product(u)+0.75_dp*pCopula(u,models(2))
  call assert_close(mixture_cdf(u,models,weights),expected,1.0e-14_dp)

  base = gumbel_copula(2.0_dp)
  call assert_close(khoudraji_cdf(u,base,[1.0_dp,1.0_dp]),pCopula(u,base),1.0e-14_dp)
  call assert_close(khoudraji_cdf(u,base,[0.0_dp,0.0_dp]),product(u),1.0e-14_dp)

  expected = (u(1)**(-3.0_dp)+u(2)**(-3.0_dp)-1.0_dp)**(-1.0_dp/3.0_dp)
  call assert_close(nested_clayton_cdf([u(1),u(2),1.0_dp],2,2.0_dp,3.0_dp),expected,1.0e-13_dp)

  call assert_true(stirling_first(5,2) == 50_i8)
  call assert_true(stirling_second(5,2) == 15_i8)
  call assert_true(eulerian_number(4,1) == 11_i8)

  total = 0.0_dp
  do k = 1, 1000
    total = total+sibuya_pmf(k,0.7_dp)
  end do
  call assert_close(total,1.0_dp,1.0e-2_dp)
  total = 0.0_dp
  do k = 1, 1000
    total = total+logseries_pmf(k,0.7_dp)
  end do
  call assert_close(total,1.0_dp,1.0e-10_dp)
  call assert_true(random_sibuya(0.7_dp,123_i8) >= 1)
  call assert_true(random_logseries(0.7_dp,321_i8) >= 1)

  print '(a)', 'test_compositions_special: PASS'
contains
  subroutine assert_close(actual, reference, tolerance)
    real(dp), intent(in) :: actual, reference, tolerance
    if (abs(actual-reference) > tolerance*(1.0_dp+abs(reference))) then
      print '(a,3es25.16)', 'mismatch: ', actual, reference, abs(actual-reference)
      error stop 1
    end if
  end subroutine assert_close
  subroutine assert_true(condition)
    logical, intent(in) :: condition
    if (.not. condition) error stop 1
  end subroutine assert_true
end program test_compositions_special
