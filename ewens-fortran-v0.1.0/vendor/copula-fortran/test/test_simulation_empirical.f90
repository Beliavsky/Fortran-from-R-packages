! SPDX-License-Identifier: GPL-3.0-or-later
program test_simulation_empirical
  use copula
  implicit none
  type(copula_model) :: model
  type(test_result) :: test
  real(dp), allocatable :: sample(:,:), ranks(:,:), transformed(:,:), recovered(:,:)
  real(dp) :: x(5,2), point(2), tau_sample
  logical :: ok

  model = clayton_copula(2.0_dp)
  call rCopula(5000,model,sample,ok,998877_i8)
  call assert_true(ok)
  tau_sample = sample_kendall_tau(sample(:,1),sample(:,2))
  call assert_close(tau_sample,0.5_dp,0.035_dp)

  call rosenblatt_transform(sample(1:500,:),model,transformed,ok)
  call assert_true(ok)
  call inverse_rosenblatt(transformed,model,recovered,ok)
  call assert_true(ok)
  call assert_close(maxval(abs(recovered-sample(1:500,:))),0.0_dp,2.0e-5_dp)
  call assert_close(sum(transformed(:,2))/real(size(transformed,1),dp),0.5_dp,0.04_dp)

  x = reshape([1.0_dp,2.0_dp,2.0_dp,4.0_dp,5.0_dp, &
               5.0_dp,4.0_dp,3.0_dp,2.0_dp,1.0_dp],[5,2])
  call pobs(x,ranks,ok)
  call assert_true(ok)
  call assert_true(all(ranks > 0.0_dp) .and. all(ranks < 1.0_dp))
  call assert_close(ranks(2,1),ranks(3,1),1.0e-14_dp)
  point = [0.6_dp,0.6_dp]
  call assert_true(empirical_copula(ranks,point) >= 0.0_dp)

  call rCopula(150,independence_copula(2),sample,ok,1234_i8)
  test = independence_test(sample,99,4321_i8)
  call assert_true(test%ok)
  call assert_true(test%p_value > 0.01_dp)
  test = exchangeability_test(sample,49,9876_i8)
  call assert_true(test%ok)
  call assert_true(test%p_value > 0.01_dp)
  test = radial_symmetry_test(sample,49,5432_i8)
  call assert_true(test%ok)
  call assert_true(test%p_value > 0.01_dp)

  print '(a)', 'test_simulation_empirical: PASS'
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
end program test_simulation_empirical
