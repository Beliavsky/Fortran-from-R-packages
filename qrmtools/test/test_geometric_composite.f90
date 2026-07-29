! SPDX-License-Identifier: GPL-3.0-or-later
program test_geometric_composite
  use qrmtools, only : dp, distribution_component, distribution_normal, &
    composite_pdf, composite_cdf, composite_quantile, geometric_var, &
    geometric_expectile, fit_result, geometric_objective
  implicit none

  type(distribution_component) :: component(2)
  type(fit_result) :: fit
  real(dp) :: cuts(1)
  real(dp) :: weights(2)
  real(dp) :: data(5,2)
  real(dp) :: level(2)
  real(dp) :: quantile

  cuts = [0.0_dp]
  weights = [0.4_dp,0.6_dp]
  component(1)%family = distribution_normal
  component(1)%parameters(1:2) = [-1.0_dp,1.0_dp]
  component(2)%family = distribution_normal
  component(2)%parameters(1:2) = [1.0_dp,1.0_dp]

  call assert_close(composite_cdf(0.0_dp,cuts,component,weights),0.4_dp,2.0e-12_dp)
  call assert_true(composite_pdf(-0.5_dp,cuts,component,weights)>0.0_dp)
  quantile = composite_quantile(0.75_dp,cuts,component,weights)
  call assert_close(composite_cdf(quantile,cuts,component,weights),0.75_dp,2.0e-10_dp)

  data = reshape([-2.0_dp,-1.0_dp,0.0_dp,1.0_dp,2.0_dp, &
                  -1.0_dp,-0.5_dp,0.0_dp,0.5_dp,1.0_dp],[5,2])
  level = [0.2_dp,0.1_dp]
  fit = geometric_var(data,level,max_iterations=3000)
  call assert_true(fit%ok)
  call assert_true(size(fit%parameters)==2)
  call assert_true(geometric_objective(fit%parameters,data,level) <= &
    geometric_objective([0.0_dp,0.0_dp],data,level)+1.0e-8_dp)

  fit = geometric_expectile(data,level,max_iterations=3000)
  call assert_true(fit%ok)
  call assert_true(size(fit%parameters)==2)

  print '(a)', 'test_geometric_composite: PASS'

contains

  subroutine assert_close(actual,expected,tolerance)
    real(dp), intent(in) :: actual,expected,tolerance
    if(abs(actual-expected)>tolerance*max(1.0_dp,abs(expected))) then
      print *, 'mismatch:',actual,expected,abs(actual-expected)
      error stop 1
    end if
  end subroutine assert_close

  subroutine assert_true(condition)
    logical, intent(in) :: condition
    if(.not.condition) error stop 1
  end subroutine assert_true

end program test_geometric_composite
