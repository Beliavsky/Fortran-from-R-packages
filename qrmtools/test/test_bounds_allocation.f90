! SPDX-License-Identifier: GPL-3.0-or-later
module test_bounds_callbacks
  use qrmtools, only : dp, qpar
  implicit none
contains
  real(dp) function pareto_quantile(p) result(value)
    real(dp), intent(in) :: p
    value = qpar(p, 2.5_dp)
  end function pareto_quantile

  real(dp) function marginal_quantile(p, margin) result(value)
    real(dp), intent(in) :: p
    integer, intent(in) :: margin
    value = qpar(p, 2.5_dp, real(margin,dp))
  end function marginal_quantile
end module test_bounds_callbacks

program test_bounds_allocation
  use qrmtools, only : dp, hierarchy_node, hierarchical_matrix, alloc_ellip, &
    alloc_np, allocation_result, rearrange_matrix, rearrangement_result, &
    bound_worst_var, crude_var_bounds_hom, pareto_var_bounds_hom, &
    qpar, normal_cdf, ra_bounds, adaptive_ra_bounds, ra_bounds_result
  use test_bounds_callbacks, only : pareto_quantile, marginal_quantile
  implicit none

  type(hierarchy_node) :: root
  type(allocation_result) :: allocation
  type(rearrangement_result) :: rearranged
  type(ra_bounds_result) :: ra
  real(dp), allocatable :: matrix(:,:)
  real(dp), allocatable :: weights(:)
  real(dp) :: x(5,2)
  real(dp) :: input(4,3)
  real(dp) :: bounds(2)
  integer :: j

  allocate(root%components(1),root%children(1))
  root%value = 0.1_dp
  root%components = [1]
  root%children(1)%value = 0.6_dp
  root%children(1)%components = [2,3]
  matrix = hierarchical_matrix(root,[1.0_dp,2.0_dp,3.0_dp])
  call assert_close(matrix(1,2),0.1_dp,1.0e-14_dp)
  call assert_close(matrix(2,3),0.6_dp,1.0e-14_dp)
  call assert_close(matrix(3,3),3.0_dp,1.0e-14_dp)

  weights = alloc_ellip(4.0_dp,[1.0_dp,2.0_dp], &
    reshape([2.0_dp,0.5_dp,0.5_dp,1.0_dp],[2,2]))
  call assert_close(sum(weights-[1.0_dp,2.0_dp]),4.0_dp,1.0e-13_dp)

  x = reshape([1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp, &
               1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp],[5,2])
  allocation = alloc_np(x,0.6_dp,1.0_dp,.true.)
  call assert_true(allocation%ok)
  call assert_true(allocation%n==2)
  call assert_array_close(allocation%allocation,[4.5_dp,1.0_dp],1.0e-13_dp)

  input = reshape([1.0_dp,2.0_dp,3.0_dp,4.0_dp, &
                   1.0_dp,2.0_dp,3.0_dp,4.0_dp, &
                   1.0_dp,2.0_dp,3.0_dp,4.0_dp],[4,3])
  rearranged = rearrange_matrix(input,bound_worst_var,tolerance=0.0_dp, &
    n_lookback=3,max_iterations=200,already_sorted=.true.)
  call assert_true(rearranged%ok)
  do j=1,3
    call assert_close(sum(rearranged%rearranged(:,j)),10.0_dp,1.0e-13_dp)
    call assert_close(sum(rearranged%rearranged(:,j)**2),30.0_dp,1.0e-13_dp)
  end do
  call assert_true(rearranged%bound>=6.0_dp-1.0e-12_dp)

  bounds = crude_var_bounds_hom(0.95_dp,3,pareto_quantile)
  call assert_close(bounds(1),3.0_dp*qpar(0.95_dp/3.0_dp,2.5_dp),1.0e-12_dp)
  call assert_true(bounds(2)>bounds(1))
  bounds = pareto_var_bounds_hom(0.95_dp,2,2.5_dp)
  call assert_true(bounds(2)>bounds(1))

  ra = ra_bounds(0.95_dp,3,32,marginal_quantile,bound_worst_var, &
    tolerance=0.0_dp,n_lookback=3,max_iterations=300, &
    sample_columns=.false.)
  call assert_true(ra%ok)
  call assert_true(ra%bounds(2)>=ra%bounds(1))
  call assert_true(ra%n_used==32)

  ra = adaptive_ra_bounds(0.95_dp,3,[4,5],marginal_quantile, &
    bound_worst_var,joint_tolerance=1.0_dp,individual_tolerance=0.0_dp, &
    n_lookback=3,max_iterations=300,sample_columns=.false.)
  call assert_true(ra%ok)
  call assert_true(ra%n_used==16 .or. ra%n_used==32)
  call assert_close(normal_cdf(0.0_dp),0.5_dp,1.0e-14_dp)

  print '(a)', 'test_bounds_allocation: PASS'

contains

  subroutine assert_close(actual,expected,tolerance)
    real(dp), intent(in) :: actual,expected,tolerance
    if(abs(actual-expected)>tolerance*max(1.0_dp,abs(expected))) then
      print *, 'mismatch:',actual,expected,abs(actual-expected)
      error stop 1
    end if
  end subroutine assert_close

  subroutine assert_array_close(actual,expected,tolerance)
    real(dp), intent(in) :: actual(:),expected(:),tolerance
    if(maxval(abs(actual-expected))>tolerance*max(1.0_dp,maxval(abs(expected)))) then
      error stop 1
    end if
  end subroutine assert_array_close

  subroutine assert_true(condition)
    logical, intent(in) :: condition
    if(.not.condition) error stop 1
  end subroutine assert_true

end program test_bounds_allocation
