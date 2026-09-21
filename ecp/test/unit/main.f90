program test_ecp_unit
  use ecp_api, only: dp, get_within, get_between, e_agglo, agglo_result
  implicit none

  real(dp) :: x(2, 1), y(2, 1), z(6, 1)
  integer :: member(6)
  type(agglo_result) :: agg

  x(:, 1) = [0.0_dp, 2.0_dp]
  y(:, 1) = [4.0_dp, 6.0_dp]
  call assert_close(get_within(1.0_dp, x), 1.0_dp, 1.0e-14_dp, 'getWithin')
  call assert_close(get_between(1.0_dp, x, y), 8.0_dp, 1.0e-14_dp, 'getBetween')

  z(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp, 10.0_dp, 10.0_dp, 10.0_dp]
  member = [1, 2, 3, 4, 5, 6]
  agg = e_agglo(z, member=member, alpha=1.0_dp)
  call assert_true(size(agg%fit) == 6, 'e.agglo fit length')
  call assert_true(size(agg%merged, 1) == 5, 'e.agglo merge count')
  call assert_true(all(agg%progression(1, :) == [1, 2, 3, 4, 5, 6, 7]), 'e.agglo initial boundaries')
  call assert_true(size(agg%cluster) == 6, 'e.agglo cluster length')

  print '(a)', 'All ecp unit tests passed.'

contains

  subroutine assert_close(actual, expected, tol, label)
    real(dp), intent(in) :: actual !! Value produced by the translated routine.
    real(dp), intent(in) :: expected !! Deterministic reference value.
    real(dp), intent(in) :: tol !! Maximum allowed absolute difference.
    character(len=*), intent(in) :: label !! Human-readable test label printed on failure.

    if (abs(actual - expected) > tol) then
      write (*, '(a,2es24.14)') 'FAILED '//trim(label)//': ', actual, expected
      error stop 1
    end if
  end subroutine assert_close

  subroutine assert_true(condition, label)
    logical, intent(in) :: condition !! Logical assertion expected to be true.
    character(len=*), intent(in) :: label !! Human-readable test label printed on failure.

    if (.not. condition) then
      write (*, '(a)') 'FAILED '//trim(label)
      error stop 1
    end if
  end subroutine assert_true

end program test_ecp_unit
