! SPDX-License-Identifier: MIT
! Copyright (c) 2014 Matthew R. Barry
program test_combinations
  use iso_fortran_env, only : int64
  use pbo, only : binomial_coefficient, generate_combinations
  implicit none
  integer, allocatable :: combos(:,:)
  integer :: expected(2,6)
  logical :: success
  character(len=:), allocatable :: message

  expected(:,1) = [1,2]
  expected(:,2) = [1,3]
  expected(:,3) = [1,4]
  expected(:,4) = [2,3]
  expected(:,5) = [2,4]
  expected(:,6) = [3,4]
  call assert_true(binomial_coefficient(4,2) == 6_int64, '4 choose 2')
  call assert_true(binomial_coefficient(10,0) == 1_int64, '10 choose 0')
  call generate_combinations(4,2,combos,success,message)
  call assert_true(success, message)
  call assert_true(all(shape(combos) == [2,6]), 'combination shape')
  call assert_true(all(combos == expected), 'combination ordering')
  print '(a)', 'test_combinations: PASS'
contains
  subroutine assert_true(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not. condition) then
      print '(a)', 'FAILED: ' // trim(label)
      error stop 1
    end if
  end subroutine assert_true
end program test_combinations
