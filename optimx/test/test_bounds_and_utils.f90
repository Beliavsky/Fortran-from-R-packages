! SPDX-License-Identifier: GPL-2.0-only
program test_bounds_and_utils
  use optimx_mod
  use optimx_example_functions, only: quadratic_callback
  implicit none
  type(optimx_problem) :: problem
  type(bounds_result) :: br
  type(scale_result) :: sr
  type(optimx_control) :: control
  real(dp) :: x(2), direction(2), stepmax, xbest(2), fbest
  character(len=24), allocatable :: names(:)
  integer :: status

  call initialize_problem(problem,2)
  problem%objective => quadratic_callback
  problem%has_gradient = .true.
  problem%has_hessian = .true.
  problem%lower = [0.0_dp,-1.0_dp]
  problem%upper = [2.0_dp, 3.0_dp]
  x = [-1.0_dp,4.0_dp]
  call bmchk(problem,x,br,.true.)
  call check(maxval(abs(br%par-[0.0_dp,3.0_dp])) < 1.0e-14_dp, 'bmchk shift')

  x = [1.0_dp,1.0_dp]
  direction = [2.0_dp,-1.0_dp]
  call bmstep(problem,x,direction,stepmax)
  call check(abs(stepmax-0.5_dp) < 1.0e-14_dp, 'bmstep')

  control = ctrldefault(2)
  control%initial_step = 0.5_dp
  call axsearch(problem,x,direction,control,xbest,fbest,status)
  call check(status == 0 .and. fbest <= 1.0_dp, 'axsearch')

  call scalechk([1.0_dp,100.0_dp],problem%lower,problem%upper,sr)
  call check(sr%parameter_ratio >= 100.0_dp, 'scalechk')
  call check(checksolver('Rvmmin') .and. checksolver('Nelder-Mead'), 'checksolver known')
  call check(.not.checksolver('imaginary-solver'), 'checksolver unknown')
  call checkallsolvers(names)
  call check(size(names) >= 10, 'checkallsolvers')

  write(*,'(a)') 'test_bounds_and_utils: PASS'
contains
  subroutine check(condition,label)
    logical,intent(in)::condition
    character(len=*),intent(in)::label
    if(.not.condition)then
      write(*,'(a,a)')'FAIL: ',label
      error stop 1
    end if
  end subroutine check
end program test_bounds_and_utils
