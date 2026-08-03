! SPDX-License-Identifier: GPL-2.0-only
program test_orchestration
  use optimx_mod
  use optimx_example_functions, only: quartic_callback
  implicit none
  type(optimx_problem) :: problem
  type(optimx_control) :: control
  type(optimx_multi_result) :: multi, sequence, starts_result, converted
  type(optimx_result) :: single
  real(dp) :: x0(2), starts(2,4)
  character(len=16) :: methods(3)
  integer :: best

  call initialize_problem(problem,2)
  problem%objective => quartic_callback
  problem%has_gradient = .true.
  problem%has_hessian = .true.
  control = ctrldefault(2)
  control%maxit = 800
  x0 = [0.2_dp,0.0_dp]
  methods = [character(len=16) :: 'Nelder-Mead','Rvmmin','snewton']

  call opm(problem,x0,methods,control,multi)
  call check(multi%best >= 1 .and. multi%best <= 3, 'opm best index')
  call check(multi%runs(multi%best)%value < 1.0e-8_dp, 'opm objective')

  starts(:,1) = [-2.0_dp,-2.0_dp]
  starts(:,2) = [-0.2_dp, 3.0_dp]
  starts(:,3) = [ 0.2_dp,-1.0_dp]
  starts(:,4) = [ 2.0_dp, 4.0_dp]
  call multistart(problem,starts,'Rvmmin',control,starts_result)
  call check(starts_result%runs(starts_result%best)%value < 1.0e-8_dp, 'multistart')

  call polyopt(problem,x0,methods,control,sequence)
  call check(sequence%runs(3)%value < 1.0e-8_dp, 'polyopt')
  best = proptimr(multi)
  call check(best == multi%best, 'proptimr')

  call opm2optimr(multi,best,single)
  call optimr2opm(single,converted)
  call check(abs(converted%runs(1)%value-single%value) < 1.0e-14_dp, 'result conversion')

  write(*,'(a)') 'test_orchestration: PASS'
contains
  subroutine check(condition,label)
    logical,intent(in)::condition
    character(len=*),intent(in)::label
    if(.not.condition)then
      write(*,'(a,a)')'FAIL: ',label
      error stop 1
    end if
  end subroutine check
end program test_orchestration
