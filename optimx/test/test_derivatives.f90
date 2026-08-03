! SPDX-License-Identifier: GPL-2.0-only
program test_derivatives
  use optimx_mod
  use optimx_example_functions, only: rosenbrock_callback
  implicit none
  type(optimx_problem) :: problem
  type(derivative_check) :: gc
  type(hessian_check) :: hc
  type(kkt_result) :: kc
  real(dp) :: x(2), g(2), h(2,2), value
  integer :: status
  logical :: ok

  call initialize_problem(problem, 2)
  problem%objective => rosenbrock_callback
  problem%has_gradient = .true.
  problem%has_hessian = .true.
  x = [0.8_dp, 0.7_dp]

  call fnchk(problem, x, value, ok, status)
  call check(ok .and. value > 0.0_dp, 'fnchk')
  call grchk(problem, x, gc, 1.0e-4_dp)
  call check(gc%ok, 'grchk')
  call hesschk(problem, x, hc, 2.0e-3_dp)
  call check(hc%ok, 'hesschk')
  call gHgen(problem, x, g, h, status)
  call check(status == 0 .and. maxval(abs(g-gc%analytic)) < 1.0e-10_dp, 'gHgen gradient')

  call kktchk(problem, [1.0_dp,1.0_dp], kc)
  call check(kc%kkt1 .and. kc%kkt2, 'KKT at Rosenbrock optimum')

  call check_pd()
  write(*,'(a)') 'test_derivatives: PASS'
contains
  subroutine check_pd()
    real(dp) :: a(2,2), bound
    logical :: positive
    a = reshape([2.0_dp,0.2_dp,0.2_dp,1.0_dp],[2,2])
    call pd_check(a,positive,bound)
    call check(positive .and. bound > 0.0_dp, 'pd_check positive')
    a(2,2) = -1.0_dp
    call pd_check(a,positive,bound)
    call check(.not. positive, 'pd_check indefinite')
  end subroutine check_pd
  subroutine check(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not. condition) then
      write(*,'(a,a)') 'FAIL: ', label
      error stop 1
    end if
  end subroutine check
end program test_derivatives
