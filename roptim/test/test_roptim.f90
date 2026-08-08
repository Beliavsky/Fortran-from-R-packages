program test_roptim
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_positive_inf
  use roptim_mod, only : dp, roptim_control_t, roptim_result_t, &
       roptim_minimize, roptim_approximate_gradient, &
       method_nelder_mead, method_bfgs, method_cg, method_lbfgsb, method_sann, &
       roptim_success, roptim_invalid_input, roptim_user_stop
  implicit none

  integer :: failures

  type :: target_data_t
    real(dp) :: target(2)
  end type target_data_t

  failures = 0
  call test_bfgs_analytic(failures)
  call test_bfgs_numeric(failures)
  call test_nelder_mead(failures)
  call test_cg_variants(failures)
  call test_lbfgsb(failures)
  call test_lbfgsb_numeric(failures)
  call test_sann(failures)
  call test_scaling_maximization(failures)
  call test_custom_proposal(failures)
  call test_user_data(failures)
  call test_invalid_input(failures)
  call test_monitor_stop(failures)
  call test_derivative(failures)

  if (failures /= 0) then
    write(*, '(a,i0)') 'FAILED tests: ', failures
    error stop 1
  end if
  write(*, '(a)') 'All roptim tests passed.'

contains

  subroutine test_bfgs_analytic(failures)
    integer, intent(inout) :: failures
    real(dp) :: x(2)
    type(roptim_control_t) :: control
    type(roptim_result_t) :: result

    x = [-1.2_dp, 1.0_dp]
    control%max_iterations = 300
    control%reltol = 1.0e-10_dp
    control%compute_hessian = .true.
    call roptim_minimize(x, rosenbrock, result, method_bfgs, &
                         gradient=rosenbrock_gradient, control=control)
    call check(result%success, 'BFGS analytic success', failures)
    call check(maxval(abs(x-[1.0_dp,1.0_dp])) < 2.0e-5_dp, &
               'BFGS analytic parameters', failures)
    call check(result%value < 1.0e-10_dp, 'BFGS analytic value', failures)
    call check(allocated(result%hessian), 'BFGS Hessian allocated', failures)
    if (allocated(result%hessian)) then
      call check(abs(result%hessian(1,1)-802.0_dp) < 2.0e-2_dp, &
                 'BFGS Hessian 11', failures)
      call check(abs(result%hessian(1,2)+400.0_dp) < 2.0e-2_dp, &
                 'BFGS Hessian 12', failures)
    end if
  end subroutine test_bfgs_analytic

  subroutine test_bfgs_numeric(failures)
    integer, intent(inout) :: failures
    real(dp) :: x(2)
    type(roptim_control_t) :: control
    type(roptim_result_t) :: result

    x = [-1.2_dp, 1.0_dp]
    control%max_iterations = 400
    control%reltol = 1.0e-8_dp
    allocate(control%ndeps(2), source=1.0e-5_dp)
    call roptim_minimize(x, rosenbrock, result, method_bfgs, control=control)
    call check(result%success, 'BFGS numeric success', failures)
    call check(maxval(abs(x-[1.0_dp,1.0_dp])) < 2.0e-4_dp, &
               'BFGS numeric parameters', failures)
  end subroutine test_bfgs_numeric

  subroutine test_nelder_mead(failures)
    integer, intent(inout) :: failures
    real(dp) :: x(2)
    type(roptim_control_t) :: control
    type(roptim_result_t) :: result

    x = [-1.2_dp, 1.0_dp]
    control%max_iterations = 2000
    control%reltol = 1.0e-9_dp
    call roptim_minimize(x, rosenbrock, result, method_nelder_mead, &
                         control=control)
    call check(result%success, 'Nelder-Mead success', failures)
    call check(maxval(abs(x-[1.0_dp,1.0_dp])) < 2.0e-4_dp, &
               'Nelder-Mead parameters', failures)
  end subroutine test_nelder_mead

  subroutine test_cg_variants(failures)
    integer, intent(inout) :: failures
    real(dp) :: x(3)
    integer :: cg_type
    type(roptim_control_t) :: control
    type(roptim_result_t) :: result

    do cg_type = 1, 3
      x = [5.0_dp, -3.0_dp, 2.0_dp]
      control = roptim_control_t()
      control%max_iterations = 500
      control%reltol = 1.0e-9_dp
      control%cg_type = cg_type
      call roptim_minimize(x, quadratic, result, method_cg, &
                           gradient=quadratic_gradient, control=control)
      call check(result%success, 'CG success', failures)
      call check(maxval(abs(x-[1.0_dp,-2.0_dp,0.5_dp])) < 2.0e-5_dp, &
                 'CG parameters', failures)
    end do
  end subroutine test_cg_variants

  subroutine test_lbfgsb(failures)
    integer, intent(inout) :: failures
    real(dp) :: x(3), lower(3), upper(3)
    type(roptim_control_t) :: control
    type(roptim_result_t) :: result

    x = [0.0_dp, 0.0_dp, 0.0_dp]
    lower = [-1.0_dp, -1.0_dp, -1.0_dp]
    upper = [1.0_dp, 1.0_dp, 1.0_dp]
    control%max_iterations = 200
    control%pgtol = 1.0e-10_dp
    call roptim_minimize(x, bounded_quadratic, result, method_lbfgsb, &
                         gradient=bounded_quadratic_gradient, lower=lower, &
                         upper=upper, control=control)
    call check(result%success, 'L-BFGS-B success', failures)
    call check(maxval(abs(x-[1.0_dp,-1.0_dp,0.25_dp])) < 1.0e-6_dp, &
               'L-BFGS-B parameters', failures)
  end subroutine test_lbfgsb


  subroutine test_lbfgsb_numeric(failures)
    integer, intent(inout) :: failures
    real(dp) :: x(3)
    type(roptim_control_t) :: control
    type(roptim_result_t) :: result

    x = [0.0_dp, 0.0_dp, 0.0_dp]
    control%max_iterations = 300
    control%pgtol = 1.0e-8_dp
    allocate(control%ndeps(3), source=1.0e-6_dp)
    call roptim_minimize(x, bounded_quadratic, result, method_lbfgsb, &
                         lower=[-1.0_dp], upper=[1.0_dp], control=control)
    call check(result%success, 'L-BFGS-B numeric success', failures)
    call check(maxval(abs(x-[1.0_dp,-1.0_dp,0.25_dp])) < 2.0e-5_dp, &
               'L-BFGS-B numeric parameters', failures)
  end subroutine test_lbfgsb_numeric

  subroutine test_sann(failures)
    integer, intent(inout) :: failures
    real(dp) :: x(1)
    type(roptim_control_t) :: control
    type(roptim_result_t) :: result

    x = [8.0_dp]
    control%max_iterations = 8000
    control%temperature = 20.0_dp
    control%tmax = 20
    control%seed = 12345
    allocate(control%parscale(1), source=5.0_dp)
    call roptim_minimize(x, multimodal, result, method_sann, control=control)
    call check(result%success, 'SANN completion', failures)
    call check(result%value < 0.2_dp, 'SANN finds low basin', failures)
  end subroutine test_sann

  subroutine test_scaling_maximization(failures)
    integer, intent(inout) :: failures
    real(dp) :: x(2)
    type(roptim_control_t) :: control
    type(roptim_result_t) :: result

    x = [20.0_dp, -20.0_dp]
    control%max_iterations = 100
    control%fnscale = -1.0_dp
    allocate(control%parscale(2))
    control%parscale = [10.0_dp, 0.5_dp]
    call roptim_minimize(x, concave, result, method_bfgs, &
                         gradient=concave_gradient, control=control)
    call check(result%success, 'maximization success', failures)
    call check(maxval(abs(x-[2.0_dp,-3.0_dp])) < 1.0e-5_dp, &
               'maximization parameters', failures)
    call check(abs(result%value-5.0_dp) < 1.0e-8_dp, &
               'maximization value', failures)
  end subroutine test_scaling_maximization

  subroutine test_custom_proposal(failures)
    integer, intent(inout) :: failures
    real(dp) :: x(4)
    type(roptim_control_t) :: control
    type(roptim_result_t) :: result

    x = [4.0_dp, 3.0_dp, 2.0_dp, 1.0_dp]
    control%max_iterations = 1000
    control%temperature = 2.0_dp
    control%tmax = 10
    control%seed = 7
    call roptim_minimize(x, permutation_cost, result, method_sann, &
                         control=control, proposal=swap_proposal)
    call check(result%success, 'custom SANN completion', failures)
    call check(result%value <= 2.0_dp, 'custom SANN proposal improves', failures)
  end subroutine test_custom_proposal

  subroutine test_user_data(failures)
    integer, intent(inout) :: failures
    real(dp) :: x(2)
    type(target_data_t) :: data
    type(roptim_control_t) :: control
    type(roptim_result_t) :: result

    data%target = [3.0_dp, -4.0_dp]
    x = 0.0_dp
    control%max_iterations = 100
    call roptim_minimize(x, data_objective, result, method_bfgs, &
                         gradient=data_gradient, control=control, user_data=data)
    call check(result%success, 'user data success', failures)
    call check(maxval(abs(x-data%target)) < 1.0e-7_dp, &
               'user data parameters', failures)
  end subroutine test_user_data

  subroutine test_invalid_input(failures)
    integer, intent(inout) :: failures
    real(dp) :: x(2), lower(2), upper(2)
    type(roptim_result_t) :: result

    x = 0.0_dp
    lower = [1.0_dp, 0.0_dp]
    upper = [0.0_dp, 1.0_dp]
    call roptim_minimize(x, quadratic2, result, method_lbfgsb, &
                         lower=lower, upper=upper)
    call check(result%convergence == roptim_invalid_input, &
               'invalid bounds status', failures)
  end subroutine test_invalid_input


  subroutine test_monitor_stop(failures)
    integer, intent(inout) :: failures
    real(dp) :: x(2)
    type(roptim_control_t) :: control
    type(roptim_result_t) :: result

    x = [-1.2_dp, 1.0_dp]
    control%max_iterations = 100
    call roptim_minimize(x, rosenbrock, result, method_bfgs, &
                         gradient=rosenbrock_gradient, control=control, &
                         monitor=stop_after_two)
    call check(result%convergence == roptim_user_stop, &
               'monitor cancellation', failures)
  end subroutine test_monitor_stop

  subroutine test_derivative(failures)
    integer, intent(inout) :: failures
    real(dp) :: x(2), g(2)

    x = [1.5_dp, -0.5_dp]
    call roptim_approximate_gradient(x, quadratic2, g, ndeps=[1.0e-6_dp,1.0e-6_dp])
    call check(maxval(abs(g-2.0_dp*x)) < 1.0e-7_dp, &
               'approximate gradient', failures)
  end subroutine test_derivative

  function rosenbrock(x, user_data) result(f)
    real(dp), intent(in) :: x(:)
    class(*), intent(inout), optional :: user_data
    real(dp) :: f
    f = 100.0_dp*(x(2)-x(1)*x(1))**2 + (1.0_dp-x(1))**2
  end function rosenbrock

  subroutine rosenbrock_gradient(x, g, user_data)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: g(:)
    class(*), intent(inout), optional :: user_data
    g(1) = -400.0_dp*x(1)*(x(2)-x(1)*x(1))-2.0_dp*(1.0_dp-x(1))
    g(2) = 200.0_dp*(x(2)-x(1)*x(1))
  end subroutine rosenbrock_gradient

  function quadratic(x, user_data) result(f)
    real(dp), intent(in) :: x(:)
    class(*), intent(inout), optional :: user_data
    real(dp) :: f
    f = (x(1)-1.0_dp)**2 + 2.0_dp*(x(2)+2.0_dp)**2 + &
        3.0_dp*(x(3)-0.5_dp)**2
  end function quadratic

  subroutine quadratic_gradient(x, g, user_data)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: g(:)
    class(*), intent(inout), optional :: user_data
    g = [2.0_dp*(x(1)-1.0_dp), 4.0_dp*(x(2)+2.0_dp), &
         6.0_dp*(x(3)-0.5_dp)]
  end subroutine quadratic_gradient

  function bounded_quadratic(x, user_data) result(f)
    real(dp), intent(in) :: x(:)
    class(*), intent(inout), optional :: user_data
    real(dp) :: f
    f = (x(1)-2.0_dp)**2 + (x(2)+3.0_dp)**2 + (x(3)-0.25_dp)**2
  end function bounded_quadratic

  subroutine bounded_quadratic_gradient(x, g, user_data)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: g(:)
    class(*), intent(inout), optional :: user_data
    g = 2.0_dp*(x-[2.0_dp,-3.0_dp,0.25_dp])
  end subroutine bounded_quadratic_gradient

  function multimodal(x, user_data) result(f)
    real(dp), intent(in) :: x(:)
    class(*), intent(inout), optional :: user_data
    real(dp) :: f
    f = 0.03_dp*(x(1)-1.5_dp)**2 + sin(2.5_dp*x(1))**2
  end function multimodal

  function concave(x, user_data) result(f)
    real(dp), intent(in) :: x(:)
    class(*), intent(inout), optional :: user_data
    real(dp) :: f
    f = 5.0_dp - (x(1)-2.0_dp)**2 - 2.0_dp*(x(2)+3.0_dp)**2
  end function concave

  subroutine concave_gradient(x, g, user_data)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: g(:)
    class(*), intent(inout), optional :: user_data
    g = [-2.0_dp*(x(1)-2.0_dp), -4.0_dp*(x(2)+3.0_dp)]
  end subroutine concave_gradient

  function permutation_cost(x, user_data) result(f)
    real(dp), intent(in) :: x(:)
    class(*), intent(inout), optional :: user_data
    real(dp) :: f
    integer :: i
    f = 0.0_dp
    do i = 1, size(x)
      f = f + abs(x(i)-real(i,dp))
    end do
  end function permutation_cost

  subroutine swap_proposal(current, candidate, scale, user_data)
    real(dp), intent(in) :: current(:)
    real(dp), intent(out) :: candidate(:)
    real(dp), intent(in) :: scale
    class(*), intent(inout), optional :: user_data
    real(dp) :: u1, u2, tmp
    integer :: i, j
    candidate = current
    call random_number(u1)
    call random_number(u2)
    i = 1 + int(u1*real(size(current),dp))
    j = 1 + int(u2*real(size(current),dp))
    i = min(size(current), i)
    j = min(size(current), j)
    tmp = candidate(i)
    candidate(i) = candidate(j)
    candidate(j) = tmp
  end subroutine swap_proposal

  function data_objective(x, user_data) result(f)
    real(dp), intent(in) :: x(:)
    class(*), intent(inout), optional :: user_data
    real(dp) :: f
    select type (user_data)
    type is (target_data_t)
      f = sum((x-user_data%target)**2)
    class default
      f = huge(1.0_dp)
    end select
  end function data_objective

  subroutine data_gradient(x, g, user_data)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: g(:)
    class(*), intent(inout), optional :: user_data
    select type (user_data)
    type is (target_data_t)
      g = 2.0_dp*(x-user_data%target)
    class default
      g = huge(1.0_dp)
    end select
  end subroutine data_gradient


  subroutine stop_after_two(x, f, iteration, evaluations, stop, user_data)
    real(dp), intent(in) :: x(:), f
    integer, intent(in) :: iteration, evaluations
    logical, intent(out) :: stop
    class(*), intent(inout), optional :: user_data
    stop = iteration >= 2
  end subroutine stop_after_two

  function quadratic2(x, user_data) result(f)
    real(dp), intent(in) :: x(:)
    class(*), intent(inout), optional :: user_data
    real(dp) :: f
    f = sum(x*x)
  end function quadratic2

  subroutine check(condition, label, failures)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    integer, intent(inout) :: failures
    if (.not. condition) then
      failures = failures + 1
      write(*, '(a)') 'FAIL: '//trim(label)
    end if
  end subroutine check

end program test_roptim
