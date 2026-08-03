! SPDX-License-Identifier: MIT
module mfgarch_optimization
  use mfgarch_kinds, only : dp
  use mfgarch_math, only : identity_matrix, vector_norm, finite_value
  use mfgarch_status, only : mfgarch_success, mfgarch_invalid_argument, &
    mfgarch_not_converged, mfgarch_numerical_error
  implicit none
  private

  abstract interface
    function scalar_objective(x) result(value)
      import :: dp
      real(dp), intent(in) :: x(:)
      real(dp) :: value
    end function scalar_objective
  end interface

  public :: scalar_objective, optimize_bfgs, optimize_nelder_mead, numerical_gradient

contains

  subroutine numerical_gradient(function, x, step, gradient, value, evaluations, status)
    procedure(scalar_objective) :: function
    real(dp), intent(in) :: x(:), step
    real(dp), allocatable, intent(out) :: gradient(:)
    real(dp), intent(out) :: value
    integer, intent(inout) :: evaluations
    integer, intent(out) :: status
    real(dp), allocatable :: xp(:), xm(:)
    real(dp) :: fp, fm, h
    integer :: i

    status = mfgarch_success
    allocate(gradient(size(x)), xp(size(x)), xm(size(x)))
    value = function(x)
    evaluations = evaluations + 1
    if (.not. finite_value(value)) then
      status = mfgarch_numerical_error
      gradient = 0.0_dp
      return
    end if
    do i = 1, size(x)
      xp = x
      xm = x
      h = step * max(1.0_dp, abs(x(i)))
      xp(i) = xp(i) + h
      xm(i) = xm(i) - h
      fp = function(xp)
      fm = function(xm)
      evaluations = evaluations + 2
      if (.not. finite_value(fp) .or. .not. finite_value(fm)) then
        status = mfgarch_numerical_error
        gradient = 0.0_dp
        return
      end if
      gradient(i) = (fp - fm) / (2.0_dp * h)
    end do
  end subroutine numerical_gradient

  subroutine optimize_bfgs(function, start, tolerance, gradient_step, max_iterations, &
      max_evaluations, trace, solution, objective, status, iterations, evaluations)
    procedure(scalar_objective) :: function
    real(dp), intent(in) :: start(:), tolerance, gradient_step
    integer, intent(in) :: max_iterations, max_evaluations
    logical, intent(in) :: trace
    real(dp), allocatable, intent(out) :: solution(:)
    real(dp), intent(out) :: objective
    integer, intent(out) :: status, iterations, evaluations
    real(dp), allocatable :: x(:), x_new(:), gradient(:), gradient_new(:)
    real(dp), allocatable :: inverse_hessian(:,:), direction(:), s(:), y(:), hy(:)
    real(dp) :: value, value_new, alpha, slope, ys, yhy, rho
    integer :: n, i, gradient_status

    n = size(start)
    allocate(solution(n), x(n), x_new(n), inverse_hessian(n,n), direction(n), s(n), y(n), hy(n))
    x = start
    inverse_hessian = identity_matrix(n)
    evaluations = 0
    iterations = 0
    call numerical_gradient(function, x, gradient_step, gradient, value, evaluations, status)
    if (status /= mfgarch_success) then
      solution = x
      objective = value
      return
    end if

    do i = 1, max_iterations
      iterations = i
      if (vector_norm(gradient) <= tolerance * (1.0_dp + abs(value))) then
        status = mfgarch_success
        exit
      end if
      if (evaluations >= max_evaluations) then
        status = mfgarch_not_converged
        exit
      end if

      direction = -matmul(inverse_hessian, gradient)
      slope = dot_product(gradient, direction)
      if (slope >= -sqrt(epsilon(1.0_dp)) * max(1.0_dp, vector_norm(gradient))) then
        inverse_hessian = identity_matrix(n)
        direction = -gradient
        slope = -dot_product(gradient, gradient)
      end if

      alpha = 1.0_dp
      do
        x_new = x + alpha * direction
        value_new = function(x_new)
        evaluations = evaluations + 1
        if (finite_value(value_new)) then
          if (value_new <= value + 1.0e-4_dp * alpha * slope) exit
        end if
        alpha = 0.5_dp * alpha
        if (alpha < 1.0e-12_dp .or. evaluations >= max_evaluations) exit
      end do
      if (alpha < 1.0e-12_dp .or. .not. finite_value(value_new)) then
        status = mfgarch_not_converged
        exit
      end if

      call numerical_gradient(function, x_new, gradient_step, gradient_new, value_new, &
        evaluations, gradient_status)
      if (gradient_status /= mfgarch_success) then
        status = gradient_status
        exit
      end if
      s = x_new - x
      y = gradient_new - gradient
      ys = dot_product(y, s)
      if (ys > 1.0e-12_dp * vector_norm(y) * max(vector_norm(s), tiny(1.0_dp))) then
        hy = matmul(inverse_hessian, y)
        yhy = dot_product(y, hy)
        rho = 1.0_dp / ys
        inverse_hessian = inverse_hessian + ((ys + yhy) * rho * rho) * &
          outer_product_local(s, s) - rho * (outer_product_local(hy, s) + outer_product_local(s, hy))
      else
        inverse_hessian = identity_matrix(n)
      end if
      x = x_new
      value = value_new
      call move_alloc(gradient_new, gradient)
      if (trace) write(*,'(a,i0,2(a,es14.6))') 'BFGS iteration ', i, ' objective=', value, &
        ' gradient norm=', vector_norm(gradient)
      if (maxval(abs(s)) <= tolerance * (1.0_dp + maxval(abs(x)))) then
        status = mfgarch_success
        exit
      end if
      status = mfgarch_not_converged
    end do
    solution = x
    objective = value
  end subroutine optimize_bfgs

  subroutine optimize_nelder_mead(function, start, initial_step, tolerance, max_iterations, &
      max_evaluations, trace, solution, objective, status, iterations, evaluations)
    procedure(scalar_objective) :: function
    real(dp), intent(in) :: start(:), initial_step, tolerance
    integer, intent(in) :: max_iterations, max_evaluations
    logical, intent(in) :: trace
    real(dp), allocatable, intent(out) :: solution(:)
    real(dp), intent(out) :: objective
    integer, intent(out) :: status, iterations, evaluations
    real(dp), allocatable :: simplex(:,:), values(:), centroid(:), reflected(:), expanded(:), contracted(:)
    real(dp) :: fr, fe, fc, scale
    integer :: n, i, j

    n = size(start)
    if (n == 0 .or. initial_step <= 0.0_dp) then
      status = mfgarch_invalid_argument
      iterations = 0
      evaluations = 0
      allocate(solution(0))
      objective = huge(1.0_dp)
      return
    end if
    allocate(simplex(n,n+1), values(n+1), centroid(n), reflected(n), expanded(n), contracted(n), solution(n))
    simplex(:,1) = start
    do j = 1, n
      simplex(:,j+1) = start
      scale = initial_step * max(1.0_dp, abs(start(j)))
      simplex(j,j+1) = simplex(j,j+1) + scale
    end do
    evaluations = 0
    do j = 1, n + 1
      values(j) = function(simplex(:,j))
      evaluations = evaluations + 1
    end do
    status = mfgarch_not_converged
    iterations = 0

    do i = 1, max_iterations
      iterations = i
      call sort_simplex(simplex, values)
      if (maxval(abs(simplex(:,2:n+1) - spread(simplex(:,1),2,n))) <= &
          tolerance * (1.0_dp + maxval(abs(simplex(:,1)))) .and. &
          maxval(abs(values(2:n+1) - values(1))) <= tolerance * (1.0_dp + abs(values(1)))) then
        status = mfgarch_success
        exit
      end if
      if (evaluations >= max_evaluations) exit

      centroid = sum(simplex(:,1:n), dim=2) / real(n, dp)
      reflected = centroid + (centroid - simplex(:,n+1))
      fr = function(reflected)
      evaluations = evaluations + 1
      if (fr < values(1)) then
        expanded = centroid + 2.0_dp * (reflected - centroid)
        fe = function(expanded)
        evaluations = evaluations + 1
        if (fe < fr) then
          simplex(:,n+1) = expanded
          values(n+1) = fe
        else
          simplex(:,n+1) = reflected
          values(n+1) = fr
        end if
      else if (fr < values(n)) then
        simplex(:,n+1) = reflected
        values(n+1) = fr
      else
        if (fr < values(n+1)) then
          contracted = centroid + 0.5_dp * (reflected - centroid)
        else
          contracted = centroid + 0.5_dp * (simplex(:,n+1) - centroid)
        end if
        fc = function(contracted)
        evaluations = evaluations + 1
        if (fc < min(fr, values(n+1))) then
          simplex(:,n+1) = contracted
          values(n+1) = fc
        else
          do j = 2, n + 1
            simplex(:,j) = simplex(:,1) + 0.5_dp * (simplex(:,j) - simplex(:,1))
            values(j) = function(simplex(:,j))
            evaluations = evaluations + 1
          end do
        end if
      end if
      if (trace .and. mod(i,10) == 0) write(*,'(a,i0,a,es14.6)') &
        'Nelder-Mead iteration ', i, ' objective=', minval(values)
    end do
    call sort_simplex(simplex, values)
    solution = simplex(:,1)
    objective = values(1)
  end subroutine optimize_nelder_mead

  pure function outer_product_local(a, b) result(matrix)
    real(dp), intent(in) :: a(:), b(:)
    real(dp) :: matrix(size(a),size(b))
    matrix = spread(a,2,size(b)) * spread(b,1,size(a))
  end function outer_product_local

  subroutine sort_simplex(simplex, values)
    real(dp), intent(inout) :: simplex(:,:), values(:)
    real(dp), allocatable :: column(:)
    real(dp) :: temp
    integer :: i, j, minimum

    allocate(column(size(simplex,1)))
    do i = 1, size(values) - 1
      minimum = i
      do j = i + 1, size(values)
        if (values(j) < values(minimum)) minimum = j
      end do
      if (minimum /= i) then
        temp = values(i)
        values(i) = values(minimum)
        values(minimum) = temp
        column = simplex(:,i)
        simplex(:,i) = simplex(:,minimum)
        simplex(:,minimum) = column
      end if
    end do
  end subroutine sort_simplex

end module mfgarch_optimization
