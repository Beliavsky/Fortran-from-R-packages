! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
module multiatsm_optimization
  use multiatsm_kinds, only : dp
  use multiatsm_linalg, only : eye, spectral_radius, cholesky_lower
  implicit none
  private

  abstract interface
    function scalar_objective(x) result(value)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp) :: value
    end function scalar_objective

    subroutine vector_map(x, value)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp), allocatable, intent(out) :: value(:)
    end subroutine vector_map
  end interface

  type, public :: optimization_result
    real(dp), allocatable :: x(:)
    real(dp) :: objective = huge(1.0_dp)
    integer :: iterations = 0
    integer :: evaluations = 0
    integer :: status = 1
    logical :: converged = .false.
  end type optimization_result

  public :: numerical_gradient, numerical_jacobian
  public :: bfgs_minimize, nelder_mead_minimize
  public :: stabilize_transition, lower_factor_to_psd, psd_to_lower_parameters
  public :: block_diagonal_psd, scale_from_jacobian

contains

  subroutine numerical_gradient(fun, x, gradient, evaluations, relative_step)
    procedure(scalar_objective) :: fun
    real(dp), intent(in) :: x(:)
    real(dp), allocatable, intent(out) :: gradient(:)
    integer, intent(out), optional :: evaluations
    real(dp), intent(in), optional :: relative_step
    real(dp), allocatable :: xp(:), xm(:)
    real(dp) :: h, rel
    integer :: i

    rel = epsilon(1.0_dp)**(1.0_dp / 3.0_dp)
    if (present(relative_step)) rel = max(relative_step, epsilon(1.0_dp))
    allocate(gradient(size(x)), xp(size(x)), xm(size(x)))
    do i = 1, size(x)
      h = rel * max(1.0_dp, abs(x(i)))
      xp = x
      xm = x
      xp(i) = xp(i) + h
      xm(i) = xm(i) - h
      gradient(i) = (fun(xp) - fun(xm)) / (2.0_dp * h)
    end do
    if (present(evaluations)) evaluations = 2 * size(x)
  end subroutine numerical_gradient

  subroutine numerical_jacobian(fun, x, jacobian, info, relative_step)
    procedure(vector_map) :: fun
    real(dp), intent(in) :: x(:)
    real(dp), allocatable, intent(out) :: jacobian(:, :)
    integer, intent(out) :: info
    real(dp), intent(in), optional :: relative_step
    real(dp), allocatable :: xp(:), xm(:), fp(:), fm(:), f0(:)
    real(dp) :: h, rel
    integer :: i, m

    call fun(x, f0)
    m = size(f0)
    if (m < 1) then
      allocate(jacobian(0, size(x)))
      info = -1
      return
    end if
    rel = epsilon(1.0_dp)**(1.0_dp / 3.0_dp)
    if (present(relative_step)) rel = max(relative_step, epsilon(1.0_dp))
    allocate(jacobian(m, size(x)), xp(size(x)), xm(size(x)))
    do i = 1, size(x)
      h = rel * max(1.0_dp, abs(x(i)))
      xp = x
      xm = x
      xp(i) = xp(i) + h
      xm(i) = xm(i) - h
      call fun(xp, fp)
      call fun(xm, fm)
      if (size(fp) /= m .or. size(fm) /= m) then
        jacobian = 0.0_dp
        info = -2
        return
      end if
      jacobian(:, i) = (fp - fm) / (2.0_dp * h)
    end do
    info = 0
  end subroutine numerical_jacobian

  function scale_from_jacobian(jacobian, floor_value) result(scale)
    real(dp), intent(in) :: jacobian(:, :)
    real(dp), intent(in), optional :: floor_value
    real(dp) :: scale(size(jacobian, 2))
    real(dp) :: floor_local
    integer :: j

    floor_local = 1.0e-8_dp
    if (present(floor_value)) floor_local = max(floor_value, tiny(1.0_dp))
    do j = 1, size(jacobian, 2)
      scale(j) = max(sqrt(sum(jacobian(:, j) * jacobian(:, j))), floor_local)
    end do
    scale = 1.0_dp / scale
  end function scale_from_jacobian

  subroutine bfgs_minimize(fun, x0, result, tolerance, max_iterations, gradient_tolerance)
    procedure(scalar_objective) :: fun
    real(dp), intent(in) :: x0(:)
    type(optimization_result), intent(out) :: result
    real(dp), intent(in), optional :: tolerance, gradient_tolerance
    integer, intent(in), optional :: max_iterations
    real(dp), allocatable :: x(:), xnew(:), g(:), gnew(:), direction(:), hessian_inv(:, :)
    real(dp), allocatable :: s(:), y(:), ident(:, :), tmp(:, :)
    real(dp) :: f, fnew, tol, gtol, step, slope, rho, ys
    integer :: n, iter, maxit, nev, ls

    n = size(x0)
    tol = 1.0e-8_dp
    if (present(tolerance)) tol = max(tolerance, epsilon(1.0_dp))
    gtol = sqrt(tol)
    if (present(gradient_tolerance)) gtol = max(gradient_tolerance, epsilon(1.0_dp))
    maxit = 1000
    if (present(max_iterations)) maxit = max(1, max_iterations)

    allocate(x(n), xnew(n), direction(n), hessian_inv(n, n), s(n), y(n), ident(n, n), tmp(n, n))
    x = x0
    ident = eye(n)
    hessian_inv = ident
    f = fun(x)
    result%evaluations = 1
    call numerical_gradient(fun, x, g, nev)
    result%evaluations = result%evaluations + nev

    do iter = 1, maxit
      if (maxval(abs(g)) <= gtol) then
        result%converged = .true.
        result%status = 0
        exit
      end if
      direction = -matmul(hessian_inv, g)
      slope = dot_product(g, direction)
      if (slope >= -epsilon(1.0_dp) * max(1.0_dp, dot_product(g, g))) then
        direction = -g
        hessian_inv = ident
        slope = -dot_product(g, g)
      end if

      step = 1.0_dp
      do ls = 1, 40
        xnew = x + step * direction
        fnew = fun(xnew)
        result%evaluations = result%evaluations + 1
        if (fnew <= f + 1.0e-4_dp * step * slope) exit
        step = 0.5_dp * step
      end do
      if (ls > 40 .or. step <= 16.0_dp * epsilon(1.0_dp)) then
        result%status = 2
        exit
      end if

      call numerical_gradient(fun, xnew, gnew, nev)
      result%evaluations = result%evaluations + nev
      s = xnew - x
      y = gnew - g
      ys = dot_product(y, s)
      if (ys > sqrt(epsilon(1.0_dp)) * sqrt(dot_product(y, y) * dot_product(s, s))) then
        rho = 1.0_dp / ys
        tmp = ident - rho * spread(s, 2, n) * spread(y, 1, n)
        hessian_inv = matmul(tmp, matmul(hessian_inv, transpose(tmp))) + &
          rho * spread(s, 2, n) * spread(s, 1, n)
      else
        hessian_inv = ident
      end if

      if (maxval(abs(s)) <= tol * max(1.0_dp, maxval(abs(x)))) then
        x = xnew
        f = fnew
        g = gnew
        result%converged = .true.
        result%status = 0
        exit
      end if
      x = xnew
      f = fnew
      g = gnew
    end do

    result%iterations = min(iter, maxit)
    allocate(result%x(n))
    result%x = x
    result%objective = f
    if (.not. result%converged .and. result%status == 1 .and. iter > maxit) result%status = 1
  end subroutine bfgs_minimize

  subroutine nelder_mead_minimize(fun, x0, result, tolerance, max_iterations, initial_step)
    procedure(scalar_objective) :: fun
    real(dp), intent(in) :: x0(:)
    type(optimization_result), intent(out) :: result
    real(dp), intent(in), optional :: tolerance, initial_step
    integer, intent(in), optional :: max_iterations
    real(dp), allocatable :: simplex(:, :), values(:), centroid(:), trial(:), expanded(:), contracted(:)
    real(dp) :: tol, delta, reflected_value, expanded_value, contracted_value
    real(dp) :: spread_x, spread_f
    integer :: n, maxit, iter, j, best, worst, second_worst

    n = size(x0)
    tol = 1.0e-8_dp
    if (present(tolerance)) tol = max(tolerance, epsilon(1.0_dp))
    delta = 0.05_dp
    if (present(initial_step)) delta = max(initial_step, sqrt(epsilon(1.0_dp)))
    maxit = 2000
    if (present(max_iterations)) maxit = max(1, max_iterations)
    allocate(simplex(n, n + 1), values(n + 1), centroid(n), trial(n), expanded(n), contracted(n))
    simplex(:, 1) = x0
    do j = 1, n
      simplex(:, j + 1) = x0
      simplex(j, j + 1) = simplex(j, j + 1) + delta * max(1.0_dp, abs(x0(j)))
    end do
    do j = 1, n + 1
      values(j) = fun(simplex(:, j))
    end do
    result%evaluations = n + 1

    do iter = 1, maxit
      best = minloc(values, dim=1)
      worst = maxloc(values, dim=1)
      second_worst = best
      do j = 1, n + 1
        if (j /= worst) then
          if (second_worst == best .or. values(j) > values(second_worst)) second_worst = j
        end if
      end do
      spread_x = 0.0_dp
      do j = 1, n + 1
        spread_x = max(spread_x, maxval(abs(simplex(:, j) - simplex(:, best))))
      end do
      spread_f = maxval(abs(values - values(best)))
      if (spread_x <= tol * max(1.0_dp, maxval(abs(simplex(:, best)))) .and. &
          spread_f <= tol * max(1.0_dp, abs(values(best)))) then
        result%converged = .true.
        result%status = 0
        exit
      end if

      centroid = 0.0_dp
      do j = 1, n + 1
        if (j /= worst) centroid = centroid + simplex(:, j)
      end do
      centroid = centroid / real(n, dp)
      trial = centroid + (centroid - simplex(:, worst))
      reflected_value = fun(trial)
      result%evaluations = result%evaluations + 1

      if (reflected_value < values(best)) then
        expanded = centroid + 2.0_dp * (trial - centroid)
        expanded_value = fun(expanded)
        result%evaluations = result%evaluations + 1
        if (expanded_value < reflected_value) then
          simplex(:, worst) = expanded
          values(worst) = expanded_value
        else
          simplex(:, worst) = trial
          values(worst) = reflected_value
        end if
      else if (reflected_value < values(second_worst)) then
        simplex(:, worst) = trial
        values(worst) = reflected_value
      else
        if (reflected_value < values(worst)) then
          contracted = centroid + 0.5_dp * (trial - centroid)
        else
          contracted = centroid + 0.5_dp * (simplex(:, worst) - centroid)
        end if
        contracted_value = fun(contracted)
        result%evaluations = result%evaluations + 1
        if (contracted_value < min(values(worst), reflected_value)) then
          simplex(:, worst) = contracted
          values(worst) = contracted_value
        else
          do j = 1, n + 1
            if (j /= best) then
              simplex(:, j) = simplex(:, best) + 0.5_dp * (simplex(:, j) - simplex(:, best))
              values(j) = fun(simplex(:, j))
              result%evaluations = result%evaluations + 1
            end if
          end do
        end if
      end if
    end do

    best = minloc(values, dim=1)
    allocate(result%x(n))
    result%x = simplex(:, best)
    result%objective = values(best)
    result%iterations = min(iter, maxit)
  end subroutine nelder_mead_minimize

  subroutine stabilize_transition(matrix, stable, info, upper_bound)
    real(dp), intent(in) :: matrix(:, :)
    real(dp), allocatable, intent(out) :: stable(:, :)
    integer, intent(out) :: info
    real(dp), intent(in), optional :: upper_bound
    real(dp) :: radius, upper

    if (size(matrix, 1) /= size(matrix, 2)) then
      allocate(stable(0, 0))
      info = -1
      return
    end if
    upper = 0.9999_dp
    if (present(upper_bound)) upper = min(max(upper_bound, 1.0e-8_dp), 1.0_dp - epsilon(1.0_dp))
    radius = spectral_radius(matrix, info)
    if (info /= 0) then
      allocate(stable(0, 0))
      return
    end if
    allocate(stable(size(matrix, 1), size(matrix, 2)))
    stable = matrix
    if (radius > upper) stable = matrix * (upper / radius)
    info = 0
  end subroutine stabilize_transition

  subroutine lower_factor_to_psd(parameters, n, covariance, lower, info, positive_diagonal)
    real(dp), intent(in) :: parameters(:)
    integer, intent(in) :: n
    real(dp), allocatable, intent(out) :: covariance(:, :), lower(:, :)
    integer, intent(out) :: info
    logical, intent(in), optional :: positive_diagonal
    logical :: positive
    integer :: i, j, k

    if (n < 1 .or. size(parameters) /= n * (n + 1) / 2) then
      allocate(covariance(0, 0), lower(0, 0))
      info = -1
      return
    end if
    positive = .false.
    if (present(positive_diagonal)) positive = positive_diagonal
    allocate(lower(n, n), covariance(n, n))
    lower = 0.0_dp
    k = 0
    do j = 1, n
      do i = j, n
        k = k + 1
        lower(i, j) = parameters(k)
        if (positive .and. i == j) lower(i, j) = exp(parameters(k))
      end do
    end do
    covariance = matmul(lower, transpose(lower))
    info = 0
  end subroutine lower_factor_to_psd

  subroutine psd_to_lower_parameters(covariance, parameters, info)
    real(dp), intent(in) :: covariance(:, :)
    real(dp), allocatable, intent(out) :: parameters(:)
    integer, intent(out) :: info
    real(dp), allocatable :: lower(:, :)
    integer :: n, i, j, k

    n = size(covariance, 1)
    if (size(covariance, 2) /= n) then
      allocate(parameters(0))
      info = -1
      return
    end if
    call cholesky_lower(covariance, lower, info, 1.0e-12_dp)
    if (info /= 0) then
      allocate(parameters(0))
      return
    end if
    allocate(parameters(n * (n + 1) / 2))
    k = 0
    do j = 1, n
      do i = j, n
        k = k + 1
        parameters(k) = lower(i, j)
      end do
    end do
  end subroutine psd_to_lower_parameters

  subroutine block_diagonal_psd(parameters, block_sizes, covariance, info)
    real(dp), intent(in) :: parameters(:)
    integer, intent(in) :: block_sizes(:)
    real(dp), allocatable, intent(out) :: covariance(:, :)
    integer, intent(out) :: info
    real(dp), allocatable :: block(:, :), lower(:, :)
    integer :: total, required, b, offset, p0, p1, count

    if (any(block_sizes < 0)) then
      allocate(covariance(0, 0))
      info = -1
      return
    end if
    total = sum(block_sizes)
    required = sum(block_sizes * (block_sizes + 1) / 2)
    if (size(parameters) /= required) then
      allocate(covariance(0, 0))
      info = -2
      return
    end if
    allocate(covariance(total, total))
    covariance = 0.0_dp
    offset = 0
    p0 = 1
    do b = 1, size(block_sizes)
      count = block_sizes(b) * (block_sizes(b) + 1) / 2
      p1 = p0 + count - 1
      if (block_sizes(b) > 0) then
        call lower_factor_to_psd(parameters(p0:p1), block_sizes(b), block, lower, info)
        if (info /= 0) return
        covariance(offset + 1:offset + block_sizes(b), offset + 1:offset + block_sizes(b)) = block
      end if
      offset = offset + block_sizes(b)
      p0 = p1 + 1
    end do
    info = 0
  end subroutine block_diagonal_psd

end module multiatsm_optimization
