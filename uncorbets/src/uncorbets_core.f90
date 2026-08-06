! SPDX-License-Identifier: MIT
module uncorbets_core
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use uncorbets_kinds, only : dp
  use uncorbets_types, only : status_type, torsion_result, effective_bets_result, &
      max_effective_bets_result, set_status, uncorbets_ok, &
      uncorbets_invalid_input, uncorbets_no_convergence
  use uncorbets_linalg, only : symmetric_eigen, symmetric_sqrt, solve_linear, &
      inverse_matrix, frobenius_norm, project_simplex, is_symmetric
  implicit none
  private

  public :: sqrtm, torsion, effective_bets, max_effective_bets
  public :: torsion_pca, torsion_minimum, effective_bets_gradient

contains

  function sqrtm(x) result(out)
    real(dp), intent(in) :: x(:, :)
    type(torsion_result) :: out

    call symmetric_sqrt(x, out%matrix, out%status)
  end function sqrtm

  function torsion(sigma, model, method, max_niter, tolerance) result(out)
    real(dp), intent(in) :: sigma(:, :)
    character(len=*), intent(in), optional :: model
    character(len=*), intent(in), optional :: method
    integer, intent(in), optional :: max_niter
    real(dp), intent(in), optional :: tolerance
    type(torsion_result) :: out

    character(len=:), allocatable :: selected_model, selected_method
    integer :: niter
    real(dp) :: tol

    selected_model = 'minimum-torsion'
    if (present(model)) selected_model = lower_ascii(trim(model))
    selected_method = 'exact'
    if (present(method)) selected_method = lower_ascii(trim(method))
    niter = 10000
    if (present(max_niter)) niter = max_niter
    tol = 1.0e-8_dp
    if (present(tolerance)) tol = tolerance

    select case (selected_model)
    case ('pca')
      out = torsion_pca(sigma)
    case ('minimum-torsion', 'minimum_torsion', 'minimum torsion')
      out = torsion_minimum(sigma, selected_method, niter, tol)
    case default
      call set_status(out%status, uncorbets_invalid_input, &
          'model must be pca or minimum-torsion')
    end select
  end function torsion

  function torsion_pca(sigma) result(out)
    real(dp), intent(in) :: sigma(:, :)
    type(torsion_result) :: out

    real(dp), allocatable :: values(:), vectors(:, :)
    type(status_type) :: status
    integer :: n, j

    if (.not. valid_covariance_shape(sigma, out%status)) return
    call symmetric_eigen(sigma, values, vectors, status, descending=.true.)
    if (.not. status%ok()) then
      out%status = status
      return
    end if
    n = size(sigma, 1)
    do j = 1, n
      if (vectors(1, j) < 0.0_dp) vectors(:, j) = -vectors(:, j)
    end do
    allocate(out%matrix(n, n))
    out%matrix = transpose(vectors)
    out%iterations = 0
    out%final_error = 0.0_dp
    out%converged = .true.
    call set_status(out%status, uncorbets_ok, 'ok')
  end function torsion_pca

  function torsion_minimum(sigma, method, max_niter, tolerance) result(out)
    real(dp), intent(in) :: sigma(:, :)
    character(len=*), intent(in) :: method
    integer, intent(in) :: max_niter
    real(dp), intent(in) :: tolerance
    type(torsion_result) :: out

    real(dp), allocatable :: cmat(:, :), croot(:, :), croot_inv(:, :)
    real(dp), allocatable :: dmat(:, :), umat(:, :), uroot(:, :)
    real(dp), allocatable :: qmat(:, :), pimat(:, :), xmat(:, :), rhs(:, :)
    real(dp), allocatable :: d(:), sigma_sd(:)
    type(status_type) :: status
    real(dp) :: current_error, previous_error, relative_change
    integer :: n, i, iter
    character(len=:), allocatable :: selected_method

    if (.not. valid_covariance_shape(sigma, out%status)) return
    n = size(sigma, 1)
    if (max_niter < 1 .or. tolerance <= 0.0_dp) then
      call set_status(out%status, uncorbets_invalid_input, &
          'max_niter and tolerance must be positive')
      return
    end if
    allocate(sigma_sd(n), cmat(n, n))
    do i = 1, n
      if (sigma(i, i) <= 0.0_dp) then
        call set_status(out%status, uncorbets_invalid_input, &
            'covariance diagonal must be positive')
        return
      end if
      sigma_sd(i) = sqrt(sigma(i, i))
    end do
    do i = 1, n
      cmat(i, :) = sigma(i, :) / (sigma_sd(i) * sigma_sd)
    end do
    cmat = 0.5_dp * (cmat + transpose(cmat))
    call symmetric_sqrt(cmat, croot, status)
    if (.not. status%ok()) then
      out%status = status
      return
    end if

    selected_method = lower_ascii(trim(method))
    select case (selected_method)
    case ('approximate')
      call inverse_matrix(croot, croot_inv, status)
      if (.not. status%ok()) then
        out%status = status
        return
      end if
      allocate(out%matrix(n, n))
      do i = 1, n
        out%matrix(i, :) = sigma_sd(i) * croot_inv(i, :)
      end do
      do i = 1, n
        out%matrix(:, i) = out%matrix(:, i) / sigma_sd(i)
      end do
      out%iterations = 0
      out%final_error = frobenius_norm(croot - identity_matrix(n))
      out%converged = .true.
      call set_status(out%status, uncorbets_ok, 'ok')

    case ('exact')
      allocate(d(n), dmat(n, n), umat(n, n), rhs(n, n))
      d = 1.0_dp
      previous_error = huge(1.0_dp)
      current_error = previous_error
      do iter = 1, max_niter
        dmat = diagonal_matrix(d)
        umat = matmul(matmul(dmat, croot), matmul(croot, dmat))
        umat = 0.5_dp * (umat + transpose(umat))
        call symmetric_sqrt(umat, uroot, status)
        if (.not. status%ok()) then
          out%status = status
          return
        end if
        rhs = matmul(dmat, croot)
        call solve_linear(uroot, rhs, qmat, status)
        if (.not. status%ok()) then
          out%status = status
          return
        end if
        do i = 1, n
          d(i) = dot_product(qmat(i, :), croot(:, i))
        end do
        pimat = matmul(diagonal_matrix(d), qmat)
        current_error = frobenius_norm(croot - pimat)
        if (iter > 1) then
          relative_change = abs(current_error - previous_error) / &
              max(current_error, tiny(1.0_dp)) / real(n, dp)
          if (relative_change <= tolerance) exit
        end if
        previous_error = current_error
      end do

      out%iterations = min(iter, max_niter)
      out%final_error = current_error
      out%converged = iter <= max_niter
      call inverse_matrix(croot, croot_inv, status)
      if (.not. status%ok()) then
        out%status = status
        return
      end if
      xmat = matmul(pimat, croot_inv)
      allocate(out%matrix(n, n))
      do i = 1, n
        out%matrix(i, :) = sigma_sd(i) * xmat(i, :)
      end do
      do i = 1, n
        out%matrix(:, i) = out%matrix(:, i) / sigma_sd(i)
      end do
      if (out%converged) then
        call set_status(out%status, uncorbets_ok, 'ok')
      else
        call set_status(out%status, uncorbets_no_convergence, &
            'minimum-torsion iteration reached max_niter')
      end if

    case default
      call set_status(out%status, uncorbets_invalid_input, &
          'method must be approximate or exact')
    end select
  end function torsion_minimum

  function effective_bets(b, sigma, tmat) result(out)
    real(dp), intent(in) :: b(:)
    real(dp), intent(in) :: sigma(:, :)
    real(dp), intent(in) :: tmat(:, :)
    type(effective_bets_result) :: out

    real(dp), allocatable :: solution(:, :), rhs(:, :), y(:), z(:), sb(:)
    type(status_type) :: status
    real(dp) :: variance, entropy
    integer :: n, i

    n = size(b)
    if (.not. valid_problem_dimensions(b, sigma, tmat, out%status)) return
    allocate(rhs(n, 1), y(n), z(n), sb(n), out%probability(n))
    rhs(:, 1) = b
    call solve_linear(transpose(tmat), rhs, solution, status)
    if (.not. status%ok()) then
      out%status = status
      return
    end if
    y = solution(:, 1)
    sb = matmul(sigma, b)
    z = matmul(tmat, sb)
    variance = dot_product(b, sb)
    if (.not. ieee_is_finite(variance) .or. variance <= 0.0_dp) then
      call set_status(out%status, uncorbets_invalid_input, &
          'portfolio variance must be positive')
      return
    end if
    out%probability = y * z / variance
    entropy = 0.0_dp
    do i = 1, n
      if (out%probability(i) > 1.0e-5_dp) then
        entropy = entropy - out%probability(i) * log(out%probability(i))
      end if
    end do
    out%enb = exp(entropy)
    call set_status(out%status, uncorbets_ok, 'ok')
  end function effective_bets

  subroutine effective_bets_gradient(b, sigma, tmat, enb, gradient, status)
    real(dp), intent(in) :: b(:)
    real(dp), intent(in) :: sigma(:, :)
    real(dp), intent(in) :: tmat(:, :)
    real(dp), intent(out) :: enb
    real(dp), allocatable, intent(out) :: gradient(:)
    type(status_type), intent(out) :: status

    real(dp), allocatable :: tinvt(:, :), y(:), z(:), sb(:), q(:), p(:)
    real(dp), allocatable :: jq(:, :), jp(:, :), dhdp(:), tsigma(:, :)
    type(status_type) :: local_status
    real(dp) :: variance, entropy
    integer :: n, i

    n = size(b)
    if (.not. valid_problem_dimensions(b, sigma, tmat, status)) then
      allocate(gradient(0))
      enb = 0.0_dp
      return
    end if
    call inverse_matrix(transpose(tmat), tinvt, local_status)
    if (.not. local_status%ok()) then
      status = local_status
      allocate(gradient(0))
      enb = 0.0_dp
      return
    end if
    allocate(y(n), z(n), sb(n), q(n), p(n), jq(n, n), jp(n, n), &
        dhdp(n), gradient(n), tsigma(n, n))
    y = matmul(tinvt, b)
    sb = matmul(sigma, b)
    z = matmul(tmat, sb)
    tsigma = matmul(tmat, sigma)
    variance = dot_product(b, sb)
    if (variance <= 0.0_dp) then
      call set_status(status, uncorbets_invalid_input, &
          'portfolio variance must be positive')
      gradient = 0.0_dp
      enb = 0.0_dp
      return
    end if
    q = y * z
    p = q / variance
    jq = 0.0_dp
    do i = 1, n
      jq(i, :) = z(i) * tinvt(i, :) + y(i) * tsigma(i, :)
    end do
    jp = jq / variance
    do i = 1, n
      jp(i, :) = jp(i, :) - q(i) * (2.0_dp * sb) / (variance * variance)
    end do
    entropy = 0.0_dp
    dhdp = 0.0_dp
    do i = 1, n
      if (p(i) > 1.0e-5_dp) then
        entropy = entropy - p(i) * log(p(i))
        dhdp(i) = -(log(p(i)) + 1.0_dp)
      end if
    end do
    enb = exp(entropy)
    gradient = enb * matmul(transpose(jp), dhdp)
    call set_status(status, uncorbets_ok, 'ok')
  end subroutine effective_bets_gradient

  function max_effective_bets(x0, sigma, tmat, tolerance, maxeval, maxiter) result(out)
    real(dp), intent(in), optional :: x0(:)
    real(dp), intent(in) :: sigma(:, :)
    real(dp), intent(in) :: tmat(:, :)
    real(dp), intent(in), optional :: tolerance
    integer, intent(in), optional :: maxeval, maxiter
    type(max_effective_bets_result) :: out

    real(dp), allocatable :: x(:), xnew(:), gradient(:), gradient_new(:)
    real(dp), allocatable :: direction(:), step_vector(:), ydiff(:)
    real(dp) :: tol, value, value_new, step, armijo, sy, ss
    integer :: n, max_evaluations, max_iterations, iter, line_iter
    type(status_type) :: status
    type(effective_bets_result) :: eb_result
    logical :: accepted

    n = size(sigma, 1)
    if (size(sigma, 2) /= n .or. size(tmat, 1) /= n .or. size(tmat, 2) /= n) then
      call set_status(out%status, uncorbets_invalid_input, &
          'sigma and t must be square matrices with equal dimensions')
      return
    end if
    tol = 1.0e-10_dp
    if (present(tolerance)) tol = max(tolerance, 10.0_dp * epsilon(1.0_dp))
    max_evaluations = 5000
    if (present(maxeval)) max_evaluations = maxeval
    max_iterations = 5000
    if (present(maxiter)) max_iterations = maxiter
    if (n < 1 .or. max_evaluations < 1 .or. max_iterations < 1) then
      call set_status(out%status, uncorbets_invalid_input, &
          'problem dimension and iteration limits must be positive')
      return
    end if

    allocate(x(n), xnew(n), direction(n), step_vector(n), ydiff(n))
    if (present(x0)) then
      if (size(x0) /= n .or. .not. all(ieee_is_finite(x0))) then
        call set_status(out%status, uncorbets_invalid_input, &
            'x0 must be a finite vector matching sigma')
        return
      end if
      call project_simplex(x0, x)
    else
      x = 1.0_dp / real(n, dp)
    end if

    call effective_bets_gradient(x, sigma, tmat, value, gradient, status)
    out%objective_evaluations = 1
    out%gradient_evaluations = 1
    if (.not. status%ok()) then
      out%status = status
      return
    end if
    step = 1.0_dp
    armijo = 1.0e-4_dp

    do iter = 1, max_iterations
      call project_simplex(x + gradient, xnew)
      direction = xnew - x
      if (maxval(abs(direction)) <= tol) then
        out%converged = .true.
        exit
      end if

      accepted = .false.
      do line_iter = 1, 50
        call project_simplex(x + step * gradient, xnew)
        step_vector = xnew - x
        eb_result = effective_bets(xnew, sigma, tmat)
        out%objective_evaluations = out%objective_evaluations + 1
        if (.not. eb_result%status%ok()) then
          out%status = eb_result%status
          return
        end if
        value_new = eb_result%enb
        if (value_new >= value + armijo * dot_product(gradient, step_vector)) then
          accepted = .true.
          exit
        end if
        step = 0.5_dp * step
        if (out%objective_evaluations >= max_evaluations) exit
      end do
      if (.not. accepted) exit

      call effective_bets_gradient(xnew, sigma, tmat, value_new, gradient_new, status)
      out%gradient_evaluations = out%gradient_evaluations + 1
      if (.not. status%ok()) then
        out%status = status
        return
      end if
      step_vector = xnew - x
      ydiff = gradient_new - gradient
      sy = dot_product(step_vector, ydiff)
      ss = dot_product(step_vector, step_vector)
      if (abs(sy) > sqrt(epsilon(1.0_dp)) * max(1.0_dp, ss)) then
        step = min(100.0_dp, max(1.0e-6_dp, abs(ss / sy)))
      else
        step = 1.0_dp
      end if
      if (abs(value_new - value) <= tol * max(1.0_dp, abs(value_new)) .and. &
          maxval(abs(step_vector)) <= sqrt(tol)) then
        x = xnew
        gradient = gradient_new
        value = value_new
        out%converged = .true.
        exit
      end if
      x = xnew
      gradient = gradient_new
      value = value_new
      if (out%objective_evaluations >= max_evaluations) exit
    end do

    out%iterations = min(iter, max_iterations)
    allocate(out%weights(n), out%gradient(n), out%lambda_lower(n), &
        out%lambda_upper(n), out%hessian(n, n))
    out%weights = x
    out%enb = value
    out%gradient = -gradient
    call numerical_hessian_negative_enb(x, sigma, tmat, out%hessian)
    call estimate_multipliers(x, out%gradient, out%lambda_equality, &
        out%lambda_lower, out%lambda_upper)
    if (out%converged) then
      call set_status(out%status, uncorbets_ok, 'ok')
    else
      call set_status(out%status, uncorbets_no_convergence, &
          'projected-gradient optimizer reached a stopping limit')
    end if
  end function max_effective_bets

  subroutine numerical_hessian_negative_enb(x, sigma, tmat, hessian)
    real(dp), intent(in) :: x(:), sigma(:, :), tmat(:, :)
    real(dp), intent(out) :: hessian(:, :)
    real(dp), allocatable :: xp(:), xm(:), gp(:), gm(:)
    type(status_type) :: status
    real(dp) :: ep, em, h
    integer :: n, j

    n = size(x)
    allocate(xp(n), xm(n))
    hessian = 0.0_dp
    do j = 1, n
      h = epsilon(1.0_dp)**(1.0_dp / 3.0_dp) * max(1.0_dp, abs(x(j)))
      xp = x
      xm = x
      xp(j) = xp(j) + h
      xm(j) = xm(j) - h
      call effective_bets_gradient(xp, sigma, tmat, ep, gp, status)
      if (.not. status%ok()) cycle
      call effective_bets_gradient(xm, sigma, tmat, em, gm, status)
      if (.not. status%ok()) cycle
      hessian(:, j) = -(gp - gm) / (2.0_dp * h)
    end do
    hessian = 0.5_dp * (hessian + transpose(hessian))
  end subroutine numerical_hessian_negative_enb

  subroutine estimate_multipliers(x, gradient, lambda_eq, lambda_lower, lambda_upper)
    real(dp), intent(in) :: x(:), gradient(:)
    real(dp), intent(out) :: lambda_eq
    real(dp), intent(out) :: lambda_lower(:), lambda_upper(:)
    logical, allocatable :: free(:)
    real(dp) :: stationarity
    integer :: i

    allocate(free(size(x)))
    free = x > 1.0e-7_dp .and. x < 1.0_dp - 1.0e-7_dp
    if (any(free)) then
      lambda_eq = -sum(gradient, mask=free) / real(count(free), dp)
    else
      lambda_eq = -sum(gradient) / real(size(x), dp)
    end if
    lambda_lower = 0.0_dp
    lambda_upper = 0.0_dp
    do i = 1, size(x)
      stationarity = gradient(i) + lambda_eq
      if (x(i) <= 1.0e-7_dp) lambda_lower(i) = max(0.0_dp, stationarity)
      if (x(i) >= 1.0_dp - 1.0e-7_dp) lambda_upper(i) = max(0.0_dp, -stationarity)
    end do
  end subroutine estimate_multipliers

  logical function valid_covariance_shape(sigma, status)
    real(dp), intent(in) :: sigma(:, :)
    type(status_type), intent(out) :: status
    integer :: n

    n = size(sigma, 1)
    valid_covariance_shape = .false.
    if (n < 1 .or. size(sigma, 2) /= n) then
      call set_status(status, uncorbets_invalid_input, &
          'sigma must be a nonempty square matrix')
      return
    end if
    if (.not. all(ieee_is_finite(sigma)) .or. .not. is_symmetric(sigma)) then
      call set_status(status, uncorbets_invalid_input, &
          'sigma must be finite and symmetric')
      return
    end if
    valid_covariance_shape = .true.
    call set_status(status, uncorbets_ok, 'ok')
  end function valid_covariance_shape

  logical function valid_problem_dimensions(b, sigma, tmat, status)
    real(dp), intent(in) :: b(:), sigma(:, :), tmat(:, :)
    type(status_type), intent(out) :: status
    integer :: n

    n = size(b)
    valid_problem_dimensions = .false.
    if (n < 1 .or. size(sigma, 1) /= n .or. size(sigma, 2) /= n .or. &
        size(tmat, 1) /= n .or. size(tmat, 2) /= n) then
      call set_status(status, uncorbets_invalid_input, &
          'b, sigma, and t have incompatible dimensions')
      return
    end if
    if (.not. all(ieee_is_finite(b)) .or. .not. all(ieee_is_finite(sigma)) .or. &
        .not. all(ieee_is_finite(tmat))) then
      call set_status(status, uncorbets_invalid_input, &
          'b, sigma, and t must be finite')
      return
    end if
    valid_problem_dimensions = .true.
    call set_status(status, uncorbets_ok, 'ok')
  end function valid_problem_dimensions

  function diagonal_matrix(d) result(a)
    real(dp), intent(in) :: d(:)
    real(dp) :: a(size(d), size(d))
    integer :: i
    a = 0.0_dp
    do i = 1, size(d)
      a(i, i) = d(i)
    end do
  end function diagonal_matrix

  function identity_matrix(n) result(a)
    integer, intent(in) :: n
    real(dp) :: a(n, n)
    integer :: i
    a = 0.0_dp
    do i = 1, n
      a(i, i) = 1.0_dp
    end do
  end function identity_matrix

  function lower_ascii(text) result(lower)
    character(len=*), intent(in) :: text
    character(len=len(text)) :: lower
    integer :: i, code
    lower = text
    do i = 1, len(text)
      code = iachar(lower(i:i))
      if (code >= iachar('A') .and. code <= iachar('Z')) &
          lower(i:i) = achar(code + iachar('a') - iachar('A'))
    end do
  end function lower_ascii

end module uncorbets_core
