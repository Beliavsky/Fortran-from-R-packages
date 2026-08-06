! SPDX-License-Identifier: GPL-3.0-or-later
module robstattm_logistic
  use robstattm_kinds, only : dp
  use robstattm_types, only : logistic_result, robstattm_success, robstattm_invalid_argument, &
    robstattm_singular, robstattm_no_convergence
  use robstattm_utils, only : add_intercept, safe_logistic, logistic_deviance_residuals, matrix_inverse
  use robustbase_bylogreg, only : by_logistic_result, by_logistic_fit
  use robustbase_linalg, only : least_squares
  use rrcov_types, only : rrcov_covariance_result => covariance_result
  use rrcov_robust, only : cov_mcd
  use rrcov_stats, only : median, mad_scale, chi_square_quantile, normal_cdf
  implicit none
  private
  public :: logreg_by, logreg_wby, logreg_wml
contains
  subroutine logreg_by(x, y, result, intercept, const, max_iter, tol)
    real(dp), intent(in) :: x(:, :), y(:)
    type(logistic_result), intent(out) :: result
    logical, intent(in), optional :: intercept
    real(dp), intent(in), optional :: const, tol
    integer, intent(in), optional :: max_iter
    real(dp), allocatable :: design(:, :)
    type(by_logistic_result) :: fit
    logical :: use_intercept
    real(dp) :: c, tt
    integer :: mi
    use_intercept = .true.
    if (present(intercept)) use_intercept = intercept
    c = 0.5_dp
    if (present(const)) c = const
    mi = 1000
    if (present(max_iter)) mi = max_iter
    tt = 1.0e-7_dp
    if (present(tol)) tt = tol
    call add_intercept(x, use_intercept, design)
    call by_logistic_fit(design, y, fit, const=c, max_iter=mi, tol=tt)
    call copy_by_result(fit, y, result)
    allocate(result%leverage_weights(size(y)))
    result%leverage_weights = 1.0_dp
    result%method = 'BY'
  end subroutine logreg_by

  subroutine logreg_wby(x, y, result, intercept, const, max_iter, tol)
    real(dp), intent(in) :: x(:, :), y(:)
    type(logistic_result), intent(out) :: result
    logical, intent(in), optional :: intercept
    real(dp), intent(in), optional :: const, tol
    integer, intent(in), optional :: max_iter
    real(dp), allocatable :: design(:, :), xfit(:, :), yfit(:)
    logical, allocatable :: keep(:)
    type(by_logistic_result) :: fit
    logical :: use_intercept
    real(dp) :: c, tt
    integer :: mi
    use_intercept = .true.
    if (present(intercept)) use_intercept = intercept
    c = 0.5_dp
    if (present(const)) c = const
    mi = 1000
    if (present(max_iter)) mi = max_iter
    tt = 1.0e-7_dp
    if (present(tol)) tt = tol
    call add_intercept(x, use_intercept, design)
    call leverage_filter(x, keep)
    call pack_rows(design, y, keep, xfit, yfit)
    if (size(yfit) < size(design, 2) + 1) then
      xfit = design
      yfit = y
      keep = .true.
    end if
    call by_logistic_fit(xfit, yfit, fit, const=c, max_iter=mi, tol=tt)
    call copy_by_result_full(fit, design, y, result)
    allocate(result%leverage_weights(size(y)))
    result%leverage_weights = merge(1.0_dp, 0.0_dp, keep)
    result%method = 'WBY'
  end subroutine logreg_wby

  subroutine logreg_wml(x, y, result, intercept, max_iter, tol)
    real(dp), intent(in) :: x(:, :), y(:)
    type(logistic_result), intent(out) :: result
    logical, intent(in), optional :: intercept
    integer, intent(in), optional :: max_iter
    real(dp), intent(in), optional :: tol
    real(dp), allocatable :: design(:, :), xfit(:, :), yfit(:), beta(:), covariance(:, :)
    logical, allocatable :: keep(:)
    logical :: use_intercept, converged
    integer :: mi, iterations, status
    real(dp) :: tt
    use_intercept = .true.
    if (present(intercept)) use_intercept = intercept
    mi = 100
    if (present(max_iter)) mi = max_iter
    tt = 1.0e-8_dp
    if (present(tol)) tt = tol
    call add_intercept(x, use_intercept, design)
    call leverage_filter(x, keep)
    call pack_rows(design, y, keep, xfit, yfit)
    if (size(yfit) < size(design, 2) + 1) then
      xfit = design
      yfit = y
      keep = .true.
    end if
    call logistic_mle(xfit, yfit, beta, covariance, iterations, converged, status, mi, tt)
    if (status /= 0) then
      result%status = robstattm_singular
      return
    end if
    call fill_logistic_result(beta, covariance, design, y, iterations, converged, result)
    allocate(result%leverage_weights(size(y)))
    result%leverage_weights = merge(1.0_dp, 0.0_dp, keep)
    result%method = 'WML'
  end subroutine logreg_wml

  subroutine logistic_mle(x, y, beta, covariance, iterations, converged, status, max_iter, tol)
    real(dp), intent(in) :: x(:, :), y(:)
    real(dp), allocatable, intent(out) :: beta(:), covariance(:, :)
    integer, intent(out) :: iterations, status
    logical, intent(out) :: converged
    integer, intent(in) :: max_iter
    real(dp), intent(in) :: tol
    real(dp), allocatable :: eta(:), fitted(:), weights(:), working(:), xw(:, :), yw(:), new_beta(:), info_matrix(:, :)
    real(dp) :: delta
    integer :: n, p, j, info
    n = size(x, 1)
    p = size(x, 2)
    if (size(y) /= n .or. any(y < 0.0_dp) .or. any(y > 1.0_dp)) then
      allocate(beta(0), covariance(0,0))
      status = 1
      converged = .false.
      iterations = 0
      return
    end if
    allocate(beta(p), new_beta(p), eta(n), fitted(n), weights(n), working(n), xw(n,p), yw(n))
    beta = 0.0_dp
    converged = .false.
    status = 0
    do iterations = 1, max_iter
      eta = matmul(x, beta)
      fitted = safe_logistic(eta)
      weights = max(fitted * (1.0_dp - fitted), 1.0e-8_dp)
      working = eta + (y - fitted) / weights
      do j = 1, p
        xw(:, j) = x(:, j) * sqrt(weights)
      end do
      yw = working * sqrt(weights)
      call least_squares(xw, yw, new_beta, info)
      if (info /= 0) then
        status = info
        exit
      end if
      delta = maxval(abs(new_beta - beta))
      beta = new_beta
      if (delta <= tol * (1.0_dp + maxval(abs(beta)))) then
        converged = .true.
        exit
      end if
    end do
    eta = matmul(x, beta)
    fitted = safe_logistic(eta)
    weights = max(fitted * (1.0_dp - fitted), 1.0e-8_dp)
    info_matrix = matmul(transpose(x), x * spread(weights, 2, p))
    covariance = matrix_inverse(info_matrix, info)
    if (info /= 0) status = info
  end subroutine logistic_mle

  subroutine leverage_filter(x, keep)
    real(dp), intent(in) :: x(:, :)
    logical, allocatable, intent(out) :: keep(:)
    integer, allocatable :: continuous(:)
    real(dp), allocatable :: xc(:, :), distances(:)
    type(rrcov_covariance_result) :: mcd
    real(dp) :: cutoff, scale
    integer :: j, n, p, pc
    logical :: binary
    n = size(x, 1)
    p = size(x, 2)
    allocate(continuous(p))
    pc = 0
    do j = 1, p
      binary = all(abs(x(:,j)) <= 1.0e-12_dp .or. abs(x(:,j)-1.0_dp) <= 1.0e-12_dp)
      if (.not. binary) then
        pc = pc + 1
        continuous(pc) = j
      end if
    end do
    allocate(keep(n))
    keep = .true.
    if (pc == 0) return
    allocate(xc(n, pc))
    do j = 1, pc
      xc(:,j) = x(:,continuous(j))
    end do
    if (pc == 1) then
      scale = mad_scale(xc(:,1))
      if (scale > tiny(1.0_dp)) then
        cutoff = inverse_normal_cdf(0.9875_dp)
        keep = abs(xc(:,1) - median(xc(:,1))) / scale <= cutoff
      end if
    else
      call cov_mcd(xc, mcd, alpha=0.75_dp, nsamp=min(500, max(100, n*5)), reweight=.true.)
      distances = mcd%distances
      cutoff = chi_square_quantile(0.975_dp, real(pc, dp))
      keep = distances <= cutoff
    end if
  end subroutine leverage_filter

  function inverse_normal_cdf(probability) result(value)
    real(dp), intent(in) :: probability
    real(dp) :: value, lo, hi, mid
    integer :: i
    lo = -10.0_dp
    hi = 10.0_dp
    do i = 1, 100
      mid = 0.5_dp * (lo + hi)
      if (normal_cdf(mid) < probability) then
        lo = mid
      else
        hi = mid
      end if
    end do
    value = 0.5_dp * (lo + hi)
  end function inverse_normal_cdf

  subroutine pack_rows(x, y, mask, xout, yout)
    real(dp), intent(in) :: x(:, :), y(:)
    logical, intent(in) :: mask(:)
    real(dp), allocatable, intent(out) :: xout(:, :), yout(:)
    integer :: j
    allocate(xout(count(mask), size(x,2)), yout(count(mask)))
    do j = 1, size(x,2)
      xout(:,j) = pack(x(:,j), mask)
    end do
    yout = pack(y, mask)
  end subroutine pack_rows

  subroutine copy_by_result(source, y, target)
    type(by_logistic_result), intent(in) :: source
    real(dp), intent(in) :: y(:)
    type(logistic_result), intent(out) :: target
    target%coefficients = source%coefficients
    target%fitted = source%fitted
    target%covariance = source%covariance
    target%standard_errors = source%standard_errors
    call logistic_deviance_residuals(y, target%fitted, target%residual_deviances)
    target%objective = source%objective
    target%iterations = source%iterations
    target%converged = source%converged
    target%status = merge(robstattm_success, robstattm_no_convergence, source%converged)
  end subroutine copy_by_result

  subroutine copy_by_result_full(source, design, y, target)
    type(by_logistic_result), intent(in) :: source
    real(dp), intent(in) :: design(:, :), y(:)
    type(logistic_result), intent(out) :: target
    target%coefficients = source%coefficients
    target%fitted = safe_logistic(matmul(design, source%coefficients))
    target%covariance = source%covariance
    target%standard_errors = source%standard_errors
    call logistic_deviance_residuals(y, target%fitted, target%residual_deviances)
    target%objective = source%objective
    target%iterations = source%iterations
    target%converged = source%converged
    target%status = merge(robstattm_success, robstattm_no_convergence, source%converged)
  end subroutine copy_by_result_full

  subroutine fill_logistic_result(beta, covariance, design, y, iterations, converged, result)
    real(dp), intent(in) :: beta(:), covariance(:, :), design(:, :), y(:)
    integer, intent(in) :: iterations
    logical, intent(in) :: converged
    type(logistic_result), intent(out) :: result
    integer :: i
    result%coefficients = beta
    result%covariance = covariance
    allocate(result%standard_errors(size(beta)))
    do i = 1, size(beta)
      result%standard_errors(i) = sqrt(max(covariance(i,i), 0.0_dp))
    end do
    result%fitted = safe_logistic(matmul(design, beta))
    call logistic_deviance_residuals(y, result%fitted, result%residual_deviances)
    result%objective = -sum(y * log(max(result%fitted,1.0e-14_dp)) + &
      (1.0_dp-y)*log(max(1.0_dp-result%fitted,1.0e-14_dp)))
    result%iterations = iterations
    result%converged = converged
    result%status = merge(robstattm_success, robstattm_no_convergence, converged)
  end subroutine fill_logistic_result
end module robstattm_logistic
