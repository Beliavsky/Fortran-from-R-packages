! SPDX-License-Identifier: GPL-2.0-or-later
module apt_regression
  use apt_kinds, only : dp
  use apt_special, only : student_t_cdf, f_cdf, chi_square_cdf, normal_cdf
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  implicit none
  private

  type, public :: regression_result
    integer :: nobs = 0
    integer :: ncoef = 0
    integer :: df_resid = 0
    integer :: status = 0
    real(dp) :: sse = 0.0_dp
    real(dp) :: sigma2 = 0.0_dp
    real(dp) :: r_squared = 0.0_dp
    real(dp) :: adjusted_r_squared = 0.0_dp
    real(dp) :: f_statistic = 0.0_dp
    real(dp) :: f_p_value = 1.0_dp
    real(dp) :: log_likelihood = 0.0_dp
    real(dp) :: aic = 0.0_dp
    real(dp) :: bic = 0.0_dp
    real(dp), allocatable :: coefficients(:)
    real(dp), allocatable :: standard_errors(:)
    real(dp), allocatable :: t_statistics(:)
    real(dp), allocatable :: p_values(:)
    real(dp), allocatable :: fitted(:)
    real(dp), allocatable :: residuals(:)
    real(dp), allocatable :: covariance(:,:)
  end type regression_result

  type, public :: hypothesis_result
    integer :: restrictions = 0
    integer :: df_denominator = 0
    integer :: status = 0
    real(dp) :: f_statistic = 0.0_dp
    real(dp) :: p_value = 1.0_dp
  end type hypothesis_result

  type, public :: residual_diagnostics
    real(dp) :: durbin_watson = 0.0_dp
    real(dp) :: durbin_watson_p = 1.0_dp
    real(dp) :: ljung_box_4_p = 1.0_dp
    real(dp) :: ljung_box_8_p = 1.0_dp
    real(dp) :: ljung_box_12_p = 1.0_dp
  end type residual_diagnostics

  public :: fit_ols, linear_f_test, zero_coefficient_f_test
  public :: ljung_box_test, durbin_watson_test, compute_residual_diagnostics

  interface
    subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
      import dp
      integer, intent(in) :: n, nrhs, lda, ldb
      integer, intent(out) :: ipiv(*)
      real(dp), intent(inout) :: a(lda, *), b(ldb, *)
      integer, intent(out) :: info
    end subroutine dgesv
    subroutine dgetrf(m, n, a, lda, ipiv, info)
      import dp
      integer, intent(in) :: m, n, lda
      real(dp), intent(inout) :: a(lda, *)
      integer, intent(out) :: ipiv(*)
      integer, intent(out) :: info
    end subroutine dgetrf
    subroutine dgetri(n, a, lda, ipiv, work, lwork, info)
      import dp
      integer, intent(in) :: n, lda, lwork
      real(dp), intent(inout) :: a(lda, *)
      integer, intent(in) :: ipiv(*)
      real(dp), intent(inout) :: work(*)
      integer, intent(out) :: info
    end subroutine dgetri
  end interface

contains

  subroutine fit_ols(y, x, fit)
    real(dp), intent(in) :: y(:)
    real(dp), intent(in) :: x(:,:)
    type(regression_result), intent(out) :: fit
    integer :: n, p, i, info, lwork
    integer, allocatable :: ipiv(:)
    real(dp), allocatable :: xtx(:,:), xty(:,:), work(:)
    real(dp) :: ymean, tss, pi, nanv, df_model, q, mse_ml
    logical :: has_intercept

    nanv = ieee_value(0.0_dp, ieee_quiet_nan)
    n = size(y)
    p = size(x, 2)
    fit%nobs = n
    fit%ncoef = p
    fit%df_resid = n - p
    allocate(fit%coefficients(p), fit%standard_errors(p), fit%t_statistics(p), &
      fit%p_values(p), fit%fitted(n), fit%residuals(n), fit%covariance(p,p))
    fit%coefficients = nanv
    fit%standard_errors = nanv
    fit%t_statistics = nanv
    fit%p_values = nanv
    fit%fitted = nanv
    fit%residuals = nanv
    fit%covariance = nanv

    if (size(x,1) /= n .or. n <= p .or. p < 1) then
      fit%status = 1
      return
    end if

    allocate(xtx(p,p), xty(p,1), ipiv(p))
    xtx = matmul(transpose(x), x)
    xty(:,1) = matmul(transpose(x), y)
    call dgesv(p, 1, xtx, p, ipiv, xty, p, info)
    if (info /= 0) then
      fit%status = 2
      return
    end if
    fit%coefficients = xty(:,1)
    fit%fitted = matmul(x, fit%coefficients)
    fit%residuals = y - fit%fitted
    fit%sse = dot_product(fit%residuals, fit%residuals)
    fit%sigma2 = fit%sse / real(fit%df_resid, dp)

    xtx = matmul(transpose(x), x)
    call dgetrf(p, p, xtx, p, ipiv, info)
    if (info /= 0) then
      fit%status = 3
      return
    end if
    lwork = max(1, 64 * p)
    allocate(work(lwork))
    call dgetri(p, xtx, p, ipiv, work, lwork, info)
    if (info /= 0) then
      fit%status = 4
      return
    end if
    fit%covariance = fit%sigma2 * xtx
    do i = 1, p
      fit%standard_errors(i) = sqrt(max(0.0_dp, fit%covariance(i,i)))
      if (fit%standard_errors(i) > 0.0_dp) then
        fit%t_statistics(i) = fit%coefficients(i) / fit%standard_errors(i)
        fit%p_values(i) = 2.0_dp * (1.0_dp - student_t_cdf(abs(fit%t_statistics(i)), &
          real(fit%df_resid, dp)))
      end if
    end do

    has_intercept = all(abs(x(:,1) - 1.0_dp) <= 100.0_dp*epsilon(1.0_dp))
    if (has_intercept) then
      ymean = sum(y) / real(n, dp)
      tss = sum((y - ymean)**2)
      df_model = real(max(0, p - 1), dp)
      if (tss > 0.0_dp) then
        fit%r_squared = 1.0_dp - fit%sse / tss
        fit%adjusted_r_squared = 1.0_dp - (1.0_dp - fit%r_squared) * &
          real(n - 1, dp) / real(fit%df_resid, dp)
      end if
    else
      tss = dot_product(y, y)
      df_model = real(p, dp)
      if (tss > 0.0_dp) then
        fit%r_squared = 1.0_dp - fit%sse / tss
        fit%adjusted_r_squared = 1.0_dp - (1.0_dp - fit%r_squared) * &
          real(n, dp) / real(fit%df_resid, dp)
      end if
    end if

    if (df_model > 0.0_dp .and. fit%sse > 0.0_dp) then
      q = max(0.0_dp, (tss - fit%sse) / df_model) / fit%sigma2
      fit%f_statistic = q
      fit%f_p_value = 1.0_dp - f_cdf(q, df_model, real(fit%df_resid, dp))
    end if

    pi = acos(-1.0_dp)
    mse_ml = max(fit%sse / real(n,dp), tiny(1.0_dp))
    fit%log_likelihood = -0.5_dp * real(n,dp) * &
      (log(2.0_dp*pi) + 1.0_dp + log(mse_ml))
    fit%aic = -2.0_dp * fit%log_likelihood + 2.0_dp * real(p + 1, dp)
    fit%bic = -2.0_dp * fit%log_likelihood + log(real(n,dp)) * real(p + 1,dp)
    fit%status = 0
  end subroutine fit_ols

  subroutine linear_f_test(fit, rmat, rhs, test)
    type(regression_result), intent(in) :: fit
    real(dp), intent(in) :: rmat(:,:)
    real(dp), intent(in) :: rhs(:)
    type(hypothesis_result), intent(out) :: test
    integer :: m, p, info
    integer, allocatable :: ipiv(:)
    real(dp), allocatable :: s(:,:), d(:,:), sol(:,:)
    real(dp) :: qform

    m = size(rmat,1)
    p = size(rmat,2)
    test%restrictions = m
    test%df_denominator = fit%df_resid
    if (fit%status /= 0 .or. p /= fit%ncoef .or. size(rhs) /= m .or. m < 1) then
      test%status = 1
      return
    end if
    allocate(s(m,m), d(m,1), sol(m,1), ipiv(m))
    d(:,1) = matmul(rmat, fit%coefficients) - rhs
    s = matmul(matmul(rmat, fit%covariance), transpose(rmat))
    sol = d
    call dgesv(m, 1, s, m, ipiv, sol, m, info)
    if (info /= 0) then
      test%status = 2
      return
    end if
    qform = dot_product(d(:,1), sol(:,1))
    test%f_statistic = max(0.0_dp, qform / real(m,dp))
    test%p_value = 1.0_dp - f_cdf(test%f_statistic, real(m,dp), &
      real(fit%df_resid,dp))
    test%status = 0
  end subroutine linear_f_test

  subroutine zero_coefficient_f_test(fit, indices, test)
    type(regression_result), intent(in) :: fit
    integer, intent(in) :: indices(:)
    type(hypothesis_result), intent(out) :: test
    real(dp), allocatable :: r(:,:), rhs(:)
    integer :: i
    allocate(r(size(indices), fit%ncoef), rhs(size(indices)))
    r = 0.0_dp
    rhs = 0.0_dp
    do i = 1, size(indices)
      if (indices(i) < 1 .or. indices(i) > fit%ncoef) then
        test%status = 1
        return
      end if
      r(i, indices(i)) = 1.0_dp
    end do
    call linear_f_test(fit, r, rhs, test)
  end subroutine zero_coefficient_f_test

  subroutine ljung_box_test(residuals, lag, statistic, p_value, status)
    real(dp), intent(in) :: residuals(:)
    integer, intent(in) :: lag
    real(dp), intent(out) :: statistic, p_value
    integer, intent(out), optional :: status
    integer :: n, k, h
    real(dp) :: meanv, denom, rho, q

    n = size(residuals)
    h = min(lag, n - 1)
    if (h < 1) then
      statistic = 0.0_dp
      p_value = 1.0_dp
      if (present(status)) status = 1
      return
    end if
    meanv = sum(residuals) / real(n,dp)
    denom = sum((residuals - meanv)**2)
    if (denom <= 0.0_dp) then
      statistic = 0.0_dp
      p_value = 1.0_dp
      if (present(status)) status = 2
      return
    end if
    q = 0.0_dp
    do k = 1, h
      rho = dot_product(residuals(1:n-k)-meanv, residuals(1+k:n)-meanv) / denom
      q = q + rho * rho / real(n-k,dp)
    end do
    statistic = real(n*(n+2),dp) * q
    p_value = 1.0_dp - chi_square_cdf(statistic, real(h,dp))
    if (present(status)) status = 0
  end subroutine ljung_box_test

  subroutine durbin_watson_test(residuals, statistic, p_value)
    real(dp), intent(in) :: residuals(:)
    real(dp), intent(out) :: statistic, p_value
    integer :: n
    real(dp) :: den, rho, z

    n = size(residuals)
    den = dot_product(residuals, residuals)
    if (n < 3 .or. den <= 0.0_dp) then
      statistic = 0.0_dp
      p_value = 1.0_dp
      return
    end if
    statistic = sum((residuals(2:n)-residuals(1:n-1))**2) / den
    rho = 1.0_dp - 0.5_dp * statistic
    z = rho * sqrt(real(n,dp))
    p_value = 2.0_dp * min(normal_cdf(z), 1.0_dp-normal_cdf(z))
    p_value = min(1.0_dp, max(0.0_dp, p_value))
  end subroutine durbin_watson_test

  subroutine compute_residual_diagnostics(residuals, diag)
    real(dp), intent(in) :: residuals(:)
    type(residual_diagnostics), intent(out) :: diag
    real(dp) :: q
    call durbin_watson_test(residuals, diag%durbin_watson, diag%durbin_watson_p)
    call ljung_box_test(residuals, 4, q, diag%ljung_box_4_p)
    call ljung_box_test(residuals, 8, q, diag%ljung_box_8_p)
    call ljung_box_test(residuals, 12, q, diag%ljung_box_12_p)
  end subroutine compute_residual_diagnostics

end module apt_regression
